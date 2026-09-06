# A folder of output for one fit, assembled from the package's existing
# reporting functions.
#
# Design, because it is not obvious from the code:
#
# 1. Every component is one existing, tested function. This file contains no
#    new statistics except .ctReportProfile(), which replaces the dead
#    llsurface() -- see the comment there. If a component needs a new
#    calculation, the calculation belongs in the function that owns it.
#
# 2. Components are numbered, and the number is the filename prefix, so an
#    alphabetical directory listing is the intended reading order. The numbers
#    are cheap to renumber and are not referred to anywhere else.
#
# 3. Components that only some backends support are *skipped with a reason*,
#    never errored. Handing a julia fit to a report that half of the package
#    cannot yet describe should still produce the half it can. Same for a
#    component whose optional dependency (LaTeX) is missing, and for one that
#    throws: the index records what happened and the run continues. A report
#    that aborts two thirds of the way through is worth less than a report with
#    two thirds of its rows marked "failed".
#
# 4. The index is written last, from a record of what actually happened, not
#    from the requested component list.
#
# 5. Interpretation lines are the point of the index. A user who can already
#    read a lagged covariance plot does not need this function; the one who
#    ran ctFit() because a colleague told them to does. Keep each to a phrase
#    -- the reasoning belongs in comments, here and in the component
#    functions.

# The component registry: the run order, and the cheapest set each component
# belongs to. The numeric filename prefix is not repeated here -- it lives in
# the component's own file names in .ctReportDo, so there is one place for it
# to be wrong rather than two that can disagree. Order here is the order they
# run in, and must match the prefixes.
.ctReportComponents <- function() c(
  summary        = "quick",
  parmatrices    = "quick",
  identification = "quick",
  profile        = "quick",
  discretepars   = "quick",
  network        = "quick",
  predictions    = "quick",
  residuals      = "quick",
  covcheck       = "default",
  postpred       = "default",
  tipredeffects  = "default",
  crossval       = "all",
  equations      = "all"
)

.ctReportSets <- function() {
  sets <- .ctReportComponents()
  list(quick = names(sets)[sets == "quick"],
    default = names(sets)[sets %in% c("quick", "default")],
    all = names(sets))
}

.ctReportResolve <- function(components) {
  sets <- .ctReportSets()
  known <- names(.ctReportComponents())
  if (length(components) == 1L && components %in% names(sets)) {
    return(sets[[components]])
  }
  bad <- setdiff(components, c(known, names(sets)))
  if (length(bad)) {
    stop("Unknown report component(s): ", paste(bad, collapse = ", "),
      ". Available: ", paste(c(names(sets), known), collapse = ", "), call. = FALSE)
  }
  out <- unique(unlist(lapply(components,
    function(x) if (x %in% names(sets)) sets[[x]] else x)))
  known[known %in% out]
}


# --- shared helpers -----------------------------------------------------------

# pdf() plus a dev.off() that survives an error in the plotting code. `expr` is
# a promise, so force() runs it in the caller's frame.
.ctReportPdf <- function(path, expr, width = 9, height = 7) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(NULL)
}

.ctReportWrite <- function(path, lines) {
  writeLines(as.character(lines), path)
  invisible(NULL)
}

# print() output as lines, with a max.print high enough that a parameter table
# is not truncated mid-matrix and a width wide enough that one is not wrapped
# into two blocks. Restores the options whatever happens.
.ctReportCapture <- function(expr) {
  mp <- options(max.print = 99999L, width = 200L)
  on.exit(options(mp), add = TRUE)
  utils::capture.output(force(expr))
}

# Two different standard errors per raw parameter, and the report needs both.
#
# `marginal` is sqrt of the covariance diagonal -- the width the summary
# reports, which already accounts for every other parameter being unknown too.
# `conditional` is 1/sqrt(-H_ii) -- the width implied by the curvature along
# that one coordinate with every other parameter held fixed.
#
# They differ by however correlated the parameter is with the rest, and that
# ratio is the ridge diagnostic: marginal/conditional of 5 means this
# parameter is individually well determined and jointly nearly free. Getting
# these the wrong way round is easy and quiet -- a slice walks one coordinate
# with the others fixed, so it probes the *conditional* curvature, and scaling
# its grid by the marginal SE makes a perfectly quadratic likelihood look
# violently non-quadratic.
.ctReportRawSE <- function(fit) {
  cov <- if (inherits(fit, "ctJuliaFit")) fit$uncertainty$cov else fit$stanfit$cov
  hess <- if (inherits(fit, "ctJuliaFit")) fit$uncertainty$hessian else fit$stanfit$uncertainty$hessian

  marginal <- NULL
  if (!is.null(cov) && is.matrix(cov)) marginal <- suppressWarnings(sqrt(diag(cov)))
  if (is.null(marginal) || any(!is.finite(marginal)) || all(marginal == 0)) {
    draws <- .ctFitRawPosterior(fit)
    if (!is.null(draws) && length(dim(draws)) == 2L && nrow(draws) > 2L) {
      marginal <- apply(draws, 2, stats::sd)
    }
  }

  conditional <- NULL
  if (!is.null(hess) && is.matrix(hess) && nrow(hess) == ncol(hess)) {
    # Hessian of the log probability, so negative definite at a maximum.
    curv <- -diag(hess)
    conditional <- ifelse(curv > 0, 1 / sqrt(pmax(curv, .Machine$double.eps)), Inf)
  }

  clean <- function(x) {
    if (is.null(x)) return(NULL)
    # as.numeric() drops the names these pick up from a dimnamed covariance.
    # Left on, they become spurious row names in the profile's data frames.
    x <- as.numeric(x)
    # A zero or non-finite width is a parameter the uncertainty step could say
    # nothing about. Give it the median rather than dropping it -- a slice
    # through exactly that parameter is what the profile is for.
    fallback <- stats::median(x[is.finite(x) & x > 0])
    if (!is.finite(fallback) || fallback <= 0) fallback <- 0.1
    x[!is.finite(x) | x <= 0] <- fallback
    x
  }
  list(marginal = clean(marginal), conditional = clean(conditional),
    conditional_raw = conditional)
}

# A log-probability function of the raw parameter vector, on either backend.
# Both of these already exist and both already return -1e100 rather than
# erroring at an invalid point, which is what makes a slice safe to walk.
.ctReportLpg <- function(fit) {
  if (inherits(fit, "ctJuliaFit")) return(.ctBackendLpgFunc(fit))
  if (is.null(fit$stanmodel) || is.null(fit$standata)) return(NULL)
  ctOptimFitLpgFunc(fit, cores = 1)$lpg
}

# Subjects worth plotting: the ones with the most rows, since a subject with
# two observations shows nothing. Returns original (user) ids.
.ctReportSubjects <- function(fit, n = 6L) {
  idmap <- .ctFitIdMap(fit)
  rows <- .ctFitRowSubject(fit)
  counts <- table(rows)
  internal <- as.integer(names(counts))[order(as.integer(counts), decreasing = TRUE)]
  internal <- internal[seq_len(min(n, length(internal)))]
  out <- idmap[match(internal, idmap[, 2]), 1]
  out[!is.na(out)]
}

# An interval grid for ctDiscretePars: zero to the 90th percentile of the
# per-subject observed span, which is the range over which a discrete-time
# coefficient from this data means anything.
.ctReportTimes <- function(fit, n = 20L) {
  tt <- .ctFitRowTime(fit)
  ss <- .ctFitRowSubject(fit)
  spans <- tapply(tt, ss, function(x) diff(range(x, na.rm = TRUE)))
  span <- suppressWarnings(as.numeric(stats::quantile(spans, .9, na.rm = TRUE)))
  if (!is.finite(span) || span <= 0) span <- 1
  seq(0, span, length.out = n)
}


# --- the profile --------------------------------------------------------------

# What llsurface() was reaching for, made into something a user can act on.
#
# The old version sampled a cloud of points around the estimate and smoothed
# the log probability against each coordinate, which conflates the effect of
# that coordinate with everything else moving at the same time. This walks one
# coordinate at a time, holding the rest at the estimate, over a grid measured
# in that parameter's *conditional* standard error -- the width its own
# curvature implies. That makes the reference curve exact rather than
# eyeballed: a quadratic log probability drops by k^2/2 at k conditional
# standard errors, so the drop at +/- 2 is 2.0 whatever the model. Four
# departures each mean something specific:
#
#   off-peak      the slice maximum is not at the estimate -- the optimizer
#                 did not arrive, whatever it reported
#   non-quadratic the drop is far from 2 -- the surface is not the shape the
#                 standard error assumes
#   asymmetric    the two sides disagree -- so a symmetric interval is the
#                 wrong shape even if its width is right
#   ridged        the reported (marginal) interval is much wider than this
#                 curvature -- the parameter is pinned on its own and nearly
#                 free jointly, i.e. it sits on a ridge with others
#
# `ridged` needs no likelihood evaluation at all; it is in this table because
# it is the thing that explains an otherwise puzzling clean slice under a wide
# reported interval, and the two belong side by side.
#
# It is a conditional slice, not a profile likelihood: no re-optimisation of
# the other parameters. That is what makes it cheap enough to run by default.
# It therefore cannot name a ridge, only detect that a parameter is on one --
# `fit$identifiability` (component 03) names the direction.
.ctReportProfile <- function(fit, npoints = 9L, maxpars = 50L, spread = 2) {
  lpg <- .ctReportLpg(fit)
  if (is.null(lpg)) stop("This fit does not carry a likelihood that can be re-evaluated.",
    call. = FALSE)
  centre <- .ctFitRawEstimate(fit)
  if (is.null(centre) || !length(centre)) {
    stop("This fit has no point estimate to profile around.", call. = FALSE)
  }
  centre <- as.numeric(centre)
  npar <- length(centre)

  ses <- .ctReportRawSE(fit)
  fit_se <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == npar) return(x)
    rep(stats::median(x), npar)
  }
  marginal <- fit_se(ses$marginal)
  conditional <- fit_se(ses$conditional)
  # Without a Hessian there is no conditional width; fall back to the marginal
  # one and say so, rather than silently reporting a ratio of 1.
  haveboth <- !is.null(marginal) && !is.null(conditional)
  scale <- if (!is.null(conditional)) conditional else marginal
  if (is.null(scale)) scale <- rep(0.1, npar)

  nms <- try(.ctFitRawParNames(fit), silent = TRUE)
  if (inherits(nms, "try-error") || length(nms) != npar) {
    nms <- paste0("raw[", seq_len(npar), "]")
  }
  # Unnamed: a named length-1 vector recycled into data.frame() becomes row
  # names, which warns and then discards them.
  nms <- unname(as.character(nms))

  # Which parameters, when there are too many to walk them all. The flattest
  # coordinates first: those are the ones a slice is most likely to say
  # something about.
  ord <- order(scale, decreasing = TRUE)
  truncated <- npar > maxpars
  which <- sort(ord[seq_len(min(maxpars, npar))])

  npoints <- max(3L, as.integer(npoints))
  if (npoints %% 2L == 0L) npoints <- npoints + 1L   # keep the centre on the grid
  k <- seq(-spread, spread, length.out = npoints)
  expected <- spread^2 / 2
  base <- as.numeric(lpg(centre))[1L]

  grid <- vector("list", length(which))
  tab <- data.frame(parameter = nms[which], estimate = centre[which],
    se_report = if (is.null(marginal)) NA_real_ else marginal[which],
    se_slice = scale[which],
    ridge = if (haveboth) marginal[which] / conditional[which] else NA_real_,
    drop_lo = NA_real_, drop_hi = NA_real_, peak = NA_real_,
    flag = "", stringsAsFactors = FALSE)

  for (i in seq_along(which)) {
    j <- which[i]
    ll <- vapply(k, function(kk) {
      if (kk == 0) return(base)
      p <- centre
      p[j] <- p[j] + kk * scale[j]
      v <- suppressWarnings(as.numeric(lpg(p))[1L])
      if (!is.finite(v)) NA_real_ else v
    }, numeric(1))
    # -1e100 is both backends' "invalid point" sentinel, not a likelihood.
    ll[is.finite(ll) & ll < -1e99] <- NA_real_
    grid[[i]] <- data.frame(parameter = nms[j], k = k,
      value = centre[j] + k * scale[j], ll = ll, stringsAsFactors = FALSE)
    tab$drop_lo[i] <- base - ll[1L]
    tab$drop_hi[i] <- base - ll[npoints]
    if (any(is.finite(ll))) tab$peak[i] <- k[which.max(ifelse(is.finite(ll), ll, -Inf))]

    flags <- character(0)
    if (is.finite(tab$peak[i]) && abs(tab$peak[i]) > 0.5) flags <- c(flags, "off-peak")
    drops <- c(tab$drop_lo[i], tab$drop_hi[i])
    if (any(is.finite(drops))) {
      worst <- suppressWarnings(max(abs(drops[is.finite(drops)] - expected)))
      if (worst > expected / 2) flags <- c(flags, "non-quadratic")
    }
    dif <- abs(tab$drop_lo[i] - tab$drop_hi[i])
    if (is.finite(dif) && dif > expected / 2) flags <- c(flags, "asymmetric")
    if (is.finite(tab$ridge[i]) && tab$ridge[i] > 3) flags <- c(flags, "ridged")
    tab$flag[i] <- paste(flags, collapse = ",")
  }

  out <- list(table = tab, grid = do.call(rbind, grid), base = base,
    spread = spread, expected = expected, npar = npar, truncated = truncated,
    haveboth = haveboth)
  class(out) <- "ctReportProfile"
  out
}

# The slices, with the quadratic each parameter's standard error implies drawn
# over them. Reading the two curves against each other is the whole diagnostic,
# so they go on the same panel.
.ctReportProfilePlot <- function(x) {
  d <- x$grid
  d$reference <- x$base - (d$k^2) / 2
  d$Curve <- "likelihood"
  ref <- d
  ref$ll <- ref$reference
  ref$Curve <- "quadratic from SE"
  d <- rbind(d, ref)
  flag <- x$table$flag[match(d$parameter, x$table$parameter)]
  d$parameter <- factor(ifelse(nzchar(flag), paste0(d$parameter, " (", flag, ")"),
    d$parameter))
  k <- ll <- Curve <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = k, y = ll, colour = Curve, linetype = Curve)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("likelihood" = "black",
      "quadratic from SE" = "red")) +
    ggplot2::labs(x = "Conditional standard errors from the estimate",
      y = "Log probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 7))
}


# --- identification and convergence ------------------------------------------

# Everything the fit already knows about whether the optimizer arrived and
# whether the data pin the parameters, in one place. Backend-specific because
# the two backends record different things -- the julia path computes an
# eigen decomposition of the Hessian at fit time and stores it, the stan path
# stores the Hessian itself.
.ctReportIdentificationLines <- function(fit, summ) {
  L <- character(0)
  add <- function(...) L <<- c(L, paste0(...))

  add("Convergence and identification")
  add("")
  add("Backend: ", if (inherits(fit, "ctJuliaFit")) "julia" else "stan")
  if (!is.null(summ$loglik)) add("Log likelihood: ", .ctReportNum(summ$loglik))
  if (!is.null(summ$npars)) add("Free parameters: ", summ$npars)
  if (!is.null(summ$aic)) add("AIC: ", .ctReportNum(summ$aic))
  add("")

  if (inherits(fit, "ctJuliaFit")) {
    e <- fit$estimate
    add("Optimizer")
    add("  converged: ", isTRUE(e$converged))
    if (!is.null(e$iterations)) add("  iterations: ", e$iterations)
    if (!is.null(e$gradient_norm)) {
      add("  gradient norm: ", .ctReportNum(e$gradient_norm),
        if (!is.null(e$gradient_tolerance)) paste0("  (tolerance ",
          .ctReportNum(e$gradient_tolerance), ")") else "")
    }
    if (isTRUE(e$stalled)) add("  stalled: TRUE  -- the line search stopped finding a step")
    if (isTRUE(e$saturated)) {
      add("  transform flat at the estimate: ",
        paste(e$saturated_parameters, collapse = ", "))
      # The two outcomes that share a zero gradient. TRUE is a failed fit;
      # FALSE is a maximum with those coordinates unidentified.
      add("  optimizer overstepped into it: ", isTRUE(e$overshot),
        if (isTRUE(is.finite(e$overshoot_gain)))
          paste0("  (best pullback gains ", .ctReportNum(e$overshoot_gain), ")")
        else "")
    }
    if (!is.null(e$linesearch)) add("  line search: ", e$linesearch)
    add("")
    idf <- fit$identifiability
    if (!is.null(idf)) {
      add("Information matrix")
      if (!is.null(idf$condition)) add("  condition number: ", .ctReportNum(idf$condition))
      add("  weak directions: ", if (is.null(idf$nweak)) "unknown" else idf$nweak)
      if (!is.null(idf$directions) && length(idf$directions)) {
        add("")
        for (di in seq_along(idf$directions)) {
          dd <- idf$directions[[di]]
          add("  direction ", di, "  eigenvalue ", .ctReportNum(dd$eigenvalue))
          if (!is.null(dd$parameters)) {
            ord <- order(abs(dd$loadings), decreasing = TRUE)
            ord <- ord[seq_len(min(6L, length(ord)))]
            add("    ", paste0(dd$parameters[ord], " ",
              .ctReportNum(dd$loadings[ord]), collapse = "   "))
          }
        }
      }
      add("")
    }
    ivc <- fit$uncertainty$intervalcheck
    if (!is.null(ivc) && nrow(ivc$table)) {
      add("Reported interval against the curvature at the estimate")
      add("  ratio > ", ivc$threshold, ": ", ivc$nflagged, " parameter(s)")
      worst <- utils::head(ivc$table, 5)
      for (i in seq_len(nrow(worst))) {
        add("  ", worst$param[i], "  se ", .ctReportNum(worst$se[i]),
          "  curvature ", .ctReportNum(worst$curvature_se[i]),
          "  ratio ", .ctReportNum(worst$ratio[i]))
      }
      add("")
    }
    if (!is.null(fit$collapsedScales) && nrow(fit$collapsedScales)) {
      add("Population SDs at the transform floor (no individual variation left):")
      L <<- c(L, .ctReportCapture(print(fit$collapsedScales, row.names = FALSE)))
      add("")
    }
  } else {
    of <- fit$stanfit$optimfit
    add("Optimizer")
    if (!is.null(of$iter)) add("  iterations: ", of$iter)
    if (!is.null(of$ginfn)) add("  gradient (max abs): ", .ctReportNum(of$ginfn))
    if (!is.null(of$terminate)) {
      # A list of criterion/value on this path, a bare string on others.
      trm <- unlist(of$terminate)
      lab <- names(trm)
      add("  termination: ", paste0(if (is.null(lab)) "" else paste0(lab, "="),
        as.character(trm), collapse = ", "))
    }
    add("")
    # The stan path stores the Hessian but not an eigen summary of it, so
    # compute the same two numbers .ctBackendIdentifiability() reports for
    # julia rather than leaving the row blank.
    h <- fit$stanfit$uncertainty$hessian
    if (!is.null(h) && is.matrix(h) && nrow(h) == ncol(h)) {
      info <- -(h + t(h)) / 2
      ev <- suppressWarnings(try(eigen(info, symmetric = TRUE, only.values = TRUE)$values,
        silent = TRUE))
      if (!inherits(ev, "try-error")) {
        pos <- ev[ev > 0]
        add("Information matrix")
        if (length(pos)) add("  condition number: ", .ctReportNum(max(pos) / min(pos)))
        add("  non-positive eigenvalues: ", sum(ev <= 0), " of ", length(ev))
        add("")
      }
    }
  }

  if (!is.null(summ$residCovStdConditioning)) {
    add("Residual covariance: ", summ$residCovStdConditioning)
    add("")
  }

  add("Read: a condition number above about 1e6, a non-positive eigenvalue, or")
  add("any weak direction means the data do not pin every parameter. Estimates")
  add("can then be moved along that direction without worsening the fit, so")
  add("their intervals describe a shape the data did not choose. The ratio")
  add("above says the same thing per parameter: how much wider the reported")
  add("interval is than that parameter's own curvature supports. Healthy fits")
  add("sit near 1; a ratio of 100 means the width is entanglement with the")
  add("other parameters, and it will not repeat between runs.")
  add("04-loglik-profile says which individual parameters are involved.")
  L
}

.ctReportNum <- function(x) {
  if (is.null(x) || !length(x)) return("NA")
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return("NA")
  format(x, digits = 6, trim = TRUE)
}


# --- the index ----------------------------------------------------------------

# One line saying what the file is, one saying how to read it. Nothing longer:
# a user who wants the mechanism has the help page, and a paragraph here would
# be skipped whole.
.ctReportInterpretation <- function() list(
  summary = c("Parameter estimates with uncertainty intervals.",
    "An interval excluding zero is the usual claim of an effect; check 03 first."),
  parmatrices = c("System matrices at each requested quantile.",
    "DRIFT diagonals negative means stable; off-diagonals are cross-effects per unit time."),
  identification = c("Whether the optimizer arrived, and whether the data pin the parameters.",
    "Read this first. A weak direction or a huge condition number invalidates the intervals in 01."),
  profile = c("Log-probability slice through each parameter against the quadratic its curvature implies.",
    "Off-peak: not converged. Ridged: pinned alone, free jointly. Non-quadratic or asymmetric: the interval is the wrong shape."),
  discretepars = c("Regression coefficients between latents as a function of time interval.",
    "This is what the continuous-time DRIFT means for a given gap between measurements."),
  network = c("The temporal and contemporaneous networks at one time interval.",
    "An arrow is a directed effect over that interval; an undirected edge is shared innovation."),
  predictions = c("Observed data with model expectations, for the most-observed subjects.",
    "Expectations conditional on covariates only show what the model predicts blind; conditional on all data, what it can reconstruct."),
  residuals = c("Standardised prior residuals, and their autocorrelation.",
    "Autocorrelation outside the interval means structure the model has not taken up."),
  covcheck = c("Lagged covariances, observed against model-generated.",
    "Empirical outside the generated interval is misfit at that lag; a whole row or column off is a variable the model gets wrong."),
  postpred = c("Distribution of generated data against observed.",
    "Observed away from the generated mass means the model cannot produce data like yours."),
  tipredeffects = c("Predicted parameter values across the range of each covariate.",
    "A band that excludes a flat line is a covariate effect on that parameter."),
  crossval = c("K-fold out-of-sample log likelihood.",
    "Compare across models on the same data and folds; the in/out gap is overfitting."),
  equations = c("The fitted model written out as equations.",
    "For a methods section, and for checking the model is the one you meant.")
)

.ctReportIndex <- function(fit, summ, record, folder) {
  interp <- .ctReportInterpretation()
  L <- character(0)
  add <- function(...) L <<- c(L, paste0(...))
  # A component with no interpretation entry would otherwise write its heading
  # and then nothing: paste0(NULL) is character(0), which appends no line and
  # raises nothing. A visible placeholder instead, so the gap is in the file
  # rather than in the reader's understanding.
  say <- function(nm, i) {
    x <- interp[[nm]][i]
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      return(paste0("(ctReport has no description for component '", nm,
        "' -- please report this.)"))
    }
    x
  }

  m <- .ctFitModelObject(fit)
  add("# ctsem report")
  add("")
  add("Backend: ", if (inherits(fit, "ctJuliaFit")) "julia" else "stan",
    ". Latents: ", m$n.latent, ". Manifests: ", m$n.manifest,
    ". Covariates: ", m$n.TIpred, ".")
  add("Subjects: ", .ctFitNsubjects(fit), ". Rows: ", length(.ctFitRowTime(fit)),
    ". Uncertainty draws: ", .ctFitNsamples(fit), ".")
  if (!is.null(summ$loglik)) {
    add("Log likelihood: ", .ctReportNum(summ$loglik),
      ". Free parameters: ", if (is.null(summ$npars)) "NA" else summ$npars,
      ". AIC: ", .ctReportNum(summ$aic), ".")
  }
  add("Written ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by ctReport().")
  add("")

  done <- record[vapply(record, function(x) x$status == "ok", logical(1))]
  if (length(done)) {
    add("## Files")
    add("")
    for (nm in names(done)) {
      r <- done[[nm]]
      add("### ", paste(r$files, collapse = ", "))
      add("")
      add(say(nm, 1))
      add("")
      add("*", say(nm, 2), "*")
      add("")
      if (!is.null(r$note)) {
        add(r$note)
        add("")
      }
    }
  }

  other <- record[vapply(record, function(x) x$status != "ok", logical(1))]
  if (length(other)) {
    add("## Not produced")
    add("")
    for (nm in names(other)) {
      r <- other[[nm]]
      add("- **", nm, "** -- ", if (r$status == "failed") "failed: " else "", r$note)
    }
    add("")
  }

  add("---")
  add("")
  add("Each file comes from one ctsem function; the help page for that function")
  add("has the arguments and the detail. Regenerate with `ctReport(fit)`.")
  .ctReportWrite(file.path(folder, "00-index.md"), L)
  invisible(L)
}


# --- the components -----------------------------------------------------------
#
# Each takes (fit, ctx) and returns a character vector of file names, or calls
# ctx$skip(reason) and returns NULL. ctx carries the folder, the shared summary
# object, and the run's arguments.

.ctReportDo <- list(

  summary = function(fit, ctx) {
    f <- "01-summary.txt"
    L <- .ctReportCapture(print(ctx$summ))
    s <- ctx$summ
    extra <- character(0)
    # The rows a reader would otherwise have to find by eye. Not a new
    # statistic -- a filter on the summary that is already printed above.
    if (!is.null(s$popmeans)) {
      pm <- s$popmeans
      qcols <- intersect(c("2.5%", "97.5%"), colnames(pm))
      if (length(qcols) == 2L) {
        keep <- sign(pm[, qcols[1]]) == sign(pm[, qcols[2]])
        keep[is.na(keep)] <- FALSE
        extra <- c(extra, "", "Population means whose interval excludes zero:",
          if (any(keep)) .ctReportCapture(print(pm[keep, , drop = FALSE])) else "  none")
      }
    }
    for (blk in c("tipreds", "rawpopcorr")) {
      b <- s[[blk]]
      if (is.null(b) || !nrow(b)) next
      qcols <- intersect(c("2.5%", "97.5%"), colnames(b))
      if (length(qcols) != 2L) next
      keep <- sign(b[, qcols[1]]) == sign(b[, qcols[2]])
      keep[is.na(keep)] <- FALSE
      lab <- if (blk == "tipreds") "Covariate effects" else "Random effect correlations"
      extra <- c(extra, "", paste0(lab, " whose interval excludes zero:"),
        if (any(keep)) .ctReportCapture(print(b[keep, , drop = FALSE])) else "  none")
    }
    .ctReportWrite(file.path(ctx$folder, f), c(L, extra))
    f
  },

  parmatrices = function(fit, ctx) {
    f <- "02-parameter-matrices.txt"
    probs <- sort(unique(ctx$quantiles))
    L <- character(0)
    for (p in probs) {
      mats <- suppressMessages(ctSummaryMatrices(fit, calcfuncargs = list(probs = p),
        timeinterval = ctx$timeinterval))
      L <- c(L, paste0("=== quantile ", p, " ==="), "",
        .ctReportCapture(print(mats)), "")
    }
    .ctReportWrite(file.path(ctx$folder, f), L)
    f
  },

  identification = function(fit, ctx) {
    f <- "03-identification.txt"
    .ctReportWrite(file.path(ctx$folder, f),
      .ctReportIdentificationLines(fit, ctx$summ))
    f
  },

  profile = function(fit, ctx) {
    if (is.null(.ctReportLpg(fit))) {
      ctx$skip("this fit does not carry a re-evaluable likelihood")
      return(NULL)
    }
    pr <- .ctReportProfile(fit, npoints = ctx$profilepoints)
    tab <- pr$table
    for (cc in c("estimate", "se_report", "se_slice")) tab[[cc]] <- round(tab[[cc]], 4)
    for (cc in c("ridge", "drop_lo", "drop_hi")) tab[[cc]] <- round(tab[[cc]], 2)
    L <- c("Log-probability slices through each raw parameter, one at a time,",
      "with the others held at the estimate.",
      "",
      "se_report  the standard error the summary reports (all else unknown).",
      "se_slice   the standard error this parameter's own curvature implies",
      "           (all else fixed). The slice grid is in these units.",
      "ridge      se_report / se_slice. Above 3, most of the reported",
      "           uncertainty comes from correlation with other parameters.",
      paste0("drop_lo    fall in log probability ", pr$spread,
        " se_slice below the estimate; a"),
      paste0("drop_hi    quadratic surface gives ", pr$expected, " on both sides."),
      "peak       where the slice actually peaks, in se_slice.",
      "",
      "Flags: off-peak (the optimizer is not at the maximum), non-quadratic",
      "(the surface is not the shape the standard error assumes), asymmetric",
      "(so a symmetric interval is the wrong shape), ridged (see above).",
      "")
    if (!pr$haveboth) {
      L <- c(L, "No Hessian on this fit, so se_slice is the reported SE and ridge",
        "cannot be computed. The drops then exceed 2 wherever a parameter is",
        "correlated with others, and only off-peak still reads cleanly. Refit",
        "with uncertainty for the rest.", "")
    }
    if (pr$truncated) {
      L <- c(L, paste0("Note: ", nrow(tab), " of ", pr$npar,
        " parameters shown, those with the flattest slices."), "")
    }
    L <- c(L, .ctReportCapture(print(tab, row.names = FALSE)))
    flagged <- tab[nzchar(tab$flag), , drop = FALSE]
    L <- c(L, "", if (nrow(flagged)) {
      c("Flagged:", .ctReportCapture(print(flagged[, c("parameter", "flag")],
        row.names = FALSE)))
    } else "Nothing flagged.")
    .ctReportWrite(file.path(ctx$folder, "04-loglik-profile.txt"), L)
    .ctReportPdf(file.path(ctx$folder, "04-loglik-profile.pdf"),
      print(.ctReportProfilePlot(pr)))
    if (nrow(flagged)) {
      ctx$note(paste0(nrow(flagged), " of ", nrow(tab), " parameters flagged: ",
        paste(unique(unlist(strsplit(flagged$flag, ","))), collapse = ", "), "."))
    }
    c("04-loglik-profile.txt", "04-loglik-profile.pdf")
  },

  discretepars = function(fit, ctx) {
    f <- "05-discrete-parameters.pdf"
    times <- if (identical(ctx$times, "auto")) .ctReportTimes(fit) else ctx$times
    dp <- suppressMessages(ctDiscretePars(fit, times = times, plot = FALSE,
      nsamples = ctx$nsamples, cores = ctx$cores))
    dpo <- suppressMessages(try(ctDiscretePars(fit, times = times, plot = FALSE,
      nsamples = ctx$nsamples, cores = ctx$cores, observational = TRUE), silent = TRUE))
    .ctReportPdf(file.path(ctx$folder, f), {
      print(ctDiscreteParsPlot(dp, quantiles = ctx$quantiles,
        title = "Effect of an intervention on a latent"))
      if (!inherits(dpo, "try-error")) {
        print(ctDiscreteParsPlot(dpo, quantiles = ctx$quantiles,
          title = "Effect of an observation of a latent"))
      }
    })
    f
  },

  network = function(fit, ctx) {
    f <- "06-network.pdf"
    # engine = 'ggplot' rather than qgraph, which is a Suggests.
    net <- suppressMessages(ctNetwork(fit, dt = ctx$timeinterval, plot = FALSE))
    .ctReportPdf(file.path(ctx$folder, f),
      print(suppressMessages(ctNetworkPlot(net, dt = ctx$timeinterval,
        engine = "ggplot"))))
    utils::write.csv(net$edges, file.path(ctx$folder, "06-network-edges.csv"),
      row.names = FALSE)
    c(f, "06-network-edges.csv")
  },

  predictions = function(fit, ctx) {
    f <- "07-predicted-trajectories.pdf"
    subs <- .ctReportSubjects(fit, ctx$nsubjects)
    if (!length(subs)) {
      ctx$skip("no subjects could be identified in the fit")
      return(NULL)
    }
    k <- suppressWarnings(suppressMessages(ctPredict(fit, subjects = subs)))
    krem <- suppressWarnings(suppressMessages(ctPredict(fit, subjects = subs,
      removeObs = TRUE)))
    .ctReportPdf(file.path(ctx$folder, f), {
      print(plot(krem, kalmanvec = "yprior", polygonsteps = FALSE, plot = FALSE) +
          ggplot2::ggtitle("Expectations given covariates only"))
      print(plot(k, kalmanvec = "ysmooth", polygonsteps = FALSE, plot = FALSE) +
          ggplot2::ggtitle("Expectations given all data"))
      print(plot(krem, kalmanvec = "etaprior", polygonsteps = FALSE, plot = FALSE) +
          ggplot2::ggtitle("Latent expectations given covariates only"))
      print(plot(k, kalmanvec = "etasmooth", polygonsteps = FALSE, plot = FALSE) +
          ggplot2::ggtitle("Latent expectations given all data"))
    })
    ctx$note(paste0("Subjects shown: ", paste(subs, collapse = ", "), "."))
    f
  },

  residuals = function(fit, ctx) {
    res <- suppressWarnings(suppressMessages(ctResiduals(fit)))
    # plotctACF() separately rather than through ctACF(plot = TRUE), because
    # its smoothing spline needs more distinct time intervals than a short
    # panel has, and ctACF() gives no way to reach the argument that turns it
    # off. plotctACF() now falls back to the raw samples itself and says so in a
    # message, so catch that message to write the note; the try() stays as
    # defence against a failure it does not handle.
    ac <- suppressWarnings(suppressMessages(
      ctACFresiduals(fit, plot = FALSE, nboot = 20)))
    .ctReportPdf(file.path(ctx$folder, "08-residuals.pdf"), {
      acfmsg <- character(0)
      ok <- try(withCallingHandlers(suppressWarnings(print(plotctACF(ac))),
        message = function(m) {
          acfmsg <<- c(acfmsg, conditionMessage(m))
          invokeRestart("muffleMessage")
        }), silent = TRUE)
      if (inherits(ok, "try-error")) {
        suppressWarnings(suppressMessages(print(plotctACF(ac, estimateSpline = FALSE))))
        ctx$note("Too few distinct time intervals for the smoothing spline; samples only.")
      } else if (any(grepl("plotting ACF samples", acfmsg, fixed = TRUE))) {
        ctx$note("Too few distinct time intervals for the smoothing spline; samples only.")
      }
    })
    vars <- setdiff(names(res), c("Subject", "Time"))
    L <- c("Standardised prior residuals: distribution per variable.",
      "Mean near 0 and SD near 1 is what the model claims.", "",
      .ctReportCapture(print(round(data.frame(
        mean = vapply(vars, function(v) mean(res[[v]], na.rm = TRUE), numeric(1)),
        sd = vapply(vars, function(v) stats::sd(res[[v]], na.rm = TRUE), numeric(1)),
        min = vapply(vars, function(v) min(res[[v]], na.rm = TRUE), numeric(1)),
        max = vapply(vars, function(v) max(res[[v]], na.rm = TRUE), numeric(1)),
        n = vapply(vars, function(v) sum(!is.na(res[[v]])), numeric(1))
      ), 3))))
    .ctReportWrite(file.path(ctx$folder, "08-residuals.txt"), L)
    c("08-residuals.pdf", "08-residuals.txt")
  },

  covcheck = function(fit, ctx) {
    cc <- suppressWarnings(suppressMessages(ctFitCheckCov(ctx$generated(),
      plot = FALSE, nsamples = ctx$nsamples, cores = ctx$cores, lags = ctx$lags)))
    utils::write.csv(cc, file.path(ctx$folder, "09-covariance-check.csv"),
      row.names = FALSE)
    .ctReportPdf(file.path(ctx$folder, "09-covariance-check.pdf"), {
      p <- ctFitCovCheckPlot(cc, cor = isTRUE(attr(cc, "ctFitCovCheck_cor")))
      if (inherits(p, "list")) lapply(p, print) else print(p)
    })
    if ("Sig" %in% names(cc)) {
      nsig <- sum(cc$Sig %in% TRUE)
      ctx$note(paste0(nsig, " of ", nrow(cc),
        " observed covariances fall outside the generated interval."))
    }
    c("09-covariance-check.csv", "09-covariance-check.pdf")
  },

  postpred = function(fit, ctx) {
    f <- "10-posterior-predictive.pdf"
    pp <- suppressWarnings(suppressMessages(ctPostPredPlots(ctx$generated())))
    .ctReportPdf(file.path(ctx$folder, f), lapply(pp, function(p) try(print(p), silent = TRUE)))
    f
  },

  tipredeffects = function(fit, ctx) {
    m <- .ctFitModelObject(fit)
    if (is.null(m$n.TIpred) || m$n.TIpred < 1) {
      ctx$skip("the model has no covariates (time independent predictors)")
      return(NULL)
    }
    f <- "11-covariate-effects.pdf"
    .ctReportPdf(file.path(ctx$folder, f), {
      for (tip in seq_len(m$n.TIpred)) {
        te <- suppressWarnings(suppressMessages(ctTIpredEffects(fit,
          whichTIpreds = tip, nsamples = ctx$nsamples, plot = FALSE,
          timeinterval = ctx$timeinterval)))
        # Parameters the covariate does not move carry no information and
        # crowd the panels out.
        keep <- apply(te$y[, , 2, drop = FALSE], 2,
          function(x) length(unique(round(x, 8))) > 1)
        if (!any(keep)) next
        te$y <- te$y[, keep, , drop = FALSE]
        mats <- unique(gsub("\\[.*", "", colnames(te$y)))
        for (mm in mats) {
          sub <- te
          sub$y <- te$y[, grep(paste0("^", mm, "\\["), colnames(te$y)), , drop = FALSE]
          if (!dim(sub$y)[2]) next
          print(ctPlotArrayGG(sub) +
              ggplot2::ggtitle(paste0(m$TIpredNames[tip], " on ", mm)))
        }
      }
    })
    f
  },

  crossval = function(fit, ctx) {
    f <- "12-cross-validation.txt"
    lo <- suppressWarnings(suppressMessages(ctLOO(fit, folds = ctx$folds,
      cores = ctx$cores)))
    keep <- c("insampleLogLik", "outsampleLogLik",
      "insampleRowwiseEntropy", "outsampleRowwiseEntropy",
      "insampleSubjectwiseEntropy", "outsampleSubjectwiseEntropy")
    keep <- keep[keep %in% names(lo)]
    L <- c(paste0(ctx$folds, "-fold cross validation."), "",
      .ctReportCapture(print(data.frame(quantity = keep,
        value = vapply(keep, function(x) as.numeric(lo[[x]])[1], numeric(1))),
        row.names = FALSE)),
      "", "Out-of-sample entropy is the comparable number across models fitted",
      "to the same data with the same folds. Lower is better.")
    .ctReportWrite(file.path(ctx$folder, f), L)
    f
  },

  equations = function(fit, ctx) {
    stem <- "13-model-equations"
    write <- function(x, compile) suppressWarnings(suppressMessages(
      ctModelLatex(x, folder = ctx$folder, filename = stem, open = FALSE,
        compile = compile)))
    notes <- character(0)
    # One reason this can fall short of "the fitted model, as a compiled pdf":
    # compiling needs a LaTeX installation. Writing the .tex does not, so a
    # machine without one still gets the equations with the estimates in them.
    out <- try(write(fit, TRUE), silent = TRUE)
    if (inherits(out, "try-error")) write(fit, FALSE)
    files <- paste0(stem, c(".tex", ".pdf"))
    files <- files[file.exists(file.path(ctx$folder, files))]
    if (!any(grepl("[.]pdf$", files))) {
      notes <- c(notes, "LaTeX was not available to compile it; the .tex source is here.")
    }
    if (length(notes)) ctx$note(paste(notes, collapse = " "))
    files
  }
)


#' Write a folder of summary, diagnostic and interpretation output for a fit
#'
#' Runs the package's reporting and fit-checking functions over a fitted model
#' and writes the results into one folder, with a markdown index that says what
#' each file is and how to read it. Works for either backend; components a
#' backend or a model does not support are recorded in the index as not
#' produced, rather than stopping the run.
#'
#' @param fit A fit from \code{\link{ctFit}} (class \code{ctStanFit} or
#'   \code{ctJuliaFit}).
#' @param folder Directory to write into; created if it does not exist. Files
#'   from a previous \code{ctReport} run in the same folder are removed first,
#'   so the folder always matches its index.
#' @param components Either one of \code{'quick'} (nothing that needs generated
#'   data), \code{'default'}, \code{'all'} (adds cross validation and LaTeX), or
#'   a character vector of component names: \code{summary}, \code{parmatrices},
#'   \code{identification}, \code{profile}, \code{discretepars},
#'   \code{network}, \code{predictions}, \code{residuals}, \code{covcheck},
#'   \code{postpred}, \code{tipredeffects}, \code{crossval},
#'   \code{equations}.
#' @param nsamples Draws used by the components that sample from the fit's
#'   uncertainty, and by data generation for \code{covcheck} and
#'   \code{postpred}.
#' @param cores Cores for the components that can use more than one
#'   (\code{discretepars}, \code{covcheck}, \code{postpred}, \code{crossval}).
#'   Defaults to 1: nothing parallelises unless asked.
#' @param quantiles Quantiles for the parameter matrices and the discrete
#'   parameter plots.
#' @param times Interval grid for \code{discretepars}. \code{'auto'} uses zero
#'   to the 90th percentile of the observed per-subject time span.
#' @param timeinterval Time interval at which the networks, the discrete-time
#'   parameter matrices and the covariate effects are reported.
#' @param nsubjects Number of subjects plotted by \code{predictions}, chosen as
#'   those with the most observations.
#' @param lags Lags for the lagged covariance check.
#' @param folds Folds for \code{crossval}.
#' @param profilepoints Grid points per parameter in the likelihood profile.
#' @param verbose Report each component as it runs.
#'
#' @return Invisibly, a list with the folder and, per component, whether it was
#'   written, skipped or failed. Called for the files it writes.
#'
#' @details
#' The expensive parts are \code{covcheck} and \code{postpred}, which generate
#' data from the fit, and \code{crossval}, which refits. \code{components =
#' 'quick'} excludes all three.
#'
#' \code{profile} walks the log probability along each raw parameter, holding
#' the others at the estimate, over a grid measured in the standard error that
#' parameter's own curvature implies. It flags parameters whose slice peaks away
#' from the estimate (the optimizer did not arrive), whose surface is not the
#' shape a standard error assumes, and whose reported interval is much wider
#' than their own curvature (they sit on a ridge with other parameters). It is a
#' conditional slice, not a profile likelihood, which is why it is cheap enough
#' to run by default; \code{identification} names the ridge directions.
#'
#' @examples
#' \donttest{
#' f <- ctReport(ctstantestfit, folder = file.path(tempdir(), 'ctReport'),
#'   components = c('summary', 'identification', 'profile'), profilepoints = 5)
#' list.files(f$folder)
#' }
#' @export
ctReport <- function(fit, folder = "ctReport", components = "default",
  nsamples = 200, cores = 1, quantiles = c(.025, .5, .975), times = "auto",
  timeinterval = 1, nsubjects = 6, lags = 0:3, folds = 10, profilepoints = 9,
  verbose = TRUE) {

  if (!inherits(fit, c("ctStanFit", "ctJuliaFit"))) {
    stop("fit must be a ctStanFit or ctJuliaFit from ctFit().", call. = FALSE)
  }
  # Default 1: a report is run once, and nothing here should take cores the
  # user did not offer. An explicit request is passed through -- silently
  # capping it would read as a bug on a machine that has the cores.
  cores <- suppressWarnings(as.integer(cores[1]))
  if (!isTRUE(is.finite(cores)) || cores < 1L) cores <- 1L
  wanted <- .ctReportResolve(components)

  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  } else {
    # Only this function's own numbered files, so pointing ctReport at a
    # directory that has other things in it does not delete them.
    old <- list.files(folder, pattern = "^[0-9][0-9]-", full.names = TRUE)
    if (length(old)) unlink(old)
  }
  folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)

  summ <- suppressWarnings(suppressMessages(summary(fit, timeinterval = timeinterval)))

  # Generated data, made once and shared: ctFitCheckCov and ctPostPredPlots
  # would each otherwise generate their own with their own defaults.
  genfit <- NULL
  generated <- function() {
    if (is.null(genfit)) {
      if (is.null(fit$generated)) {
        genfit <<- suppressWarnings(suppressMessages(
          ctGenerateFromFit(fit, nsamples = nsamples, cores = cores)))
      } else genfit <<- fit
    }
    genfit
  }

  record <- list()
  for (nm in wanted) {
    state <- new.env(parent = emptyenv())
    state$skipped <- NULL
    state$note <- NULL
    ctx <- list(folder = folder, summ = summ, nsamples = nsamples, cores = cores,
      quantiles = quantiles, times = times, timeinterval = timeinterval,
      nsubjects = nsubjects, lags = lags, folds = folds,
      profilepoints = profilepoints, generated = generated,
      skip = function(reason) state$skipped <- reason,
      note = function(text) state$note <- text)
    t0 <- Sys.time()
    # Messages, not warnings: the plotting stack is chatty (loess formulas,
    # unknown aesthetics) and a report should print one line per component.
    # A warning still reaches the user, and a failure is recorded either way.
    files <- try(suppressMessages(.ctReportDo[[nm]](fit, ctx)), silent = TRUE)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(files, "try-error")) {
      record[[nm]] <- list(status = "failed",
        note = trimws(conditionMessage(attr(files, "condition"))))
      if (verbose) message(sprintf("%-15s failed: %s", nm, record[[nm]]$note))
    } else if (!is.null(state$skipped)) {
      record[[nm]] <- list(status = "skipped", note = state$skipped)
      if (verbose) message(sprintf("%-15s skipped: %s", nm, state$skipped))
    } else {
      record[[nm]] <- list(status = "ok", files = files, note = state$note)
      if (verbose) message(sprintf("%-15s %.0fs", nm, el))
    }
  }

  .ctReportIndex(fit, summ, record, folder)
  if (verbose) message("Written to ", folder, ". Start at 00-index.md.")
  invisible(list(folder = folder, components = record))
}
