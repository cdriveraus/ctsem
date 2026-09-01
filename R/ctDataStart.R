# Starting values read off the data ------------------------------------------
#
# ctsem optimises on a raw, unconstrained scale, and every fit has started from
# the same fixed point on it -- `rnorm(npar, 0, 0.01)`, i.e. essentially raw
# zero. What raw zero *means* is decided by each parameter's transform, and for
# the variance family it means a lot: `10 * log1p_exp(2 * param)` puts DIFFUSION
# at 6.93 and `5 * log1p_exp(2 * param)` puts T0VAR and MANIFESTVAR at 3.47,
# whatever the data are measured in.
#
# That is a guess about the scale of the user's data, and it fails in one
# direction. Measured on a three-indicator factor model, 50 subjects x 10
# occasions, with the data multiplied by a constant:
#
#   scale x0.01   true DIFFUSION 0.015   fails: loadings -844 and -1242, |raw| 248
#   scale x1      true DIFFUSION 1.5     converges, loadings 0.808 / 1.150
#   scale x100    true DIFFUSION 150     converges, same loadings, |raw| 16
#
# Starting far above the data's scale lets the optimiser collapse the latent to
# zero variance and explain everything as measurement error, which is a real
# local optimum 216 log units below the answer. Seeding DIFFUSION alone at the
# data's scale takes that same fit from -1993.43 to -1777.27, converged, and
# exactly onto the optimum a DIFFUSION-anchored parameterisation finds.
#
# Recentring the transform instead was tried and rejected: shifting raw zero to
# 1.0 by an inner offset leaves the x0.01 case failing identically (the start is
# still 67x too large) and makes the x100 case land 186 log units worse. Any
# fixed constant only moves which scales are unlucky, because what matters is
# the ratio between the start and the data.
#
# What the manifest scale says about the latent scale depends on the link, and
# using the raw standard deviation for every type would be worse than the fixed
# default for two of them. Measured on one known latent process of sd 1.128:
#
#   indicator     sd(y)   the statistic that recovers the latent scale
#   continuous    1.268   sd(y)                        -> 1.015
#   censored      1.089   sd(y), attenuated            -> 0.871
#   count         3.919   sd(log1p(y)) / loading       -> 1.295
#   binary        0.497   nothing: the logit sets it
#   ordinal       1.192   nothing: the logit sets it
#
# A count's observations live on the rate scale while its latent lives on the
# log rate, so `sd(y)` overstates the latent scale by a factor of 6.9 here --
# in the one direction that fails. Binary and ordinal observations carry no
# scale information at all: the link fixes the latent scale, so a process
# measured only by them gets a constant instead. The ordinal figure above
# agreeing with the truth is an artefact of four categories, not signal.
#
# Deliberately modest in scope: diagonals only, of the four matrices whose
# defaults are scale guesses, by inverting each cell's own transform
# numerically so nothing here needs to know the transform's shape. Everything
# is clamped, and any cell that cannot be solved keeps the old default rather
# than getting an invented value.

# Where an unsolvable or absurd answer gets cut off. The raw bound is well
# inside the region where the standard transforms are well conditioned; the
# scale bound only excludes values that cannot be a measurement.
.ctDataStartRawBound <- 3
.ctDataStartScaleRange <- c(1e-3, 1e3)

# The share of an indicator's variance handed to the process rather than to
# measurement error. 0.8 and 0.6 because 0.8^2 + 0.6^2 = 1: a starting split of
# the observed variance, not an estimate of anything.
.ctDataStartSignal <- 0.8
.ctDataStartNoise <- 0.6

# A process measured only through a logit link has no scale in the data. The
# logistic distribution's own standard deviation is pi/sqrt(3) = 1.81; 1.5 is
# the same order and slightly conservative, which is the safe side here.
.ctDataStartLinkScale <- 1.5

#' The scale statistic for one manifest variable, by its type.
#'
#' NA means "this indicator says nothing about the latent scale", which is the
#' honest answer for binary and ordinal rather than a failure.
#' @keywords internal
.ctDataStartManifestScale <- function(y, manifesttype) {
  y <- y[is.finite(y)]
  if (length(y) < 3L) return(NA_real_)
  s <- switch(as.character(manifesttype),
    # Gaussian, and censored which is Gaussian inside its limits. Censoring
    # attenuates the spread, which errs low, which is the safe direction.
    "0" = stats::sd(y),
    "4" = stats::sd(y),
    # The latent is the log rate, so the spread has to be taken there.
    "3" = stats::sd(log1p(pmax(y, 0))),
    NA_real_)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  s
}

#' Numeric LAMBDA, with free and state-dependent cells taken as 1.
#'
#' Only used to divide an indicator's spread back to its process, so a loading
#' that is not yet a number is better treated as unity than dropped.
#' @keywords internal
.ctDataStartLambda <- function(spec, nmanifest, nlatent) {
  out <- matrix(0, nmanifest, nlatent)
  rows <- spec$parameter_table[spec$parameter_table$matrix %in% "LAMBDA", , drop = FALSE]
  for (i in seq_len(nrow(rows))) {
    r <- rows$row[i]; cc <- rows$col[i]
    if (is.na(r) || is.na(cc) || r > nmanifest || cc > nlatent) next
    v <- rows$value[i]
    out[r, cc] <- if (is.na(v)) 1 else v
  }
  out
}

#' Pooled within-subject lag-1 correlation, and the median interval.
#'
#' Successive pairs within a subject, pooled. Crude by design: it only has to
#' put the drift in the right decade, and it is clamped afterwards.
#' @keywords internal
.ctDataStartAutocor <- function(datalong, column, idname, timename) {
  d <- datalong[order(datalong[[idname]], datalong[[timename]]), , drop = FALSE]
  y <- suppressWarnings(as.numeric(d[[column]]))
  id <- d[[idname]]
  t <- suppressWarnings(as.numeric(d[[timename]]))
  same <- c(FALSE, id[-1] == id[-length(id)])
  prev <- c(NA, utils::head(y, -1))
  dt <- c(NA, diff(t))
  keep <- same & is.finite(y) & is.finite(prev) & is.finite(dt) & dt > 0
  if (sum(keep) < 5L) return(NULL)
  r <- suppressWarnings(stats::cor(y[keep], prev[keep]))
  if (!is.finite(r)) return(NULL)
  list(r = r, dt = stats::median(dt[keep]))
}

#' Solve transform(raw) = target for raw, without knowing the transform.
#'
#' The transform arrives as text using `param[k]`, so it is evaluated rather
#' than inverted symbolically -- which means a transform this function has
#' never seen still works, and one it cannot bracket returns NA rather than a
#' guess.
#' @keywords internal
.ctDataStartInvert <- function(text, target, bound = .ctDataStartRawBound) {
  if (is.na(text) || !nzchar(text) || !is.finite(target)) return(NA_real_)
  expr <- try(parse(text = gsub("param\\s*\\[\\s*\\d+\\s*\\]", "param",
    as.character(text))), silent = TRUE)
  if (inherits(expr, "try-error")) return(NA_real_)
  env <- new.env(parent = baseenv())
  env$log1p_exp <- function(x) ifelse(x > 30, x, log1p(exp(x)))
  env$inv_logit <- function(x) 1 / (1 + exp(-x))
  f <- function(p) {
    env$param <- p
    v <- try(eval(expr[[1L]], envir = env), silent = TRUE)
    if (inherits(v, "try-error") || length(v) != 1L || !is.finite(v)) return(NA_real_)
    v - target
  }
  lo <- f(-bound); hi <- f(bound)
  if (!is.finite(lo) || !is.finite(hi) || is.na(lo) || is.na(hi)) return(NA_real_)
  # No sign change means the target is outside what this transform reaches
  # within the bound. Clamping to the nearer end is right: it moves the start
  # as far towards the data as the parameterisation allows.
  if (lo * hi > 0) return(if (abs(lo) < abs(hi)) -bound else bound)
  root <- try(stats::uniroot(f, lower = -bound, upper = bound, tol = 1e-8),
    silent = TRUE)
  if (inherits(root, "try-error")) return(NA_real_)
  max(-bound, min(bound, root$root))
}

#' Starting values derived from the data, on the raw scale.
#'
#' @return A numeric vector of length npar, NA wherever nothing was derived, or
#'   NULL when the data give nothing to work with.
#' @keywords internal
.ctDataStart <- function(datalong, model, spec, npar) {
  pt <- spec$parameter_table
  if (is.null(pt) || !npar || !all(c("matrix", "row", "col", "parnumber",
    "transform") %in% names(pt))) return(NULL)

  nlatent <- model$n.latent
  mnames <- model$manifestNames
  nmanifest <- length(mnames)
  idname <- if (!is.null(model$subjectIDname)) model$subjectIDname else "id"
  timename <- if (!is.null(model$timeName)) model$timeName else "time"
  if (!all(c(idname, timename) %in% names(datalong))) return(NULL)
  present <- mnames[mnames %in% names(datalong)]
  if (!length(present)) return(NULL)

  mtype <- spec$manifesttype
  if (is.null(mtype) || length(mtype) != nmanifest) mtype <- rep(0L, nmanifest)
  continuous <- !identical(spec$continuoustime, FALSE)

  # Per indicator: its spread on the scale its latent lives on.
  mscale <- rep(NA_real_, nmanifest)
  for (j in seq_len(nmanifest)) {
    if (!mnames[j] %in% names(datalong)) next
    mscale[j] <- .ctDataStartManifestScale(
      suppressWarnings(as.numeric(datalong[[mnames[j]]])), mtype[j])
  }

  lambda <- .ctDataStartLambda(spec, nmanifest, nlatent)

  # Per process: the spread of the indicators that load on it, divided back
  # through their loadings, then split with measurement error.
  latentsd <- rep(NA_real_, nlatent)
  anchor <- rep(NA_character_, nlatent)
  for (i in seq_len(nlatent)) {
    load <- lambda[, i]
    use <- which(load != 0 & is.finite(mscale))
    if (length(use)) {
      latentsd[i] <- stats::median(mscale[use] / abs(load[use])) * .ctDataStartSignal
      # Prefer a Gaussian indicator to read the autocorrelation from; a count
      # would have to be logged first and a censored one is clipped.
      pick <- use[order(mtype[use] != 0, -abs(load[use]))][1L]
      anchor[i] <- mnames[pick]
    } else {
      # Measured only through a link that sets the scale itself.
      latentsd[i] <- .ctDataStartLinkScale
    }
    latentsd[i] <- min(max(latentsd[i], .ctDataStartScaleRange[1L]),
      .ctDataStartScaleRange[2L])
  }

  # Per process: how much of itself survives one interval, and how long that
  # interval is. Clamped well away from a random walk and from white noise.
  ar <- rep(NA_real_, nlatent); dtmed <- rep(NA_real_, nlatent)
  for (i in seq_len(nlatent)) {
    if (is.na(anchor[i])) next
    a <- .ctDataStartAutocor(datalong, anchor[i], idname, timename)
    if (is.null(a)) next
    ar[i] <- min(max(a$r, 0.05), 0.95)
    dtmed[i] <- if (is.finite(a$dt) && a$dt > 0) a$dt else 1
  }

  target <- rep(NA_real_, npar)
  set <- function(parnumber, value) {
    if (is.na(parnumber) || parnumber < 1 || parnumber > npar) return(invisible())
    if (is.finite(value)) target[parnumber] <<- value
  }

  free <- pt[!is.na(pt$parnumber) & !is.na(pt$transform), , drop = FALSE]
  # `JAx` and `Jy` repeat DRIFT and LAMBDA under another name; setting a
  # parameter twice from two views of it is the same value, but reading the
  # model matrices only keeps this to the cells actually meant.
  free <- free[free$matrix %in% c("DRIFT", "DIFFUSION", "T0VAR", "MANIFESTVAR"), ,
    drop = FALSE]
  free <- free[!is.na(free$row) & !is.na(free$col) & free$row == free$col, ,
    drop = FALSE]
  if (!nrow(free)) return(NULL)

  for (k in seq_len(nrow(free))) {
    mat <- as.character(free$matrix[k]); i <- free$row[k]
    if (mat %in% "MANIFESTVAR") {
      # Only the types that have a residual at all; the rest are fixed to zero
      # and never reach here.
      if (i > nmanifest || !is.finite(mscale[i])) next
      set(free$parnumber[k], mscale[i] * .ctDataStartNoise)
      next
    }
    if (i > nlatent) next          # a carrier state, not a process
    if (mat %in% "T0VAR") { set(free$parnumber[k], latentsd[i]); next }
    if (is.na(ar[i])) next
    if (mat %in% "DRIFT") {
      value <- if (continuous) min(max(log(ar[i]) / dtmed[i], -3), -0.02) else ar[i]
      set(free$parnumber[k], value)
      next
    }
    if (mat %in% "DIFFUSION") {
      # Whatever makes the process's own stationary spread match the data:
      # sigma^2/(2|a|) in continuous time, sigma^2/(1-phi^2) in discrete.
      value <- if (continuous) {
        a <- min(max(log(ar[i]) / dtmed[i], -3), -0.02)
        latentsd[i] * sqrt(2 * abs(a))
      } else latentsd[i] * sqrt(max(1 - ar[i]^2, 1e-4))
      set(free$parnumber[k], value)
    }
  }

  if (!any(is.finite(target))) return(NULL)
  out <- rep(NA_real_, npar)
  for (k in seq_len(nrow(free))) {
    p <- free$parnumber[k]
    if (is.na(p) || p < 1 || p > npar || !is.finite(target[p])) next
    out[p] <- .ctDataStartInvert(free$transform[k], target[p])
  }
  if (!any(is.finite(out))) return(NULL)
  out
}
