# Improve a fitted model a step at a time.
#
# The search signal and the decision are kept apart. Score tests
# (`ctScoreTest`) say which relaxation the data point to and roughly how much
# it would add, at the cost of a few gradients; a refit with the relaxation
# decides whether it is kept, by its likelihood. Posterior predictive adequacy
# (`ctFitAdequacy`) is the stopping rule and the fallback: the score tests see
# only directions near the current model, and a model that is adequate by
# every check it was not tuned on is the aim. When the screen is clean and the
# checks still fail, structural changes the score tests cannot reach -- a
# higher-order process, a model the user supplies -- are tried by refit.
#
# Greedy, one change a step, retested after each: a missing relaxation echoes
# in the tests of related ones (a varying drift moves the diffusion's
# heterogeneity test; a missing moderation moves a term in the other equation),
# and those echoes go once the source is in the model.

#' Improve a model step by step, from score tests and predictive checks
#'
#' Starting from a fitted model, repeatedly: screens it with
#' \code{\link{ctScoreTest}}; refits it with the most strongly indicated
#' relaxations, one at a time, warm started from the current estimates; keeps
#' the best of them if the refit's likelihood confirms it; and stops when no
#' relaxation is indicated and the model is adequate by
#' \code{\link{ctFitAdequacy}}. When nothing is indicated but the model is not
#' adequate, raising the order of a latent process
#' (\code{\link{ctModelHigherOrder}}) and any models given in
#' \code{candidates} are tried the same way.
#'
#' The relaxations, in the language of the model:
#' \describe{
#'   \item{heterogeneity}{\code{indvarying = TRUE} for a parameter.}
#'   \item{covariates}{a covariate added as a time independent predictor, with
#'     effects on the parameters its tests single out.}
#'   \item{statedependence}{a term \code{d * eta_k * eta_m} in a latent's
#'     drift, written into \code{CINT} with \code{d} a \code{PARS} parameter:
#'     a drift coefficient that varies with a state.}
#'   \item{fixed}{a cell fixed at zero freed: a drift cross effect, a
#'     diffusion or measurement error covariance, a cross loading.}
#'   \item{order}{a latent process made second order.}
#' }
#'
#' When the search stops, each change kept is tried out again, by refitting
#' without it (\code{prune}), and dropped if the others make it redundant.
#'
#' A relaxation is kept when its refit improves on the current model by
#' \code{criterion}: a likelihood ratio test at \code{alpha} divided by the
#' number of refits tried at that step (\code{'LR'}), or a lower AIC or BIC
#' (subjects as the sample size). The likelihood ratio test of a new random
#' effect against the chi-square is conservative, as its variance starts on a
#' boundary.
#'
#' @section Error rates: The screen and the refit read the same data, so a
#'   refit confirms most of what the screen proposes, and the search's error
#'   rate is close to \code{alpha} per decision. Measured on simulated data (one
#'   latent, 80 subjects, 20 irregular rows, 40 searches each): with nothing
#'   missing, \code{alpha = 0.05} kept something in 5 per cent of searches and
#'   \code{alpha = 0.01} in none; an omitted covariate shifting a random
#'   intercept by half its sd was added in 75 and 55 per cent; an unrelated
#'   covariate never. \code{criterion = 'BIC'} kept something in 8 per cent,
#'   so at that size it is not the stricter choice. A drift that varied over
#'   subjects was never added as such at that size; a random diffusion or
#'   measurement error variance stood in for it in 10 per cent, since what the
#'   data showed was subjects differing in variance.
#'
#' The result is a model chosen by the data it is fitted to, and its
#' standard errors and tests do not account for that choice. Treat the
#' selected structure as a hypothesis to check on other data, or hold out
#' subjects for that from the start. Every step, and every refit tried and
#' rejected, is in \code{history}.
#'
#' @param fit A fit from \code{\link{ctFit}} with \code{backend = 'julia'}.
#' @param covariates Optional subject-level covariates the model may add, as
#'   for \code{\link{ctScoreTest}}.
#' @param relax Which relaxations may be made.
#' @param candidates Optional list of alternative models (from
#'   \code{\link{ctModel}}) to try when the screen is clean but the checks fail.
#' @param alpha Level for the screen (Benjamini-Hochberg within each family),
#'   for the likelihood ratio decision, and for adequacy.
#' @param maxsteps Most refit rounds.
#' @param ntry Most refits tried per round.
#' @param criterion \code{'LR'}, \code{'AIC'} or \code{'BIC'}.
#' @param adequacy Use \code{\link{ctFitAdequacy}} to decide when to stop.
#' @param prune After the search, try dropping each kept change in turn and
#'   drop it if the likelihood does not need it given the others: a change
#'   kept early can be made redundant by one kept later.
#' @param nsamples Generated datasets for the adequacy checks.
#' @param cores Cores for fitting and testing.
#' @param verbose Report progress.
#' @return An object of class \code{ctModelImprove}: \code{fit} (the final
#'   fit), \code{history} (every refit tried: the round, the relaxation, its
#'   screening p-value, the change in log likelihood and parameters, the
#'   criterion's verdict, and whether it was kept; a pruning row, \code{family}
#'   \code{'prune'}, says whether the change was dropped), \code{screens} (the score
#'   tests of each round), \code{adequacy} (the final checks), and
#'   \code{stopped} (why the search ended).
#' @seealso \code{vignette('score-tests')}.
#' @export
ctModelImprove <- function(fit, covariates = NULL,
  relax = c('heterogeneity', 'covariates', 'statedependence', 'fixed', 'order'),
  candidates = NULL, alpha = 0.05, maxsteps = 10L, ntry = 3L,
  criterion = c('LR', 'AIC', 'BIC'), adequacy = TRUE, prune = TRUE,
  nsamples = 100L, cores = 2L, verbose = TRUE) {
  if (!.ctFitIsJulia(fit)) stop("ctModelImprove() needs a fit from ",
    "ctFit(backend = 'julia').", call. = FALSE)
  if (is.null(.ctFitBaseModel(fit))) stop("This fit does not carry the model it ",
    "was built from; refit it with this version of ctsem.", call. = FALSE)
  .ctBackendResolveThreads(cores)
  relax <- match.arg(relax, several.ok = TRUE)
  criterion <- match.arg(criterion)
  say <- function(...) if (isTRUE(verbose)) message(...)
  covariates <- if (is.null(covariates)) NULL else as.data.frame(covariates)
  idname <- .ctFitModelObject(fit)$subjectIDname
  if (!is.null(covariates) && !any(c(idname, 'id') %in% names(covariates)))
    stop("covariates needs a '", idname, "' column.", call. = FALSE)
  if (!is.null(covariates) && !idname %in% names(covariates))
    names(covariates)[names(covariates) == 'id'] <- idname
  data <- as.data.frame(.ctFitLongData(fit))
  current <- fit
  startdata <- data
  kept <- list()
  history <- list(); screens <- list(); rejected <- character()
  ad <- NULL
  stopped <- 'maxsteps'
  for (step in seq_len(maxsteps)) {
    model <- .ctFitBaseModel(current)
    available <- if (is.null(covariates)) NULL else
      covariates[, c(idname, setdiff(names(covariates), c(idname, model$TIpredNames))),
        drop = FALSE]
    if (!is.null(available) && ncol(available) < 2L) available <- NULL
    tests <- intersect(relax, c('heterogeneity', 'covariates', 'statedependence', 'fixed'))
    if (is.null(available)) tests <- setdiff(tests, 'covariates')
    say("Round ", step, ": screening")
    st <- if (length(tests)) suppressMessages(ctScoreTest(current, covariates = available,
      tests = tests, cores = cores, verbose = FALSE)) else NULL
    screens[[step]] <- st
    proposals <- .ctImproveProposals(st, model, alpha)
    proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
    structural <- FALSE
    if (!length(proposals)) {
      if (isTRUE(adequacy)) {
        say("Round ", step, ": nothing indicated; checking adequacy")
        ad <- ctFitAdequacy(current, covariates = available, nsamples = nsamples,
          alpha = alpha, cores = cores)
        if (ad$adequate) { stopped <- 'adequate'; break }
      } else { stopped <- 'nothing indicated'; break }
      proposals <- .ctImproveStructural(model, relax, candidates, ad)
      proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
      structural <- TRUE
      if (!length(proposals)) { stopped <- 'not adequate, nothing left to try'; break }
    }
    # Ranked by each test's own asymptotic p-value: the bootstrap's floors at
    # 1 / (nboot + 1), and strong signals would tie there.
    proposals <- utils::head(proposals[order(vapply(proposals, `[[`, 0, 'rank'))], ntry)
    tried <- list()
    for (prop in proposals) {
      say("Round ", step, ": refitting with ", prop$label)
      edited <- try(prop$edit(model, data, covariates), silent = TRUE)
      refit <- if (inherits(edited, 'try-error')) edited else
        try(.ctImproveRefit(current, edited$model, edited$data, cores, prop$start),
          silent = TRUE)
      verdict <- .ctImproveCompare(current, refit, criterion, alpha / length(proposals))
      tried[[length(tried) + 1L]] <- list(prop = prop, refit = refit, edited = edited,
        verdict = verdict)
      history[[length(history) + 1L]] <- data.frame(round = step, family = prop$family,
        change = prop$label, screen_p = prop$p, dloglik = verdict$dloglik,
        dpar = verdict$dpar, criterion = verdict$value, kept = FALSE,
        note = verdict$note, stringsAsFactors = FALSE)
    }
    ok <- vapply(tried, function(t) isTRUE(t$verdict$accept), logical(1))
    rejected <- c(rejected, vapply(tried[!ok], function(t) t$prop$key, ''))
    if (!any(ok)) {
      if (structural) { stopped <- 'not adequate, no candidate improved'; break }
      next
    }
    best <- which(ok)[which.max(vapply(tried[ok], function(t) t$verdict$score, 0))]
    current <- tried[[best]]$refit
    data <- tried[[best]]$edited$data
    kept[[length(kept) + 1L]] <- tried[[best]]$prop
    history[[length(history) - length(tried) + best]]$kept <- TRUE
    ad <- NULL
    say("Round ", step, ": kept ", tried[[best]]$prop$label)
  }
  # Backward: a change kept early can be made redundant by one kept later --
  # a random effect on a parameter whose variation a covariate then explains.
  # Each kept change in turn, the model is rebuilt from the start with all the
  # others and refitted; if the fuller model is not better by the criterion,
  # the change goes.
  if (isTRUE(prune) && length(kept) > 1L) {
    k <- 1L
    while (k <= length(kept) && length(kept) > 1L) {
      say("Pruning: refitting without ", kept[[k]]$label)
      m <- .ctFitBaseModel(fit); d <- startdata
      rebuilt <- try({
        for (p in kept[-k]) { e <- p$edit(m, d, covariates); m <- e$model; d <- e$data }
        .ctImproveRefit(current, m, d, cores)
      }, silent = TRUE)
      verdict <- if (inherits(rebuilt, 'try-error')) list(accept = TRUE,
        dloglik = NA_real_, dpar = NA_integer_, value = NA_real_,
        note = 'refit without it failed; kept') else
        .ctImproveCompare(rebuilt, current, criterion, alpha)
      drop <- !isTRUE(verdict$accept)
      history[[length(history) + 1L]] <- data.frame(round = NA_integer_,
        family = 'prune', change = paste0('drop ', kept[[k]]$label), screen_p = NA_real_,
        dloglik = if (is.finite(verdict$dloglik)) -verdict$dloglik else NA_real_,
        dpar = if (is.finite(verdict$dpar)) -verdict$dpar else NA_integer_,
        criterion = verdict$value, kept = drop, note = verdict$note,
        stringsAsFactors = FALSE)
      if (drop) {
        say("Pruning: dropped ", kept[[k]]$label)
        current <- rebuilt; data <- d; kept <- kept[-k]; ad <- NULL
      } else k <- k + 1L
    }
  }
  if (isTRUE(adequacy) && is.null(ad)) {
    model <- .ctFitBaseModel(current)
    available <- if (is.null(covariates)) NULL else
      covariates[, c(idname, setdiff(names(covariates), c(idname, model$TIpredNames))),
        drop = FALSE]
    if (!is.null(available) && ncol(available) < 2L) available <- NULL
    ad <- ctFitAdequacy(current, covariates = available, nsamples = nsamples,
      alpha = alpha, cores = cores)
  }
  history <- if (length(history)) do.call(rbind, history) else data.frame()
  structure(list(fit = current, history = history, screens = screens,
    adequacy = ad, stopped = stopped, start = fit), class = 'ctModelImprove')
}

# What the screen indicates, as edits: one per heterogeneity parameter, state
# term and fixed cell significant within its family, and one per covariate
# whose joint test is significant, adding it with effects on the parameters
# its tests single out.
.ctImproveProposals <- function(st, model, alpha) {
  out <- list()
  if (is.null(st)) return(out)
  h <- st$heterogeneity
  if (!is.null(h)) for (r in which(h$p_adj < alpha)) {
    par <- h$parameter[r]
    out[[length(out) + 1L]] <- list(key = paste0('het:', par), family = 'heterogeneity',
      label = paste0(par, ' individually varying'), p = h$p[r],
      rank = stats::pnorm(h$z[r], lower.tail = FALSE),
      start = .ctImprovePopsdStart(model, par, h$sd[r]),
      edit = local({ par <- par; function(model, data, covs) {
        model$pars$indvarying[model$pars$param %in% par] <- TRUE
        list(model = model, data = data)
      } }))
  }
  for (cv in st$covariates) {
    if (!isTRUE(cv$joint$p < alpha)) next
    z <- cv$joint$covariate
    pt <- cv$parameters
    pt <- pt[!grepl('^popsd_', pt$parameter), , drop = FALSE]
    pars <- pt$parameter[pt$p_adj < alpha]
    if (!length(pars)) next
    out[[length(out) + 1L]] <- list(key = paste0('cov:', z, ':', paste(sort(pars), collapse = ',')),
      family = 'covariates', label = paste0(z, ' as a predictor of ', paste(pars, collapse = ', ')),
      p = cv$joint$p, rank = stats::pchisq(cv$joint$LM, cv$joint$df, lower.tail = FALSE), edit = local({ z <- z; pars <- pars; function(model, data, covs) {
        idname <- model$subjectIDname
        add <- unique(covs[, c(idname, z)])
        data <- merge(data[, setdiff(names(data), z), drop = FALSE], add, by = idname,
          all.x = TRUE, sort = FALSE)
        data <- data[order(data[[idname]], data[[model$timeName]]), , drop = FALSE]
        model$TIpredNames <- c(model$TIpredNames, z)
        model$n.TIpred <- length(model$TIpredNames)
        model$pars[[paste0(z, '_effect')]] <- ifelse(model$pars$param %in% pars, 'TRUE', 'FALSE')
        list(model = model, data = data)
      } }))
  }
  s <- st$statedependence$terms
  if (!is.null(s)) for (r in which(s$p_adj < alpha)) {
    mono <- as.integer(strsplit(s$monomial[r], ',')[[1]])
    i <- s$i[r]
    name <- paste0('sdep_', s$equation[r], '_', gsub('*', '_', s$term[r], fixed = TRUE))
    out[[length(out) + 1L]] <- list(key = paste0('sdep:', i, ':', s$monomial[r]),
      family = 'statedependence', label = paste0(s$term[r], ' in the drift of ',
        s$equation[r]), p = s$p[r], rank = s$p[r],
      edit = local({ i <- i; mono <- mono; name <- name; function(model, data, covs) {
        list(model = .ctAddCintTerms(model, i, name, list(mono)), data = data)
      } }))
  }
  f <- st$fixed$cells
  if (!is.null(f)) for (r in which(f$p_adj < alpha)) {
    cell <- f[r, ]
    out[[length(out) + 1L]] <- list(key = paste0('fixed:', cell$cell), family = 'fixed',
      label = paste0(cell$cell, ' freed'), p = cell$p, rank = cell$p,
      edit = local({ cell <- cell; function(model, data, covs) {
        k <- which(model$pars$matrix == cell$matrix & model$pars$row == cell$row &
          model$pars$col == cell$col)
        name <- cell$suggested
        while (name %in% model$pars$param) name <- paste0(name, '_')
        model$pars$param[k] <- name; model$pars$value[k] <- NA
        model$pars$transform[k] <- cell$transform; model$pars$indvarying[k] <- FALSE
        for (e in grep('_effect$', names(model$pars), value = TRUE)) model$pars[[e]][k] <- 'FALSE'
        list(model = model, data = data)
      } }))
  }
  out
}

# Changes the score tests cannot reach, for when the screen is clean and the
# checks fail: each latent made second order, and the user's own candidates.
.ctImproveStructural <- function(model, relax, candidates, ad) {
  out <- list()
  if ('order' %in% relax && isTRUE(model$continuoustime)) {
    for (i in seq_len(model$n.latent)) {
      if (paste0('d', model$latentNames[i]) %in% model$latentNames) next
      out[[length(out) + 1L]] <- list(key = paste0('order:', model$latentNames[i]),
        family = 'order', label = paste0(model$latentNames[i], ' made second order'),
        p = NA_real_, rank = 1, edit = local({ i <- i; function(model, data, covs)
          list(model = ctModelHigherOrder(model, i), data = data) }))
    }
  }
  for (k in seq_along(candidates)) {
    out[[length(out) + 1L]] <- list(key = paste0('candidate:', k), family = 'candidate',
      label = paste0('candidate model ', k), p = NA_real_, rank = 1,
      edit = local({ k <- k; function(model, data, covs) {
        list(model = candidates[[k]], data = data) } }))
  }
  out
}

# Refit with an edited model, prepared as the current fit was and warm started
# from its estimates, matched by name; new parameters start at zero, or where
# the proposal says (`start`, by raw name). A warm start can still fail -- a
# new random effect at an arbitrary scale can leave a Laplace unit's inner mode
# unconverged at the start, which ended two of three refits in testing -- so a
# refit that fails from there is tried once from ctsem's own start.
.ctImproveRefit <- function(current, model, data, cores, start = NULL) {
  args <- .ctFitReplayArgs(current)
  args$model <- model
  args$datalong <- data
  args$verbose <- 0
  args$fit <- FALSE
  args$inits <- NULL
  spec <- suppressMessages(suppressWarnings(do.call(ctFit, args)))
  names2 <- as.character(.ctBackendRawParameterNames(list(model_spec = spec),
    .ctBackendNpar(spec)))
  est <- setNames(.ctFitRawEstimate(current), .ctFitRawParNames(current))
  init <- setNames(numeric(length(names2)), names2)
  shared <- intersect(names(est), names2)
  init[shared] <- est[shared]
  start <- start[names(start) %in% names2]
  init[names(start)] <- start
  args$fit <- TRUE
  args$inits <- unname(init)
  out <- try(suppressMessages(suppressWarnings(do.call(ctFit, args))), silent = TRUE)
  if (!inherits(out, 'try-error')) return(out)
  args$inits <- NULL
  suppressMessages(suppressWarnings(do.call(ctFit, args)))
}

# The raw start for a new random effect's population sd: the heterogeneity
# test's one-step sd, through the inverse of the model's sd transform
# (`log1p_exp(2 * raw - 1) * sdscale`), floored so the optimiser does not start
# on the boundary where the sd's gradient vanishes.
.ctImprovePopsdStart <- function(model, par, sd) {
  sdscale <- suppressWarnings(as.numeric(model$pars$sdscale[match(par, model$pars$param)]))
  if (!is.finite(sdscale) || sdscale <= 0) sdscale <- 1
  s <- max(if (is.finite(sd)) sd else 0, 0.05) / sdscale
  setNames((log(expm1(s)) + 1) / 2, paste0('popsd_', par))
}

.ctImproveCompare <- function(current, refit, criterion, alpha) {
  bad <- function(note) list(accept = FALSE, dloglik = NA_real_, dpar = NA_integer_,
    value = NA_real_, score = -Inf, note = note)
  if (inherits(refit, 'try-error')) return(bad(paste('failed:',
    trimws(conditionMessage(attr(refit, 'condition'))))))
  ll0 <- current$estimate$loglik; ll1 <- refit$estimate$loglik
  k0 <- length(.ctFitRawEstimate(current)); k1 <- length(.ctFitRawEstimate(refit))
  dll <- ll1 - ll0; dk <- k1 - k0
  if (!is.finite(dll)) return(bad('log likelihood not finite'))
  n <- .ctFitNsubjects(current)
  if (criterion == 'LR' && dk > 0) {
    p <- stats::pchisq(2 * dll, dk, lower.tail = FALSE)
    return(list(accept = p < alpha, dloglik = dll, dpar = dk, value = p,
      score = -log(p), note = sprintf('LR p = %.3g against %.3g', p, alpha)))
  }
  penalty <- if (criterion == 'BIC') log(n) else 2
  delta <- -2 * dll + penalty * dk
  list(accept = delta < 0, dloglik = dll, dpar = dk, value = delta, score = -delta,
    note = sprintf('%s change %.3g', if (criterion == 'BIC') 'BIC' else 'AIC', delta))
}

#' @export
print.ctModelImprove <- function(x, digits = 3, ...) {
  cat("Model improvement: stopped because ", x$stopped, "\n", sep = "")
  cat("Log likelihood ", signif(x$start$estimate$loglik, 7), " -> ",
    signif(x$fit$estimate$loglik, 7), " (", length(.ctFitRawEstimate(x$start)), " -> ",
    length(.ctFitRawEstimate(x$fit)), " parameters)\n\n", sep = "")
  if (nrow(x$history)) {
    h <- x$history
    h[c('screen_p', 'dloglik', 'criterion')] <- lapply(h[c('screen_p', 'dloglik', 'criterion')],
      signif, digits)
    h$round <- ifelse(is.na(h$round), 'prune', as.character(h$round))
    names(h)[names(h) == 'kept'] <- 'done'
    print(h[, c('round', 'change', 'screen_p', 'dloglik', 'dpar', 'done', 'note')],
      row.names = FALSE)
  }
  if (!is.null(x$adequacy)) { cat("\n"); print(x$adequacy) }
  invisible(x)
}
