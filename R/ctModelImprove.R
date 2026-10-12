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
# slow process driving a latent, a higher-order process, a model the user
# supplies -- are tried by refit.
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
#' adequate, a slow process driving a latent
#' (\code{\link{ctModelSlowComponent}}), a cycle such as a daily rhythm
#' (\code{\link{ctModelPeriodicComponent}}), a process of its own for an item
#' (\code{\link{ctModelUniqueProcess}}), raising the order of a latent process
#' (\code{\link{ctModelHigherOrder}}) and any models given in
#' \code{candidates} are tried the same way, each ranked by the adequacy
#' check's evidence for it: a slow process for a latent whose updates stay
#' correlated over time gaps (the \code{latent} family), a cycle at the period
#' the \code{periodic} family points to, a process of its own for an item
#' whose errors stay correlated with its own (the \code{residuals} family).
#'
#' The relaxations, in the language of the model:
#' \describe{
#'   \item{heterogeneity}{\code{indvarying = TRUE} for a parameter.}
#'   \item{covariates}{a covariate added as a time independent predictor, with
#'     effects on the parameters its tests single out.}
#'   \item{statedependence}{a term \code{d * eta_k * eta_m} in a latent's
#'     drift, written into \code{CINT} with \code{d} a \code{PARS} parameter:
#'     a drift coefficient that varies with a state.}
#'   \item{volatility}{a latent's innovation sd scaled by \code{exp(v *
#'     eta_m)}: innovations whose size depends on a state.}
#'   \item{fixed}{a cell fixed at zero freed: a drift cross effect, a
#'     diffusion or measurement error covariance, a cross loading.}
#'   \item{slow}{a slow process driving a latent through its drift, its effect
#'     started at one and its persistence at a tenth of the latent's.}
#'   \item{periodic}{a cycle of the period the adequacy check points to,
#'     driving every latent with innovations of its own, each subject's
#'     phase set from the clock at its first row and free to deviate; the
#'     phase predictors are added to the data.}
#'   \item{unique}{a process of its own, persistent, loading on one item.}
#'   \item{order}{a latent process made second order.}
#' }
#'
#' On a reduced-rank level (\code{poprank}) a heterogeneity result for a
#' parameter that already varies through the shared dimensions gives it a
#' variance of its own (\code{popunique}); without reduced rank, a
#' significant \code{dimension} test makes the parameters its direction
#' loads on vary together. Relaxations the score tests cannot tell apart (the same
#' \code{class}) are not refitted in the same round.
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
#'   data showed was subjects differing in variance. A search makes a decision
#'   per family it screens, so its error rate grows with the families: on two
#'   processes with three items each (60 subjects, about 48 rows each), where
#'   the screen also tests state dependence and volatility in each process,
#'   \code{alpha = 0.05} kept something in 4 of 20 searches with nothing
#'   missing (a volatility term twice, state dependence once, a random effect
#'   once), and \code{alpha = 0.01} is the stricter choice. With an omitted
#'   cross effect, a process of an item's own and a daily cycle in the same
#'   data, the search kept the cross effect in 19 of 20 searches, the cycle in
#'   16 and the item's process in all 20.
#'
#' The result is a model chosen by the data it is fitted to, and its
#' standard errors and tests do not account for that choice (MacCallum,
#' Roznowski & Necowitz, 1992). \code{holdout} sets aside a share of the
#' subjects: the search runs on the rest, and each change it keeps is then
#' tested on the held-out subjects alone, by refitting the final model without
#' it there; \code{confirmation} says which held, and the final model, refitted
#' to every subject, keeps only those. Every step, and every refit tried and
#' rejected, is in \code{history}.
#'
#' With many observations a model can be rejected by the checks for
#' departures too small to matter. With \code{close}, the search also stops
#' when every check that fails is a residual or within-row correlation family
#' whose size is not significantly above \code{close}
#' (\code{\link{ctFitAdequacy}}).
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
#' @param ntry Refits tried at once: a round refits the strongest \code{ntry}
#'   proposals. When those come from the adequacy check and none improves the
#'   fit, the next \code{ntry} are tried, and so on until one does or none is
#'   left, since the check does not change until the model does.
#' @param criterion \code{'LR'}, \code{'AIC'} or \code{'BIC'}.
#' @param adequacy Use \code{\link{ctFitAdequacy}} to decide when to stop.
#' @param prune After the search, try dropping each kept change in turn and
#'   drop it if the likelihood does not need it given the others: a change
#'   kept early can be made redundant by one kept later.
#' @param close Optional root mean square residual correlation regarded as a
#'   close enough fit to stop at; see Details.
#' @param holdout Share of subjects, between 0 and 0.5, held out of the search
#'   to confirm its changes on.
#' @param periods Periods for the adequacy check's \code{periodic} family, as
#'   for \code{\link{ctFitAdequacy}}; by default it scans.
#' @param form How covariates are screened, as for \code{\link{ctScoreTest}};
#'   \code{'threshold'} lets the search add a covariate as a step.
#' @param nsamples Generated datasets for the adequacy checks.
#' @param cores Cores for fitting and testing.
#' @param verbose Report progress.
#' @return An object of class \code{ctModelImprove}: \code{fit} (the final
#'   fit), \code{history} (every refit tried: the round, the relaxation, its
#'   screening p-value, the change in log likelihood and parameters, the
#'   criterion's verdict, and whether it was kept; a pruning row, \code{family}
#'   \code{'prune'}, says whether the change was dropped), \code{screens} (the score
#'   tests of each round), \code{adequacy} (the final checks),
#'   \code{confirmation} (with \code{holdout}: each kept change's likelihood
#'   ratio test on the held-out subjects), and \code{stopped} (why the search
#'   ended). With \code{holdout}, \code{fit} is the model with the confirmed
#'   changes refitted to every subject, and \code{selected} the search's own
#'   choice fitted to the subjects it searched on.
#' @seealso \code{vignette('score-tests')}.
#' @export
ctModelImprove <- function(fit, covariates = NULL,
  relax = c('heterogeneity', 'covariates', 'statedependence', 'volatility', 'fixed', 'slow',
    'unique', 'periodic', 'order'),
  candidates = NULL, alpha = 0.05, maxsteps = 10L, ntry = 3L,
  criterion = c('LR', 'AIC', 'BIC'), adequacy = TRUE, prune = TRUE,
  close = NULL, holdout = 0, form = 'auto', periods = NULL, nsamples = 100L, cores = 2L,
  verbose = TRUE) {
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
  selected <- NULL
  if (!is.numeric(holdout) || holdout < 0 || holdout > 0.5) stop("holdout is a ",
    "share of subjects between 0 and 0.5.", call. = FALSE)
  held <- NULL
  if (holdout > 0) {
    # A fixed stream, so a search is reproducible, leaving the caller's alone.
    if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
      oldseed <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
      on.exit(assign('.Random.seed', oldseed, envir = globalenv()), add = TRUE)
    } else on.exit(rm('.Random.seed', envir = globalenv()), add = TRUE)
    set.seed(20261011L)
    ids <- unique(data[[idname]])
    held <- sample(ids, max(2L, round(holdout * length(ids))))
    alldata <- data
    data <- data[!data[[idname]] %in% held, , drop = FALSE]
    say("Holding out ", length(held), " of ", length(ids), " subjects; refitting on the rest")
    fit <- .ctImproveRefit(fit, .ctFitBaseModel(fit), data, cores)
  }
  current <- fit
  startdata <- data
  kept <- list()
  extraargs <- list()
  history <- list(); screens <- list(); rejected <- character()
  ad <- NULL
  stopped <- 'maxsteps'
  for (step in seq_len(maxsteps)) {
    model <- .ctFitBaseModel(current)
    available <- if (is.null(covariates)) NULL else
      covariates[, c(idname, setdiff(names(covariates), c(idname, model$TIpredNames))),
        drop = FALSE]
    if (!is.null(available) && ncol(available) < 2L) available <- NULL
    tests <- intersect(relax, c('heterogeneity', 'covariates', 'statedependence', 'volatility',
      'fixed'))
    if (is.null(available)) tests <- setdiff(tests, 'covariates')
    say("Round ", step, ": screening")
    st <- if (length(tests)) suppressMessages(ctScoreTest(current, covariates = available,
      tests = tests, form = form, cores = cores, verbose = FALSE)) else NULL
    screens[[step]] <- st
    proposals <- .ctImproveProposals(st, model, alpha, current)
    proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
    structural <- FALSE
    if (!length(proposals)) {
      if (isTRUE(adequacy)) {
        say("Round ", step, ": nothing indicated; checking adequacy")
        ad <- ctFitAdequacy(current, covariates = available, periods = periods,
          nsamples = nsamples, alpha = alpha, close = if (is.null(close)) 0.05 else close, cores = cores)
        if (ad$adequate) { stopped <- 'adequate'; break }
        if (!is.null(close) && isTRUE(ad$close_fit)) { stopped <- 'close fit'; break }
      } else { stopped <- 'nothing indicated'; break }
      proposals <- .ctImproveStructural(model, relax, candidates, ad, current)
      proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
      structural <- TRUE
      if (!length(proposals)) { stopped <- 'not adequate, nothing left to try'; break }
    }
    # Ranked by each test's own asymptotic p-value: the bootstrap's floors at
    # 1 / (nboot + 1), and strong signals would tie there. Of relaxations the
    # tests cannot tell apart (one class), only the strongest is refitted this
    # round; if it fails, the next round's screen says whether another is
    # still indicated. Structural proposals come from a check that cannot
    # change until the model does, so when a batch of them fails the next
    # batch follows in the same round, until one improves or none is left;
    # stopping after the first left a daily cycle untried on the Fisher data
    # behind three item processes. Their likelihood ratio level is divided
    # by every structural proposal the round might try, not by the batch.
    proposals <- proposals[order(vapply(proposals, `[[`, 0, 'rank'))]
    cls <- vapply(proposals, function(p) if (is.null(p$class) || is.na(p$class)) NA_character_ else
      paste(p$family, p$class), '')
    proposals <- proposals[is.na(cls) | !duplicated(cls)]
    if (!structural) proposals <- utils::head(proposals, ntry)
    level <- alpha / length(proposals)
    tried <- list()
    ok <- logical()
    while (length(proposals)) {
      batch <- utils::head(proposals, ntry)
      proposals <- proposals[-seq_along(batch)]
      for (prop in batch) {
        say("Round ", step, ": refitting with ", prop$label)
        edited <- try(prop$edit(model, data, covariates), silent = TRUE)
        refit <- if (inherits(edited, 'try-error')) edited else
          try(.ctImproveRefit(current, edited$model, edited$data, cores, prop$start,
            edited$args), silent = TRUE)
        verdict <- .ctImproveCompare(current, refit, criterion, level)
        tried[[length(tried) + 1L]] <- list(prop = prop, refit = refit, edited = edited,
          verdict = verdict)
        history[[length(history) + 1L]] <- data.frame(round = step, family = prop$family,
          change = prop$label, screen_p = prop$p, dloglik = verdict$dloglik,
          dpar = verdict$dpar, criterion = verdict$value, kept = FALSE,
          note = verdict$note, stringsAsFactors = FALSE)
      }
      ok <- vapply(tried, function(t) isTRUE(t$verdict$accept), logical(1))
      if (any(ok)) break
    }
    rejected <- c(rejected, vapply(tried[!ok], function(t) t$prop$key, ''))
    if (!any(ok)) {
      if (structural) { stopped <- 'not adequate, no candidate improved'; break }
      next
    }
    best <- which(ok)[which.max(vapply(tried[ok], function(t) t$verdict$score, 0))]
    current <- tried[[best]]$refit
    data <- tried[[best]]$edited$data
    extraargs <- utils::modifyList(extraargs, as.list(tried[[best]]$edited$args))
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
      rebuilt <- try({
        r <- .ctImproveRebuild(.ctFitBaseModel(fit), startdata, kept[-k], covariates)
        .ctImproveRefit(current, r$model, r$data, cores, args = r$args)
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
        r <- .ctImproveRebuild(.ctFitBaseModel(fit), startdata, kept[-k], covariates)
        current <- rebuilt; data <- r$data; extraargs <- r$args; kept <- kept[-k]; ad <- NULL
      } else k <- k + 1L
    }
  }
  if (isTRUE(adequacy) && is.null(ad)) {
    model <- .ctFitBaseModel(current)
    available <- if (is.null(covariates)) NULL else
      covariates[, c(idname, setdiff(names(covariates), c(idname, model$TIpredNames))),
        drop = FALSE]
    if (!is.null(available) && ncol(available) < 2L) available <- NULL
    ad <- ctFitAdequacy(current, covariates = available, periods = periods,
      nsamples = nsamples, alpha = alpha, close = if (is.null(close)) 0.05 else close, cores = cores)
  }
  confirmation <- NULL
  if (!is.null(held) && length(kept)) {
    say("Confirming ", length(kept), " change(s) on the ", length(held), " held-out subjects")
    holddata <- alldata[alldata[[idname]] %in% held, , drop = FALSE]
    full <- .ctImproveRebuild(.ctFitBaseModel(fit), holddata, kept, covariates)
    fullfit <- try(.ctImproveRefit(current, full$model, full$data, cores, args = full$args),
      silent = TRUE)
    confirmation <- do.call(rbind, lapply(seq_along(kept), function(k) {
      without <- .ctImproveRebuild(.ctFitBaseModel(fit), holddata, kept[-k], covariates)
      w <- try(.ctImproveRefit(current, without$model, without$data, cores,
        args = without$args), silent = TRUE)
      verdict <- if (inherits(fullfit, 'try-error')) list(dloglik = NA_real_,
        dpar = NA_integer_, value = NA_real_, accept = NA, note = 'refit failed') else
        .ctImproveCompare(w, fullfit, 'LR', alpha)
      data.frame(change = kept[[k]]$label, dloglik = verdict$dloglik, dpar = verdict$dpar,
        p = verdict$value, confirmed = isTRUE(verdict$accept), note = verdict$note,
        stringsAsFactors = FALSE)
    }))
    # The final model keeps what the held-out subjects confirmed; the search's
    # own choice stays in `selected`.
    selected <- current
    confirmedkept <- kept[confirmation$confirmed]
    say("Refitting the final model, with the ", length(confirmedkept),
      " confirmed change(s), to every subject")
    final <- .ctImproveRebuild(.ctFitBaseModel(fit), alldata, confirmedkept, covariates)
    current <- .ctImproveRefit(current, final$model, final$data, cores, args = final$args)
  }
  history <- if (length(history)) do.call(rbind, history) else data.frame()
  structure(list(fit = current, history = history, screens = screens,
    adequacy = ad, confirmation = confirmation, stopped = stopped, start = fit,
    selected = selected,
    holdout = held), class = 'ctModelImprove')
}

# The base model with a sequence of kept changes applied, for some data: the
# model, the data the changes need (added covariates), and the ctFit
# arguments they change (popunique, poprank).
.ctImproveRebuild <- function(model, data, kept, covariates) {
  args <- list()
  for (p in kept) {
    e <- p$edit(model, data, covariates, args)
    model <- e$model; data <- e$data
    if (!is.null(e$args)) args <- utils::modifyList(args, e$args)
  }
  list(model = model, data = data, args = args)
}

# What the screen indicates, as edits: one per heterogeneity parameter, state
# term and fixed cell significant within its family, and one per covariate
# whose joint test is significant, adding it with effects on the parameters
# its tests single out.
.ctImproveProposals <- function(st, model, alpha, current = NULL) {
  out <- list()
  if (is.null(st)) return(out)
  # The ctFit arguments the current fit was made with, for edits that change
  # a reduced-rank level rather than the model.
  replay <- if (is.null(current)) list() else .ctFitReplayArgs(current)
  h <- st$heterogeneity
  if (!is.null(h)) for (r in which(h$p_adj < alpha)) {
    par <- h$parameter[r]
    varying <- par %in% model$pars$param[as.logical(model$pars$indvarying) %in% TRUE]
    out[[length(out) + 1L]] <- if (!varying) list(key = paste0('het:', par),
      family = 'heterogeneity', label = paste0(par, ' individually varying'), p = h$p[r],
      rank = stats::pnorm(h$z[r], lower.tail = FALSE), class = h$class[r],
      start = .ctImprovePopsdStart(model, par, h$sd[r]),
      edit = local({ par <- par; function(model, data, covs, args = list()) {
        model$pars$indvarying[model$pars$param %in% par] <- TRUE
        list(model = model, data = data)
      } })) else list(key = paste0('own:', par), family = 'heterogeneity',
      label = paste0(par, ' given a variance of its own'), p = h$p[r],
      rank = stats::pnorm(h$z[r], lower.tail = FALSE), class = h$class[r],
      edit = local({ par <- par; pu <- replay$popunique; function(model, data, covs, args = list()) {
        old <- if (!is.null(args$popunique)) args$popunique else pu
        if (isTRUE(old)) stop('every parameter already has a variance of its own')
        old <- if (is.null(old) || isFALSE(old)) character() else as.character(old)
        list(model = model, data = data, args = list(popunique = unique(c(old, par))))
      } }))
  }
  # The dimension test runs only where the level is not of reduced rank
  # (ctScoreTest), so a significant one makes the parameters its direction
  # loads on vary together.
  d <- st$dimension
  if (!is.null(d) && isTRUE(d$p < alpha)) {
    pars <- names(d$direction)[abs(d$direction) >= 0.5]
    pars <- setdiff(pars, model$pars$param[as.logical(model$pars$indvarying) %in% TRUE])
    if (length(pars) > 1L) out[[length(out) + 1L]] <- list(
      key = paste0('dim:', paste(sort(pars), collapse = ',')), family = 'dimension',
      label = paste0(paste(pars, collapse = ', '), ' individually varying together'),
      p = d$p, rank = d$p, edit = local({ pars <- pars; function(model, data, covs,
        args = list()) {
        model$pars$indvarying[model$pars$param %in% pars] <- TRUE
        list(model = model, data = data)
      } }))
  }
  for (cv in st$covariates) {
    if (!isTRUE(cv$joint$p < alpha)) next
    z <- cv$joint$covariate
    cut <- cv$joint$cut
    pt <- cv$parameters
    # Effects on spreads (popsd_, scale()) are reported, not added: a time
    # independent predictor acts on a parameter's value.
    pt <- pt[!grepl('^popsd_|^scale\\(', pt$parameter), , drop = FALSE]
    pars <- pt$parameter[pt$p_adj < alpha]
    if (!length(pars)) next
    # A column name the model can take: the sign and point as letters, since
    # 'age_above_-0.297' was not found in the data it had just been added to.
    zname <- if (is.null(cut)) z else paste0(z, '_above_',
      gsub('.', 'p', sub('^-', 'm', format(signif(cut, 3))), fixed = TRUE))
    out[[length(out) + 1L]] <- list(key = paste0('cov:', zname, ':', paste(sort(pars), collapse = ',')),
      family = 'covariates', label = paste0(if (is.null(cut)) z else paste0(z, ' above ',
        signif(cut, 3)), ' as a predictor of ', paste(pars, collapse = ', ')),
      p = cv$joint$p, rank = stats::pchisq(cv$joint$LM, cv$joint$df, lower.tail = FALSE),
      edit = local({ z <- z; pars <- pars; cut <- cut; zname <- zname
        function(model, data, covs, args = list()) {
        idname <- model$subjectIDname
        add <- unique(covs[, c(idname, z)])
        if (!is.null(cut)) {
          add[[zname]] <- as.numeric(add[[z]] > cut)
          add <- add[, c(idname, zname)]
          z <- zname
        }
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
        s$equation[r]), p = s$p[r], rank = s$p[r], class = s$class[r],
      edit = local({ i <- i; mono <- mono; name <- name; function(model, data, covs, args = list()) {
        list(model = .ctAddCintTerms(model, i, name, list(mono)), data = data)
      } }))
  }
  v <- st$volatility$terms
  if (!is.null(v)) for (r in which(v$p_adj < alpha)) {
    name <- paste0('vol_', v$equation[r], '_', v$by[r])
    out[[length(out) + 1L]] <- list(key = paste0('vol:', v$i[r], ':', v$m[r]),
      family = 'volatility', label = paste0('the innovation sd of ', v$equation[r],
        ' varying with ', v$by[r]), p = v$p[r], rank = v$p[r], class = v$class[r],
      edit = local({ i <- v$i[r]; m <- v$m[r]; name <- name
        function(model, data, covs, args = list()) {
          while (name %in% model$pars$param) name <- paste0(name, '_')
          list(model = .ctAddDiffusionTerms(model, i, name, m), data = data)
        } }))
  }
  f <- st$fixed$cells
  if (!is.null(f)) for (r in which(f$p_adj < alpha)) {
    cell <- f[r, ]
    out[[length(out) + 1L]] <- list(key = paste0('fixed:', cell$cell), family = 'fixed',
      label = paste0(cell$cell, ' freed'), p = cell$p, rank = cell$p, class = cell$class,
      edit = local({ cell <- cell; function(model, data, covs, args = list()) {
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
# checks fail: a slow process driving each latent, each latent made second
# order, and the user's own candidates. A slow process is ranked by the
# adequacy check's evidence for it, the latent's own standardised updates
# still positively correlated over some time gap (the smallest such p), and
# goes first where there is any; the rest rank alike.
.ctImproveStructural <- function(model, relax, candidates, ad, current = NULL) {
  out <- list()
  # A cycle the periodic family found, at the period its strongest positive
  # statistic names, driving every latent that has innovations of its own.
  if ('periodic' %in% relax && isTRUE(model$continuoustime)) {
    cyc <- .ctImprovePeriodicEvidence(ad)
    if (!is.null(cyc)) {
      name <- 'cycle'
      while (any(startsWith(c(model$latentNames, model$TIpredNames), paste0(name, '_'))))
        name <- paste0(name, 'x')
      drives <- .ctImproveNoisyLatents(model)
      out[[length(out) + 1L]] <- list(key = paste0('periodic:', signif(cyc$period, 3)),
        family = 'periodic', label = paste0('a cycle of period ', signif(cyc$period, 3)),
        p = cyc$p, rank = cyc$p,
        edit = local({ period <- cyc$period; name <- name; drives <- drives
          function(model, data, covs, args = list()) list(
            model = ctModelPeriodicComponent(model, period, intersect(drives,
              seq_len(model$n.latent)), name = name),
            data = ctPeriodicPhases(data, period, name, id = model$subjectIDname,
              time = model$timeName)) }))
    }
  }
  # An item whose errors stay correlated with its own over time, given a
  # process of its own; ranked by that evidence.
  if ('unique' %in% relax && isTRUE(model$continuoustime)) {
    for (item in model$manifestNames) {
      if (paste0('u_', item) %in% model$latentNames) next
      evidence <- .ctImproveOwnEvidence(ad, 'residuals', item)
      if (!is.finite(evidence)) next
      out[[length(out) + 1L]] <- list(key = paste0('unique:', item), family = 'unique',
        label = paste0(item, ' given a process of its own'), p = evidence, rank = evidence,
        start = .ctImproveUniqueStart(model, current, item),
        edit = local({ item <- item; function(model, data, covs, args = list())
          list(model = ctModelUniqueProcess(model, item), data = data) }))
    }
  }
  if ('slow' %in% relax && isTRUE(model$continuoustime)) {
    lat <- model$latentNames
    for (i in seq_along(lat)) {
      if (grepl('^slow_', lat[i]) || paste0('slow_', lat[i]) %in% lat ||
          lat[i] %in% paste0('d', lat)) next
      evidence <- .ctImproveSlowEvidence(ad, lat[i])
      out[[length(out) + 1L]] <- list(key = paste0('slow:', lat[i]), family = 'slow',
        label = paste0(lat[i], ' driven by a slow process'), p = evidence,
        rank = if (is.finite(evidence)) evidence else 1,
        start = .ctImproveSlowStart(model, current, i),
        edit = local({ i <- i; function(model, data, covs, args = list())
          list(model = ctModelSlowComponent(model, i), data = data) }))
    }
  }
  if ('order' %in% relax && isTRUE(model$continuoustime)) {
    for (i in seq_len(model$n.latent)) {
      if (paste0('d', model$latentNames[i]) %in% model$latentNames) next
      out[[length(out) + 1L]] <- list(key = paste0('order:', model$latentNames[i]),
        family = 'order', label = paste0(model$latentNames[i], ' made second order'),
        p = NA_real_, rank = 1, edit = local({ i <- i; function(model, data, covs,
          args = list()) list(model = ctModelHigherOrder(model, i), data = data) }))
    }
  }
  for (k in seq_along(candidates)) {
    out[[length(out) + 1L]] <- list(key = paste0('candidate:', k), family = 'candidate',
      label = paste0('candidate model ', k), p = NA_real_, rank = 1,
      edit = local({ k <- k; function(model, data, covs, args = list()) {
        list(model = candidates[[k]], data = data) } }))
  }
  out
}

# The adequacy check's evidence for a slow process in latent `name`: the
# smallest p of its updates' positive correlations with its own at other times,
# over every gap bin. NA when the check did not run the latent family.
.ctImproveSlowEvidence <- function(ad, name) .ctImproveOwnEvidence(ad, 'latent', name)

# The smallest p of a variable's positive correlations with itself at other
# times, over every gap bin, in one family of the adequacy check (`latent`
# for a latent's updates, `residuals` for an item's errors), times the number
# of gap bins it is the smallest of. Structural proposals are ranked against
# each other by this, and a cycle by its strongest period on the same terms
# (`.ctImprovePeriodicEvidence`). NA when the check did not run the family.
.ctImproveOwnEvidence <- function(ad, family, name) {
  s <- if (is.null(ad)) NULL else ad$statistics
  if (is.null(s) || !nrow(s)) return(NA_real_)
  all <- s$family == family & startsWith(s$statistic, paste0(name, '~', name, ' gap ')) &
    is.finite(s$z)
  own <- all & s$z > 0
  if (!any(own)) return(NA_real_)
  min(1, min(s$p[own]) * sum(all))
}

# The periodic family's evidence for a cycle: the period of its strongest
# positive statistic, and that statistic's p times the number of periods
# scanned -- on the terms an item's or a latent's own correlations are ranked
# by (`.ctImproveOwnEvidence`). The family's p, which the bootstrap floors at
# 1 / (nboot + 1), put a daily rhythm at z = 11 behind every item process on
# the Fisher data. NULL when the family did not run or found no positive one.
.ctImprovePeriodicEvidence <- function(ad) {
  if (is.null(ad) || is.null(ad$families)) return(NULL)
  f <- ad$families[ad$families$family == 'periodic', , drop = FALSE]
  s <- ad$statistics[ad$statistics$family == 'periodic' & is.finite(ad$statistics$z) &
    ad$statistics$z > 0, , drop = FALSE]
  if (!nrow(f) || !nrow(s)) return(NULL)
  best <- which.max(s$z)
  period <- as.numeric(sub('.* period ', '', s$statistic[best]))
  if (!is.finite(period)) return(NULL)
  scanned <- sum(ad$statistics$family == 'periodic' & is.finite(ad$statistics$z))
  list(p = min(1, s$p[best] * scanned), period = period)
}

# The latents with innovations of their own: not a cycle's clock, not the
# level of a higher order process.
.ctImproveNoisyLatents <- function(model) {
  p <- model$pars
  value <- suppressWarnings(as.numeric(p$value))
  quiet <- p$row[p$matrix %in% 'DIFFUSION' & p$row == p$col &
    (is.na(p$param) | !nzchar(p$param)) & !is.na(value) & value == 0]
  setdiff(seq_len(model$n.latent), quiet)
}

# Raw starts for a slow process added to latent `i`: its effect at one (zero
# is a stationary point of both it and the slow process), its persistence a
# tenth of the latent's own through the inverse of ctsem's default drift
# transform, -(1e-6 + 2 log1p_exp(-2 raw)). The latent's drift is read off the
# current fit; where it cannot be, -0.2, the pilot's start.
.ctImproveSlowStart <- function(model, current, i) {
  name <- model$latentNames[i]
  slow <- paste0('slow_', name)
  target <- -0.2
  cell <- model$pars[model$pars$matrix == 'DRIFT' & model$pars$row == i &
    model$pars$col == i, , drop = FALSE]
  fast <- NA_real_
  if (nrow(cell) == 1L) {
    if (is.na(cell$param)) fast <- as.numeric(cell$value) else if (!is.null(current)) {
      est <- setNames(.ctFitRawEstimate(current), .ctFitRawParNames(current))
      raw <- est[cell$param]
      fast <- tryCatch(as.numeric(eval(parse(text = cell$transform),
        list(param = raw, log1p_exp = log1p_exp))), error = function(e) NA_real_)
    }
  }
  if (is.finite(fast) && fast < 0) target <- fast / 10
  a <- max(-target - 1e-6, 1e-8) / 2
  c(setNames(-log(expm1(a)) / 2, paste0('drift_', slow)),
    setNames(1, paste0('drift_', name, '_', slow)))
}

# Raw starts for a process of its own for `item`. At raw zero ctsem's drift
# transform gives -1.39 an hour, a process that forgets within the hour, so on
# prompts hours apart it is indistinguishable from measurement error and the
# likelihood is flat in its persistence: on the Fisher data, processes for
# items whose errors correlate 0.2-0.3 with their own over a day gained 0.01
# to 2 log likelihood units from there, and none was kept. So it starts with
# a half-life of the median time between a subject's rows, and a stationary
# variance of a quarter of the item's error variance at the current fit, with
# the innovation and initial sds that implies, each through the inverse of its
# transform. NULL, and raw zero, when they cannot be worked out.
.ctImproveUniqueStart <- function(model, current, item) {
  if (is.null(current)) return(NULL)
  j <- match(item, model$manifestNames)
  m1 <- try(ctModelUniqueProcess(model, item), silent = TRUE)
  if (is.na(j) || inherits(m1, 'try-error')) return(NULL)
  value <- function(transform, raw) tryCatch(as.numeric(eval(parse(text = transform),
    list(param = raw, log1p_exp = log1p_exp))), error = function(e) NA_real_)
  inverse <- function(transform, target) tryCatch(stats::uniroot(
    function(r) value(transform, r) - target, c(-20, 20), tol = 1e-10)$root,
    error = function(e) NA_real_)
  d <- try(.ctFitLongData(current), silent = TRUE)
  if (inherits(d, 'try-error') || is.null(d)) return(NULL)
  gaps <- unlist(lapply(split(d[[model$timeName]], d[[model$subjectIDname]]),
    function(t) diff(sort(t))))
  gap <- stats::median(gaps[is.finite(gaps) & gaps > 0])
  cell <- model$pars[model$pars$matrix == 'MANIFESTVAR' & model$pars$row == j &
    model$pars$col == j, , drop = FALSE]
  sd <- NA_real_
  if (nrow(cell) == 1L) {
    if (is.na(cell$param)) sd <- as.numeric(cell$value) else {
      est <- setNames(.ctFitRawEstimate(current), .ctFitRawParNames(current))
      sd <- value(cell$transform, est[[cell$param]])
    }
  }
  if (!is.finite(gap) || !is.finite(sd) || sd <= 0) return(NULL)
  rate <- log(2) / gap
  v <- sd^2 / 4
  target <- c(DRIFT = -rate, DIFFUSION = sqrt(2 * rate * v), T0VAR = sqrt(v))
  u <- length(m1$latentNames)
  out <- numeric()
  for (mat in names(target)) {
    r <- m1$pars[m1$pars$matrix == mat & m1$pars$row == u & m1$pars$col == u, , drop = FALSE]
    if (nrow(r) != 1L || is.na(r$param)) next
    raw <- inverse(r$transform, target[[mat]])
    if (is.finite(raw)) out[r$param] <- raw
  }
  if (length(out)) out else NULL
}

# Refit with an edited model, prepared as the current fit was and warm started
# from its estimates, matched by name; new parameters start at zero, or where
# the proposal says (`start`, by raw name). A warm start can still fail -- a
# new random effect at an arbitrary scale can leave a Laplace unit's inner mode
# unconverged at the start, which ended two of three refits in testing -- so a
# refit that fails from there is tried once from ctsem's own start.
.ctImproveRefit <- function(current, model, data, cores, start = NULL, args = NULL) {
  extra <- args
  args <- .ctFitReplayArgs(current)
  if (length(extra)) args <- utils::modifyList(args, extra)
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
  # With held-out subjects the start and the search's result are fits to the
  # subjects searched on, and the final fit is to every subject: compare the
  # first two.
  end <- if (is.null(x$selected)) x$fit else x$selected
  cat("Log likelihood ", signif(x$start$estimate$loglik, 7), " -> ",
    signif(end$estimate$loglik, 7), " (", length(.ctFitRawEstimate(x$start)), " -> ",
    length(.ctFitRawEstimate(end)), " parameters)",
    if (!is.null(x$selected)) paste0(", on the ", .ctFitNsubjects(x$start),
      " subjects searched"), "\n\n", sep = "")
  if (nrow(x$history)) {
    h <- x$history
    h[c('screen_p', 'dloglik', 'criterion')] <- lapply(h[c('screen_p', 'dloglik', 'criterion')],
      signif, digits)
    h$round <- ifelse(is.na(h$round), 'prune', as.character(h$round))
    names(h)[names(h) == 'kept'] <- 'done'
    print(h[, c('round', 'change', 'screen_p', 'dloglik', 'dpar', 'done', 'note')],
      row.names = FALSE)
  }
  if (!is.null(x$confirmation)) {
    cat("\nConfirmed on ", length(x$holdout), " held-out subjects:\n", sep = "")
    cf <- x$confirmation
    cf[c('dloglik', 'p')] <- lapply(cf[c('dloglik', 'p')], signif, digits)
    print(cf[, c('change', 'dloglik', 'dpar', 'p', 'confirmed')], row.names = FALSE)
  }
  if (!is.null(x$adequacy)) { cat("\n"); print(x$adequacy) }
  invisible(x)
}
