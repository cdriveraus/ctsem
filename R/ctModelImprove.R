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
#'     diffusion covariance, a cross loading.}
#'   \item{slow}{a slow process driving a latent through its drift, its effect
#'     started at one and its persistence at a tenth of the latent's.}
#'   \item{periodic}{a cycle of the period the adequacy check points to,
#'     driving every latent with innovations of its own, each subject's
#'     phase set from the clock at its first row and free to deviate; the
#'     phase predictors are added to the data.}
#'   \item{unique}{a process of its own, persistent, loading on one item;
#'     and, where the \code{pairs} family finds two items' errors correlated
#'     within an occasion, a process the two share
#'     (\code{ctModelUniqueProcess(shared = TRUE)}): an error is meant to be
#'     random, and co-movement within an occasion that the latents' innovation
#'     covariances do not explain is read as a process faster than the
#'     occasions. Like the other structural changes it is tried once the
#'     score tests indicate nothing, ranked against them by the check's
#'     evidence, so that it does not absorb misfit that belongs to a cycle or
#'     to one item's own process.}
#'   \item{errorcov}{not by default: a measurement error covariance the score
#'     test flags is freed, for designs with a reason for correlated error,
#'     such as a method the items share, and no shared process is proposed.}
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
#' Each random effect the search adds brings its correlations with the others,
#' and with many effects for few subjects the data no longer determine them.
#' When a kept change leaves the fit's curvature weak or negative along the
#' population covariance (\code{fit$identifiability}), the covariance is
#' refitted at one rank lower (\code{poprank}), starting from its leading
#' dimensions, and the lower rank kept unless the likelihood ratio test needs
#' the dropped dimension; it steps down until the fit is identified. Later
#' random effects then load on the shared dimensions, and later heterogeneity
#' results give own variances. Single-level models only.
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
#' @param ntry Proposals from the score tests refitted in a round, at most:
#'   one at a time, the strongest first, keeping the first that improves the
#'   fit. Proposals from the adequacy check are tried the same way until one
#'   improves or none is left, since the check does not change until the model
#'   does. Pruning and a held-out confirmation refit this many side by side.
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
#' @param projection How the adequacy checks remove the estimation, as for
#'   \code{\link{ctFitAdequacy}}.
#' @param cores Cores for fitting and testing. With two or more, once a
#'   refit has taken longer than five minutes, the refits of pruning and of a
#'   held-out confirmation run side by side in worker processes, \code{ntry}
#'   at a time, as sampling chains do, each with an equal share of
#'   \code{cores}: they are independent, one fit uses its cores only partly,
#'   and a worker takes about half a minute to start.
#'   \code{options(ctsem.improve.processes = FALSE)} keeps them in this
#'   session.
#' @param verbose Report progress.
#' @return An object of class \code{ctModelImprove}: \code{fit} (the final
#'   fit), \code{history} (every refit tried: the round, the relaxation, its
#'   screening p-value, the change in log likelihood and parameters, the
#'   criterion's verdict, and whether it was kept; a pruning row, \code{family}
#'   \code{'prune'}, says whether the change was dropped, and a \code{'rank'} row
#'   whether the lower rank was kept), \code{screens} (the score
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
  close = NULL, holdout = 0, form = 'auto', periods = NULL, nsamples = 100L,
  projection = 'auto', cores = getOption('mc.cores', 2L),
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
  slowest <- 0
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
    proposals <- .ctImproveProposals(st, model, alpha, current, relax)
    proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
    structural <- FALSE
    if (!length(proposals)) {
      if (isTRUE(adequacy)) {
        say("Round ", step, ": nothing indicated; checking adequacy")
        ad <- ctFitAdequacy(current, covariates = available, periods = periods,
          nsamples = nsamples, alpha = alpha, close = if (is.null(close)) 0.05 else close,
          projection = projection, cores = cores)
        if (ad$adequate) { stopped <- 'adequate'; break }
        if (!is.null(close) && isTRUE(ad$close_fit)) { stopped <- 'close fit'; break }
      } else { stopped <- 'nothing indicated'; break }
      proposals <- .ctImproveStructural(model, relax, candidates, ad, current, alpha)
      proposals <- proposals[!vapply(proposals, `[[`, '', 'key') %in% rejected]
      structural <- TRUE
      if (!length(proposals)) { stopped <- 'not adequate, nothing left to try'; break }
      tries <- length(proposals)
      say("Round ", step, ": scoring ", tries, " structural proposal(s)")
      scored <- .ctImproveScoreStructural(current, proposals, data, covariates, cores)
      proposals <- scored$proposals
      for (sk in scored$skipped) {
        history[[length(history) + 1L]] <- data.frame(round = step, family = sk$family,
          change = sk$label, screen_p = sk$p, dloglik = NA_real_, dpar = NA_integer_,
          criterion = NA_real_, kept = FALSE, note = sk$note, stringsAsFactors = FALSE)
        rejected <- c(rejected, sk$key)
      }
      if (!length(proposals)) { stopped <- 'not adequate, no proposal scored'; break }
    }
    # Ranked by each test's own asymptotic p-value: the bootstrap's floors at
    # 1 / (nboot + 1), and strong signals would tie there. Of relaxations the
    # tests cannot tell apart (one class), only the strongest is refitted this
    # round; if it fails, the next round's screen says whether another is
    # still indicated. Structural proposals come from a check that cannot
    # change until the model does, so when one fails the next follows in the
    # same round, until one improves or none is left; stopping after the
    # first left a daily cycle untried on the Fisher data behind three item
    # processes. Their likelihood ratio level is divided by every structural
    # proposal the round might try.
    #
    # One refit at a time, the strongest first, and the first that improves
    # is kept. Refitting three side by side and keeping the best spent two
    # refits a round for little: on the Fisher search every structural round
    # that kept a change kept the first tried, and in the score-test rounds
    # every one of the three passed by tens of log likelihood units, so
    # keeping the first reorders additions the next rounds make anyway. Each
    # refit then has all the cores.
    proposals <- proposals[order(vapply(proposals, `[[`, 0, 'rank'))]
    cls <- vapply(proposals, function(p) if (is.null(p$class) || is.na(p$class)) NA_character_ else
      paste(p$family, p$class), '')
    proposals <- proposals[is.na(cls) | !duplicated(cls)]
    if (!structural) proposals <- utils::head(proposals, ntry)
    level <- alpha / if (structural) tries else length(proposals)
    tried <- list()
    ok <- logical()
    while (length(proposals)) {
      batch <- utils::head(proposals, 1L)
      proposals <- proposals[-seq_along(batch)]
      edits <- lapply(batch, function(prop) try(prop$edit(model, data, covariates),
        silent = TRUE))
      for (prop in batch) say("Round ", step, ": refitting with ", prop$label)
      refits <- .ctImproveRefits(current, edits, lapply(batch, `[[`, 'start'), cores,
        slowest)
      slowest <- max(slowest, attr(refits, 'seconds'))
      for (b in seq_along(batch)) {
        prop <- batch[[b]]
        edited <- edits[[b]]
        refit <- refits[[b]]
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
    # A kept random effect can leave the population covariance more than the
    # data determine: then the rank comes down until the likelihood needs it.
    ranked <- .ctImproveRankStep(.ctFitBaseModel(fit), current, kept, startdata, covariates,
      criterion, alpha, cores, say, label = paste0("Round ", step, ": "), round = step)
    history <- c(history, ranked$history)
    if (ranked$changed) {
      current <- ranked$current; kept <- ranked$kept; data <- ranked$data
      extraargs <- ranked$args
    }
  }
  # Backward: a change kept early can be made redundant by one kept later --
  # a random effect on a parameter whose variation a covariate then explains.
  # Each kept change in turn, the model is rebuilt from the start with all the
  # others and refitted; if the fuller model is not better by the criterion,
  # the change goes.
  if (isTRUE(prune) && sum(vapply(kept, .ctImprovePrunable, logical(1))) > 1L) {
    pruned <- .ctImprovePrune(.ctFitBaseModel(fit), current, kept, startdata, covariates,
      criterion, alpha, cores, ntry, slowest, say)
    history <- c(history, pruned$history)
    slowest <- pruned$slowest
    if (length(pruned$kept) < length(kept)) {
      current <- pruned$current; data <- pruned$data; extraargs <- pruned$args; ad <- NULL
    }
    kept <- pruned$kept
  }
  if (isTRUE(adequacy) && is.null(ad)) {
    model <- .ctFitBaseModel(current)
    available <- if (is.null(covariates)) NULL else
      covariates[, c(idname, setdiff(names(covariates), c(idname, model$TIpredNames))),
        drop = FALSE]
    if (!is.null(available) && ncol(available) < 2L) available <- NULL
    ad <- ctFitAdequacy(current, covariates = available, periods = periods,
      nsamples = nsamples, alpha = alpha, close = if (is.null(close)) 0.05 else close,
      projection = projection, cores = cores)
  }
  confirmation <- NULL
  if (!is.null(held) && length(kept)) {
    say("Confirming ", length(kept), " change(s) on the ", length(held), " held-out subjects")
    holddata <- alldata[alldata[[idname]] %in% held, , drop = FALSE]
    # The model with every change and the model without each one, refitted
    # to the held-out subjects: independent refits, side by side when refits
    # go to workers, `ntry` at a time.
    rebuilt <- lapply(c(list(seq_along(kept)), lapply(seq_along(kept), function(k)
      seq_along(kept)[-k])), function(idx) try(.ctImproveRebuild(.ctFitBaseModel(fit),
        holddata, kept[idx], covariates), silent = TRUE))
    refits <- .ctImproveRefits(current, rebuilt, vector('list', length(rebuilt)), cores,
      slowest, width = ntry)
    fullfit <- refits[[1]]
    confirmation <- do.call(rbind, lapply(seq_along(kept), function(k) {
      w <- refits[[k + 1L]]
      verdict <- if (inherits(fullfit, 'try-error') || inherits(w, 'try-error'))
        list(dloglik = NA_real_, dpar = NA_integer_, value = NA_real_, accept = NA,
          note = 'refit failed') else
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

# Whether the current fit plainly needs a kept change, so pruning can keep it
# without refitting the model rebuilt without it (`rebuilt`): a Wald test of
# the change's own parameters at their absent values -- a random effect's sd at
# zero, a term's coefficient at zero, a freed cell at the value it was fixed
# at -- on the natural scale, by the delta method from the fit's raw
# covariance, past twice the bar the criterion sets for the parameters the
# change adds. NULL where the change has no such parameters to test: a process
# or cycle added (its persistence has no absent value), loadings on a
# reduced-rank level, covariate effects. On the Fisher search every one of
# eleven pruning refits found its change needed, 12 to 48 log likelihood units
# down, most of them random effects whose sd was many standard errors from
# zero. Near the boundary the Wald test is unreliable, so a change short of
# the margin is refitted.
.ctImprovePruneNeeded <- function(current, prop, rebuilt, cores, criterion, alpha) {
  if (!isTRUE(prop$family %in% c('heterogeneity', 'dimension', 'statedependence',
    'volatility', 'fixed'))) return(NULL)
  cov <- current$estimate$cov
  if (is.null(cov)) return(NULL)
  nm <- as.character(.ctFitRawParNames(current))
  est <- stats::setNames(.ctFitRawEstimate(current), nm)
  reduced <- .ctImproveSpec(current, rebuilt$model, rebuilt$data, cores, rebuilt$args)$names
  added <- setdiff(nm, reduced)
  dpar <- length(nm) - length(reduced)
  if (dpar < 1L || !length(added)) return(NULL)
  model <- .ctFitBaseModel(current)
  # Natural value of each tested parameter as a function of its raw value,
  # and its absent value.
  value <- function(transform, raw) as.numeric(eval(parse(text = transform),
    list(param = raw, log1p_exp = log1p_exp)))
  tested <- list()
  for (a in added) {
    if (startsWith(a, 'popsd_')) {
      tested[[a]] <- list(f = function(r) log1p_exp(2 * r - 1), null = 0)
    } else if (a %in% model$pars$param) {
      row <- model$pars[match(a, model$pars$param), , drop = FALSE]
      tr <- as.character(row$transform)
      if (is.na(tr) || !nzchar(tr)) tr <- 'param'
      cell <- rebuilt$model$pars[rebuilt$model$pars$matrix == row$matrix &
        rebuilt$model$pars$row == row$row & rebuilt$model$pars$col == row$col, , drop = FALSE]
      null <- if (nrow(cell) == 1L && is.na(cell$param)) suppressWarnings(as.numeric(cell$value)) else 0
      if (!is.finite(null)) return(NULL)
      tested[[a]] <- list(f = local({ tr <- tr; function(r) value(tr, r) }), null = null)
    }
  }
  if (!length(tested)) return(NULL)
  idx <- match(names(tested), nm)
  g <- vapply(names(tested), function(a) tested[[a]]$f(est[a]) - tested[[a]]$null, 0)
  J <- vapply(names(tested), function(a) {
    h <- 1e-6 * max(1, abs(est[a]))
    (tested[[a]]$f(est[a] + h) - tested[[a]]$f(est[a] - h)) / (2 * h)
  }, 0)
  V <- cov[idx, idx, drop = FALSE] * outer(J, J)
  W <- tryCatch(as.numeric(t(g) %*% solve(V, g)), error = function(e) NA_real_)
  if (!is.finite(W)) return(NULL)
  bar <- if (criterion == 'LR') stats::qchisq(1 - alpha, dpar) else
    (if (criterion == 'BIC') log(.ctFitNsubjects(current)) else 2) * dpar
  list(needed = W > 2 * bar, W = W, bar = 2 * bar, dpar = dpar)
}

# Whether pruning may try dropping a kept change: a rank change may not
# (`.ctImproveRankChange`).
.ctImprovePrunable <- function(p) !isFALSE(p$prunable)

# The pruning pass: each kept change in turn, the base model rebuilt without
# it and refitted, and the change dropped if the fuller model is not better
# by the criterion. Until a change is dropped every refit is compared against
# the same model, so they are independent: `ntry` at a time run side by side
# when refits go to workers (`.ctImproveSideBySide`), and after a drop the
# rest of that batch, compared against the model it replaced, is refitted
# again. Otherwise one at a time, since there a batch would make refits a
# drop then discards. The decisions are the one-at-a-time pass's either way.
.ctImprovePrune <- function(base, current, kept, startdata, covariates, criterion, alpha,
    cores, ntry, slowest = 0, say = function(...) invisible()) {
  history <- list(); data <- NULL; args <- NULL
  k <- 1L
  while (k <= length(kept) && sum(vapply(kept, .ctImprovePrunable, logical(1))) > 1L) {
    width <- if (.ctImproveSideBySide(current, ntry, cores, slowest)) as.integer(ntry) else 1L
    # The next `width` changes that need a refit to be judged: a change the
    # current fit plainly needs is kept without one (`.ctImprovePruneNeeded`).
    ks <- integer(); rebuilt <- list(); last <- k - 1L
    for (j in seq.int(k, length(kept))) {
      if (length(ks) == width) break
      last <- j
      if (!.ctImprovePrunable(kept[[j]])) next
      rb <- try(.ctImproveRebuild(base, startdata, kept[-j], covariates), silent = TRUE)
      needed <- if (inherits(rb, 'try-error')) NULL else tryCatch(.ctImprovePruneNeeded(
        current, kept[[j]], rb, cores, criterion, alpha), error = function(e) NULL)
      if (isTRUE(needed$needed)) {
        say("Pruning: kept ", kept[[j]]$label, " without a refit (Wald ",
          signif(needed$W, 3), ")")
        history[[length(history) + 1L]] <- data.frame(round = NA_integer_,
          family = 'prune', change = paste0('drop ', kept[[j]]$label), screen_p = NA_real_,
          dloglik = NA_real_, dpar = -needed$dpar, criterion = NA_real_, kept = FALSE,
          note = sprintf('kept without a refit: Wald %.3g against %.3g', needed$W,
            needed$bar), stringsAsFactors = FALSE)
        next
      }
      ks <- c(ks, j); rebuilt[[length(rebuilt) + 1L]] <- rb
    }
    if (!length(ks)) break
    for (j in ks) say("Pruning: refitting without ", kept[[j]]$label)
    refits <- .ctImproveRefits(current, rebuilt, vector('list', length(ks)), cores, slowest)
    slowest <- max(slowest, attr(refits, 'seconds'))
    for (i in seq_along(ks)) {
      j <- ks[i]
      verdict <- if (inherits(refits[[i]], 'try-error')) list(accept = TRUE,
        dloglik = NA_real_, dpar = NA_integer_, value = NA_real_,
        note = 'refit without it failed; kept') else
        .ctImproveCompare(refits[[i]], current, criterion, alpha)
      drop <- !isTRUE(verdict$accept)
      history[[length(history) + 1L]] <- data.frame(round = NA_integer_,
        family = 'prune', change = paste0('drop ', kept[[j]]$label), screen_p = NA_real_,
        dloglik = if (is.finite(verdict$dloglik)) -verdict$dloglik else NA_real_,
        dpar = if (is.finite(verdict$dpar)) -verdict$dpar else NA_integer_,
        criterion = verdict$value, kept = drop, note = verdict$note,
        stringsAsFactors = FALSE)
      k <- if (i == length(ks)) last + 1L else j + 1L
      if (drop) {
        say("Pruning: dropped ", kept[[j]]$label)
        current <- refits[[i]]; data <- rebuilt[[i]]$data; args <- rebuilt[[i]]$args
        kept <- kept[-j]
        k <- j
        break
      }
    }
  }
  list(current = current, kept = kept, data = data, args = args, history = history,
    slowest = slowest)
}

# Whether a fit leaves its population covariance unidentified: a weak or
# negative direction of its curvature (`fit$identifiability$directions`) made
# mostly of the covariance's own parameters. With eleven random effects for 40
# people the Fisher example's search reached a fit whose information along 55
# correlations was about 2.5 against 4e6 elsewhere, with one direction negative:
# the estimation correction then took more than each row-level statistic's
# whole variance, and nothing downstream of that fit could be trusted.
.ctImproveWeakCovariance <- function(fit) {
  dirs <- fit$identifiability$directions
  if (!length(dirs)) return(FALSE)
  any(vapply(dirs, function(d) {
    p <- as.character(d$parameters)
    length(p) > 0L && mean(grepl('^(popsd_|rawcor_|poploading_|popunique_)', p)) >= 0.5
  }, logical(1)))
}

# Rank-r loadings with L L' the best rank-r approximation of the covariance S:
# its leading eigen-directions, rotated so that L[p, q] = 0 for p < q, the
# echelon form the engine reads (`_laplace_poploading`).
.ctImproveEchelonLoadings <- function(S, r) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  B <- e$vectors[, seq_len(r), drop = FALSE] %*%
    diag(sqrt(pmax(e$values[seq_len(r)], 0)), r)
  Q <- qr.Q(qr(t(B[seq_len(r), , drop = FALSE])))
  L <- B %*% Q
  L[upper.tri(L)] <- 0
  L
}

# Raw starting values for a fit of the current model at population rank r:
# the current fit's raw-scale population covariance (`ctsem_laplace_popcov`)
# by echelon loadings, each row divided by its effect's sdscale, named as the
# engine names them. Without it every loading would start at zero, the random
# effects with no spread at all.
.ctImproveRankStart <- function(current, r) {
  spec <- current$model_spec
  module <- .ctJuliaModule(spec$project)
  est <- .ctFitRawEstimate(current)
  S <- as.matrix(.ctBackendJuliaValue(module$ctsem_laplace_popcov(
    .ctJuliaObjective(current), .ctJuliaNumericVector(est), 1L)))
  nm <- as.character(.ctFitRawParNames(current))
  pars <- sub('^popsd_', '', grep('^popsd_', nm, value = TRUE))
  if (length(pars) != nrow(S)) return(NULL)
  model <- .ctFitBaseModel(current)
  scale <- suppressWarnings(as.numeric(model$pars$sdscale[match(pars, model$pars$param)]))
  scale[!is.finite(scale) | scale <= 0] <- 1
  L <- .ctImproveEchelonLoadings(S, r) / scale
  out <- unlist(lapply(seq_len(r), function(q) stats::setNames(L[q:nrow(L), q],
    paste0('poploading_', pars[q:nrow(L)], '_dim', q))))
  out
}

# The number of parameters varying over subjects in a model.
.ctImproveNvarying <- function(model) {
  length(unique(model$pars$param[as.logical(model$pars$indvarying) %in% TRUE &
    !is.na(model$pars$param)]))
}

# A change that sets the population covariance's rank: kept like any other, so
# a rebuild carries it, applied after every other change (`.ctImproveRebuild`)
# so that it sees how many parameters vary by then, and never pruned -- undoing
# it adds parameters rather than removing them. A rank no longer below the
# number varying is the full covariance.
.ctImproveRankChange <- function(r) list(key = 'rank', family = 'rank', p = NA_real_,
  label = paste0('population covariance of rank ', r), prunable = FALSE,
  edit = local({ r <- r; function(model, data, covs, args = list()) {
    k <- .ctImproveNvarying(model)
    list(model = model, data = data, args = list(poprank = if (r < k) r else NA))
  } }))

# The rank a fit's population covariance suggests: the smallest whose dropped
# dimensions carry under `share` of the correlation matrix's trace. On the
# Fisher search's eighth round (nine effects) that was 5, where the
# one-rank-at-a-time step had stopped after four refits, its likelihood ratio
# tests accepting each.
.ctImproveRankTarget <- function(current, share = 0.05) {
  module <- .ctJuliaModule(current$model_spec$project)
  S <- as.matrix(.ctBackendJuliaValue(module$ctsem_laplace_popcov(
    .ctJuliaObjective(current), .ctJuliaNumericVector(.ctFitRawEstimate(current)), 1L)))
  e <- eigen(stats::cov2cor((S + t(S)) / 2), symmetric = TRUE, only.values = TRUE)$values
  which(1 - cumsum(e) / sum(e) < share)[1]
}

# While the current fit leaves its population covariance unidentified, refit
# it at a lower rank, starting from the current covariance's leading
# dimensions, and keep the lower rank unless the likelihood says the dropped
# dimensions are needed: the fuller fit better by the criterion. The first
# refit goes straight to the rank the covariance suggests
# (`.ctImproveRankTarget`); if the likelihood refuses that, the ranks above it
# are tried in turn, up to one below the current. Single-level models only,
# and not below the varying T0MEANS a rank must keep. Returns the current fit
# and kept changes, and a history row for each refit.
.ctImproveRankStep <- function(base, current, kept, startdata, covariates, criterion,
    alpha, cores, say = function(...) invisible(), label = '', round = NA_integer_) {
  history <- list(); data <- NULL; args <- NULL; changed <- FALSE
  refused <- NA_integer_
  repeat {
    if (!.ctImproveWeakCovariance(current)) break
    model <- .ctFitBaseModel(current)
    if (length(model$groupIDnames)) break
    k <- .ctImproveNvarying(model)
    replay <- .ctFitReplayArgs(current)
    r <- suppressWarnings(as.integer(replay$poprank))
    if (!length(r) || is.na(r) || r >= k) r <- k
    t0 <- sum(as.logical(model$pars$indvarying) %in% TRUE & model$pars$matrix %in% 'T0MEANS')
    if (r - 1L < max(1L, t0)) break
    target <- r - 1L
    if (is.na(refused)) {
      suggested <- tryCatch(.ctImproveRankTarget(current), error = function(e) NA_integer_)
      if (isTRUE(suggested < target)) target <- max(suggested, t0, 1L)
    } else target <- refused + 1L
    change <- .ctImproveRankChange(target)
    say(label, "the population covariance is not identified; refitting at rank ", target)
    nonrank <- Filter(function(p) !identical(p$family, 'rank'), kept)
    rebuilt <- try(.ctImproveRebuild(base, startdata, c(nonrank, list(change)), covariates),
      silent = TRUE)
    start <- if (inherits(rebuilt, 'try-error')) NULL else
      tryCatch(.ctImproveRankStart(current, target), error = function(e) NULL)
    refit <- if (inherits(rebuilt, 'try-error')) rebuilt else try(.ctImproveRefit(current,
      rebuilt$model, rebuilt$data, cores, start, rebuilt$args), silent = TRUE)
    verdict <- if (inherits(refit, 'try-error')) list(accept = TRUE, dloglik = NA_real_,
      dpar = NA_integer_, value = NA_real_, note = 'refit at the lower rank failed') else
      .ctImproveCompare(refit, current, criterion, alpha)
    lower <- !isTRUE(verdict$accept)
    history[[length(history) + 1L]] <- data.frame(round = round, family = 'rank',
      change = change$label, screen_p = NA_real_,
      dloglik = if (is.finite(verdict$dloglik)) -verdict$dloglik else NA_real_,
      dpar = if (is.finite(verdict$dpar)) -verdict$dpar else NA_integer_,
      criterion = verdict$value, kept = lower, note = verdict$note, stringsAsFactors = FALSE)
    if (!lower) {
      if (target < r - 1L) { refused <- target; next }
      break
    }
    say(label, "kept ", change$label)
    current <- refit; kept <- c(nonrank, list(change)); changed <- TRUE
    data <- rebuilt$data; args <- rebuilt$args
    refused <- NA_integer_
  }
  list(current = current, kept = kept, data = data, args = args, history = history,
    changed = changed)
}

# The base model with a sequence of kept changes applied, for some data: the
# model, the data the changes need (added covariates), and the ctFit
# arguments they change (popunique, poprank).
.ctImproveRebuild <- function(model, data, kept, covariates) {
  args <- list()
  rank <- vapply(kept, function(p) identical(p$family, 'rank'), logical(1))
  for (p in c(kept[!rank], kept[rank])) {
    e <- p$edit(model, data, covariates, args)
    model <- e$model; data <- e$data
    if (!is.null(e$args)) args <- utils::modifyList(args, e$args)
  }
  list(model = model, data = data, args = args)
}

# What the screen indicates, as edits: one per heterogeneity parameter, state
# term and fixed cell significant within its family, and one per covariate
# whose joint test is significant, adding it with effects on the parameters
# its tests single out. A measurement error covariance is freed only with
# 'errorcov' in `relax`; otherwise the pairs family's evidence proposes a
# process the items share, at the structural stage (`.ctImproveStructural`).
.ctImproveProposals <- function(st, model, alpha, current = NULL,
  relax = c('unique', 'fixed')) {
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
    if (cell$matrix == 'MANIFESTVAR' && !'errorcov' %in% relax) next
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
.ctImproveStructural <- function(model, relax, candidates, ad, current = NULL,
  alpha = 0.05) {
  out <- list()
  # A cycle the periodic family found, at the period its strongest positive
  # statistic names, driving every latent that has innovations of its own.
  if ('periodic' %in% relax && isTRUE(model$continuoustime)) {
    cyc <- .ctImprovePeriodicEvidence(ad, .ctModelCyclePeriods(model))
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
  # Two items whose errors are correlated within an occasion (the pairs
  # family), given a process they share: the reading of a measurement error
  # covariance unless 'errorcov' frees the covariance at the screen. On the
  # own-evidence scale, the pair's p times the pairs tested, and only where
  # that is below alpha, since pairs grow with the square of the items.
  if ('unique' %in% relax && !'errorcov' %in% relax && isTRUE(model$continuoustime)) {
    s <- if (is.null(ad)) NULL else ad$statistics
    pr <- if (is.null(s)) NULL else s[s$family == 'pairs' & is.finite(s$z), , drop = FALSE]
    for (r in seq_len(if (is.null(pr)) 0L else nrow(pr))) {
      items <- strsplit(pr$statistic[r], ':', fixed = TRUE)[[1]]
      if (length(items) != 2L || !all(items %in% model$manifestNames)) next
      items <- model$manifestNames[sort(match(items, model$manifestNames))]
      if (paste0('u_', paste(items, collapse = '_')) %in% model$latentNames) next
      evidence <- min(1, pr$p[r] * nrow(pr))
      if (!isTRUE(evidence < alpha)) next
      out[[length(out) + 1L]] <- list(key = paste0('shared:', paste(items, collapse = ',')),
        family = 'unique', label = paste0(items[1], ' and ', items[2],
          ' given a process they share'), p = evidence, rank = evidence,
        start = .ctImproveUniqueStart(model, current, items, shared = TRUE,
          sign = if (pr$z[r] < 0) -1 else 1),
        edit = local({ items <- items; function(model, data, covs, args = list())
          list(model = ctModelUniqueProcess(model, items, shared = TRUE), data = data) }))
    }
  }
  # A slow process or a second order for a latent with innovations of its
  # own: not for a cycle's clock, which turns as written, nor for the level
  # of a process already second order.
  noisy <- .ctImproveNoisyLatents(model)
  if ('slow' %in% relax && isTRUE(model$continuoustime)) {
    lat <- model$latentNames
    for (i in seq_along(lat)) {
      if (!i %in% noisy || grepl('^slow_', lat[i]) || paste0('slow_', lat[i]) %in% lat ||
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
      if (!i %in% noisy || paste0('d', model$latentNames[i]) %in% model$latentNames) next
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

# Structural proposals scored at the current fit, before any is refitted: a
# score test of each, from the model extended by it at the point where it has
# no effect. A cycle's clock turns as written, so its couplings into the
# processes have an ordinary score at zero: their joint robust test, with the
# estimation removed exactly from the few couplings by the Hessian's columns
# along them. A process added at no effect -- a slow process uncoupled, a
# process of an item's own with no variance -- has a zero first derivative
# there, so its test is of the second, along the proposal's own start (its
# effect parameters scaled by kappa, the rest at their starts): each unit's
# log likelihood at kappa = 0.1 against the fit's gives the curvature
# q_i = 2 (l_i(0.1) - l_i(0)) / 0.01, and the one-sided test is of their sum,
# the score for kappa^2, on their spread over units (Chesher, 1984). The start
# itself gives the gain the refit begins with, a lower bound on its likelihood
# ratio. Proposals are ranked by these p-values; one whose statistic points
# the wrong way (p >= 0.5) is not refitted. Second order and candidate models
# have no such point and keep their place. Each test costs a model prepared
# and a few evaluations, where a refit of the Fisher example's late models
# took the better part of an hour.
.ctImproveScoreStructural <- function(current, proposals, data, covariates, cores) {
  skipped <- list()
  if (!.ctFitIsJulia(current)) return(list(proposals = proposals, skipped = skipped))
  est <- .ctFitRawEstimate(current)
  nm1 <- as.character(.ctFitRawParNames(current))
  base <- tryCatch(as.numeric(ctJuliaEvaluate(current, est, gradient = FALSE,
    contributions = TRUE, cores = cores)$subject_loglik), error = function(e) NULL)
  keep <- rep(TRUE, length(proposals))
  for (k in seq_along(proposals)) {
    prop <- proposals[[k]]
    if (!isTRUE(prop$family %in% c('periodic', 'slow', 'unique'))) next
    sc <- tryCatch(.ctImproveStructuralTest(current, prop, data, covariates, cores, est,
      nm1, base), error = function(e) NULL)
    if (is.null(sc) || !is.finite(sc$p)) next
    proposals[[k]]$p <- sc$p
    proposals[[k]]$rank <- sc$p
    proposals[[k]]$scored <- sc
    if (sc$p >= 0.5) {
      keep[k] <- FALSE
      skipped[[length(skipped) + 1L]] <- list(key = prop$key, family = prop$family,
        label = prop$label, p = sc$p, note = sprintf(
          'not refitted: score test p = %.3g%s', sc$p, if (is.finite(sc$gain))
            sprintf(', gain at its start %.3g', sc$gain) else ''))
    }
  }
  list(proposals = proposals[keep], skipped = skipped)
}

# One structural proposal's score test (`.ctImproveScoreStructural`): its p,
# the statistic, and for an added process the gain at its start. NULL where
# the extended model cannot be prepared or does not reproduce the fit.
.ctImproveStructuralTest <- function(current, prop, data, covariates, cores, est, nm1, base) {
  model <- .ctFitBaseModel(current)
  e <- prop$edit(model, data, covariates)
  sp <- .ctImproveSpec(current, e$model, e$data, cores, e$args)
  nm2 <- sp$names
  if (!all(nm1 %in% nm2)) return(NULL)
  x0 <- stats::setNames(numeric(length(nm2)), nm2)
  x0[nm1] <- est
  pars <- e$model$pars
  module <- .ctJuliaModule(current$model_spec$project)
  newlat <- which(!e$model$latentNames %in% model$latentNames)
  if (prop$family == 'periodic') {
    # The couplings into the processes, at zero, where the extended model must
    # reproduce the fit -- up to the normal(0, 1) prior's constant for each new
    # parameter, with priors on every parameter.
    new <- setdiff(nm2, nm1)
    at <- pars[match(new, pars$param), , drop = FALSE]
    test <- new[at$matrix %in% 'DRIFT' & !at$row %in% newlat & at$col %in% newlat]
    if (!length(test)) return(NULL)
    value1 <- ctJuliaEvaluate(current, est, gradient = FALSE, cores = cores)$value
    offset <- if (isTRUE(current$args$resolved$priors)) length(new) * stats::dnorm(0, log = TRUE) else 0
    ev <- ctJuliaEvaluate(sp$spec, x0, gradient = TRUE, cores = cores)
    if (!isTRUE(abs(ev$value - value1 - offset) <= 1e-6 * max(1, abs(value1)))) return(NULL)
    objective <- .ctJuliaObjective(sp$spec, hybrid = TRUE)
    units <- as.matrix(.ctBackendWithMaxChunks(cores, .ctJuliaGet(module$ctsem_subject_gradients(
      objective, .ctJuliaVector(as.numeric(x0)))))$scores)
    theta <- match(nm1, nm2); delta <- match(test, nm2)
    Hc <- .ctBackendWithMaxChunks(cores, as.matrix(.ctBackendJuliaValue(
      module$ctsem_hessian_columns(objective, .ctJuliaNumericVector(as.numeric(x0)),
        .ctJuliaVector(as.integer(delta))))))
    Iinv <- .ctScorePinv(-.ctBackendHessian(current, est))$inv
    A <- -t(Hc[theta, , drop = FALSE]) %*% Iinv
    grad <- as.numeric(ev$gradient)
    U <- grad[delta] - as.numeric(A %*% grad[theta])
    C <- units[, delta, drop = FALSE] - units[, theta, drop = FALSE] %*% t(A)
    j <- .ctScorePinv(.ctScoreSpread(C))
    if (!j$rank) return(NULL)
    stat <- as.numeric(t(U) %*% j$inv %*% U)
    return(list(p = .ctScoreBoot(C, 999L), statistic = stat, df = j$rank, gain = NA_real_))
  }
  # An added process: its effect parameters scaled by kappa along the start.
  start <- prop$start
  if (is.null(base) || !length(start) || !all(names(start) %in% nm2)) return(NULL)
  value <- function(transform, raw) as.numeric(eval(parse(text = transform),
    list(param = raw, log1p_exp = log1p_exp)))
  inverse <- function(transform, target) stats::uniroot(function(r)
    value(transform, r) - target, c(-30, 30), tol = 1e-12)$root
  cells <- pars[match(names(start), pars$param), , drop = FALSE]
  effect <- cells$matrix %in% c('DIFFUSION', 'T0VAR') |
    (cells$matrix %in% 'DRIFT' & !cells$row %in% newlat & cells$col %in% newlat)
  if (!any(effect)) return(NULL)
  at <- function(kappa) {
    x <- x0
    x[names(start)[!effect]] <- start[!effect]
    for (q in which(effect)) {
      tr <- as.character(cells$transform[q])
      if (is.na(tr) || !nzchar(tr)) tr <- 'param'
      x[names(start)[q]] <- inverse(tr, kappa * value(tr, start[[q]]))
    }
    x
  }
  units <- function(x) as.numeric(ctJuliaEvaluate(sp$spec, x, gradient = FALSE,
    contributions = TRUE, cores = cores)$subject_loglik)
  kappa <- 0.1
  q <- 2 * (units(at(kappa)) - base) / kappa^2
  if (length(q) != length(base) || any(!is.finite(q))) return(NULL)
  c2 <- q / 2
  v <- stats::var(c2) * length(c2)
  if (!is.finite(v) || v <= 0) return(NULL)
  z <- sum(c2) / sqrt(v)
  gain <- tryCatch(sum(units(at(1)) - base), error = function(e) NA_real_)
  list(p = stats::pnorm(z, lower.tail = FALSE), statistic = z, df = 1L, gain = gain)
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
# the Fisher data. Periods within one per cent of one in `exclude` -- a cycle
# the model has -- are passed over: a second 24 hour clock was proposed after
# the first was kept. NULL when the family did not run or found no positive
# statistic at another period.
.ctImprovePeriodicEvidence <- function(ad, exclude = numeric()) {
  if (is.null(ad) || is.null(ad$families)) return(NULL)
  f <- ad$families[ad$families$family == 'periodic', , drop = FALSE]
  s <- ad$statistics[ad$statistics$family == 'periodic' & is.finite(ad$statistics$z) &
    ad$statistics$z > 0, , drop = FALSE]
  at <- as.numeric(sub('.* period ', '', s$statistic))
  s <- s[is.finite(at) & !vapply(at, function(a) any(abs(a - exclude) <= 0.01 * a),
    logical(1)), , drop = FALSE]
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
# transform. NULL, and raw zero, when they cannot be worked out. A process
# `shared` by `item` (two or more) is placed on the first; each other loading
# starts where the item gets the same share of its own error variance, with
# the `sign` of the covariance the screen pointed to.
.ctImproveUniqueStart <- function(model, current, item, shared = FALSE, sign = 1) {
  if (is.null(current)) return(NULL)
  j <- match(item, model$manifestNames)
  m1 <- try(ctModelUniqueProcess(model, item, shared = shared), silent = TRUE)
  if (anyNA(j) || inherits(m1, 'try-error')) return(NULL)
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
  est <- setNames(.ctFitRawEstimate(current), .ctFitRawParNames(current))
  errsd <- function(jj) {
    cell <- model$pars[model$pars$matrix == 'MANIFESTVAR' & model$pars$row == jj &
      model$pars$col == jj, , drop = FALSE]
    if (nrow(cell) != 1L) return(NA_real_)
    if (is.na(cell$param)) as.numeric(cell$value) else value(cell$transform, est[[cell$param]])
  }
  sd <- vapply(j, errsd, 0)
  if (!is.finite(gap) || !all(is.finite(sd)) || any(sd <= 0)) return(NULL)
  rate <- log(2) / gap
  v <- sd[1]^2 / 4
  u <- length(m1$latentNames)
  target <- data.frame(matrix = c('DRIFT', 'DIFFUSION', 'T0VAR'), row = u,
    value = c(-rate, sqrt(2 * rate * v), sqrt(v)))
  if (isTRUE(shared)) target <- rbind(target, data.frame(matrix = 'LAMBDA', row = j[-1],
    value = sign * sd[-1] / sd[1]))
  out <- numeric()
  for (k in seq_len(nrow(target))) {
    r <- m1$pars[m1$pars$matrix == target$matrix[k] & m1$pars$row == target$row[k] &
      m1$pars$col == u, , drop = FALSE]
    if (nrow(r) != 1L || is.na(r$param)) next
    raw <- inverse(r$transform, target$value[k])
    if (is.finite(raw)) out[r$param] <- raw
  }
  if (length(out)) out else NULL
}

# Each edited model refitted (`.ctImproveRefit`), a try-error where the edit or
# the refit failed, with the longest refit made in this session as attribute
# `seconds`. Side by side in worker processes when there are two or more, the
# cores to share, and a refit has already taken `slowest` seconds above
# `.ctWorkerFitSeconds`: the refits of a round are independent, and one fit of
# the Fisher example's model with random dynamics kept about 6 of its 10 cores
# busy. The threshold is the one every worker decision shares: a worker takes
# 30-60 s to start, and three refits of a few seconds each took 46 s in fresh
# workers against 27 s here. Each worker takes `cores %/% workers` threads. The
# pool is the sampling chains' (`.ctBackendWarmWorkers`), which persists, so
# later rounds find it started; a fit made there rebuilds its objective in this
# session when it is next evaluated, the cache being keyed on the specification.
# At most `width` workers: past that a refit waits for a free one (`future()`
# blocks), and each worker's Julia holds its own objective, several GB on the
# Fisher example.
.ctImproveRefits <- function(current, edits, starts, cores, slowest = 0, width = Inf) {
  seconds <- 0
  here <- function(i) {
    if (inherits(edits[[i]], 'try-error')) return(edits[[i]])
    started <- proc.time()[['elapsed']]
    out <- try(.ctImproveRefit(current, edits[[i]]$model, edits[[i]]$data, cores,
      starts[[i]], edits[[i]]$args), silent = TRUE)
    seconds <<- max(seconds, proc.time()[['elapsed']] - started)
    out
  }
  inhere <- function() structure(lapply(seq_along(edits), here), seconds = seconds)
  ok <- which(!vapply(edits, inherits, logical(1), 'try-error'))
  if (!.ctImproveSideBySide(current, min(length(ok), width), cores, slowest)) return(inhere())
  workers <- min(length(ok), as.integer(cores), width)
  threads <- max(1L, as.integer(cores) %/% workers)
  handles <- .ctBackendWarmWorkers(current, workers, values = .ctFitRawEstimate(current),
    threads = threads)
  if (.ctBackendWarmWait(handles) < 1L) return(inhere())
  futures <- lapply(seq_along(edits), function(i) if (!i %in% ok) NULL else
    tryCatch(future::future(utils::getFromNamespace('.ctImproveRefit', 'ctsem')(current,
      edits[[i]]$model, edits[[i]]$data, threads, starts[[i]], edits[[i]]$args),
      seed = TRUE), error = function(e) NULL))
  out <- lapply(seq_along(edits), function(i) {
    if (!i %in% ok) return(edits[[i]])
    if (is.null(futures[[i]])) return(here(i))
    tryCatch(future::value(futures[[i]]), error = function(e)
      structure(conditionMessage(e), class = 'try-error', condition = e))
  })
  structure(out, seconds = seconds)
}

# Whether `.ctImproveRefits` would send `n` refits to worker processes rather
# than make them here: two or more of them and of the cores, a refit here has
# taken longer than `.ctWorkerFitSeconds`, and the fit is on julia with a pool
# to hand.
.ctImproveSideBySide <- function(current, n, cores, slowest) {
  min(n, as.integer(cores)) >= 2L && isTRUE(slowest > .ctWorkerFitSeconds) &&
    isTRUE(getOption('ctsem.improve.processes', TRUE)) &&
    .ctFitIsJulia(current) && .ctBackendCanWarm()
}

# Refit with an edited model, prepared as the current fit was but on `cores`,
# and warm started from its estimates, matched by name; new parameters start at
# zero, or where the proposal says (`start`, by raw name). A worker's refit gets
# its share of the cores; replaying the fit's own `cores` gave each worker all
# of them. A warm start can still fail -- a
# new random effect at an arbitrary scale can leave a Laplace unit's inner mode
# unconverged at the start, which ended two of three refits in testing -- so a
# refit that fails from there is tried once from ctsem's own start.
.ctImproveRefit <- function(current, model, data, cores, start = NULL, args = NULL) {
  prepared <- .ctImproveSpec(current, model, data, cores, args)
  args <- prepared$args
  names2 <- prepared$names
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

# An edited model prepared as the current fit was, on `cores`, and not fitted:
# the ctFit arguments and the raw parameter names.
.ctImproveSpec <- function(current, model, data, cores, args = NULL) {
  extra <- args
  args <- .ctFitReplayArgs(current)
  if (length(extra)) args <- utils::modifyList(args, extra)
  args$model <- model
  args$datalong <- data
  args$cores <- cores
  args$verbose <- 0
  args$fit <- FALSE
  args$inits <- NULL
  spec <- suppressMessages(suppressWarnings(do.call(ctFit, args)))
  list(args = args, spec = spec, names = as.character(.ctBackendRawParameterNames(
    list(model_spec = spec), .ctBackendNpar(spec))))
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
