# A model parameter a fit left run off where the likelihood is flat in it.
#
# A local optimiser settles in the basin its start leads to. Where that basin
# takes a parameter out to a flat end of its range -- an auto-effect so fast
# that its process is white noise at the observation interval, or so slow that
# it is a random walk -- the gradient no longer brings it back, its value says
# nothing, and the basin on the other side was never compared. On
# follhim/dis_ctsem (2026-11-21) a trend's auto-effect ran off to -64 per hour,
# the trend collapsing into each person's level, 24 log likelihood units below
# the linear trend the same model reaches with that auto-effect slow; between
# the two the profile dips by 1.6, which no local step crosses.
#
# So where a fit ends with such a parameter, it is refitted with that
# parameter held at the other side of its range (raw `-sign(raw) *
# .ctDataStartRawBound`, as far out as the data start goes) until the rest
# settle, without a Hessian; then the model is refitted from there with the
# parameter free, and the better of the two fits is kept. A parameter qualifies
# when the fit's report has it undetermined (most of its variation in the flat
# directions, as the reduction reads them: `.ctReduceFindings()`) or its
# transform saturated, and it has gone past the data start's raw bound. Each
# parameter is tried once a fit, so one that runs off to the other end then is
# not sent back. `optimcontrol$flatretry = FALSE` skips it; the refits run with
# it off.

# The model parameters `fit` left run off where the likelihood is flat.
.ctFitFlatRetryParameters <- function(fit) {
  model <- .ctFitBaseModel(fit)
  if (is.null(model)) return(character())
  found <- tryCatch(.ctReduceFindings(fit, candidates = FALSE), error = function(e) list())
  flat <- unlist(lapply(Filter(function(f) identical(f$kind, 'parameters'), found),
    `[[`, 'parameters'))
  saturated <- as.character(fit$optim$saturated_parameters)
  params <- unique(c(flat, saturated))
  params <- params[!is.na(.ctReduceCells(model, params)$matrix)]
  if (!length(params)) return(character())
  raw <- tryCatch(.ctFitNameRaw(fit, fit$estimate$raw), error = function(e) NULL)
  if (is.null(raw)) return(character())
  params <- params[params %in% names(raw)]
  params[abs(raw[params]) > .ctDataStartRawBound]
}

# The value `param` takes, on its own scale, at raw `rawvalue` with every
# other coordinate at the estimate (its cells evaluated as `.ctReduceValues()`
# evaluates them).
.ctFitFlatRetryValue <- function(fit, param, rawvalue) {
  raw <- .ctFitNameRaw(fit, fit$estimate$raw)
  k <- match(param, names(raw))
  m <- matrix(raw, nrow = 1L)
  m[1L, k] <- rawvalue
  mapped <- .ctBackendPopMeanSamples(fit, samples = m)
  column <- match(k, as.integer(mapped$parnumber))
  as.numeric(as.matrix(mapped$values)[1L, column])
}

# One refit of `model` from `current`'s estimate (by name, `start` over it),
# with `optimcontrol` added to the fit's own and the retry off. Its warnings
# are kept, not raised, for whichever fit is kept to raise.
.ctFitFlatRetryRefit <- function(current, model, data, cores, start = NULL,
    optimcontrol = list()) {
  oc <- c(list(flatretry = FALSE, checkpoint = NULL), optimcontrol)
  prepared <- .ctImproveSpec(current, model, data, cores, list(optimcontrol = oc))
  args <- prepared$args
  names2 <- prepared$names
  est <- stats::setNames(.ctFitRawEstimate(current), .ctFitRawParNames(current))
  init <- stats::setNames(numeric(length(names2)), names2)
  shared <- intersect(names(est), names2)
  init[shared] <- est[shared]
  start <- start[names(start) %in% names2]
  init[names(start)] <- start
  args$fit <- TRUE
  args$inits <- unname(init)
  warned <- list()
  fit <- withCallingHandlers(suppressMessages(do.call(ctFit, args)),
    warning = function(w) {
      warned[[length(warned) + 1L]] <<- w
      invokeRestart('muffleWarning')
    })
  list(fit = fit, warnings = warned)
}

# Hold-then-free for each parameter `fit` left run off where the likelihood is
# flat, keeping the better fit each time. `warned` are the warnings of `fit`,
# held back; the result carries the kept fit's. A record of what was tried is
# on the kept fit as `$flatretry`.
.ctFitFlatRetry <- function(fit, warned, cores, verbose = 0) {
  data <- tryCatch(as.data.frame(.ctFitLongData(fit)), error = function(e) NULL)
  base <- .ctFitBaseModel(fit)
  if (is.null(data) || is.null(base)) return(list(fit = fit, warnings = warned))
  call <- fit$args
  estonly <- isTRUE(.ctFitReplayArgs(fit)$optimcontrol$estonly)
  record <- list()
  tried <- character()
  repeat {
    params <- setdiff(.ctFitFlatRetryParameters(fit), tried)
    if (!length(params)) break
    p <- params[1L]
    tried <- c(tried, p)
    raw <- .ctFitNameRaw(fit, fit$estimate$raw)
    rawheld <- -sign(raw[[p]]) * .ctDataStartRawBound
    value <- tryCatch(.ctFitFlatRetryValue(fit, p, rawheld), error = function(e) NA_real_)
    if (!is.finite(value)) next
    from <- tryCatch(.ctFitFlatRetryValue(fit, p, raw[[p]]), error = function(e) NA_real_)
    before <- .ctFitOptimisedLoglik(fit)
    message(sprintf(paste0("%s ran off to %s, where the likelihood is flat in it: ",
      "refitting with it held at %s, then freed."), p, format(from, digits = 3),
      format(value, digits = 3)))
    edit <- if (p %in% .ctImproveVarying(base)) .ctReduceStateEdit(p, value) else
      .ctReduceFixEdit(p, value)
    held <- tryCatch(.ctFitFlatRetryRefit(fit, edit(base, data, NULL)$model, data,
      cores, optimcontrol = list(estonly = TRUE)), error = function(e) NULL)
    freed <- if (is.null(held)) NULL else tryCatch({
      # Coordinates the hold removed (the parameter, and its random and
      # predictor effects where it was fixed) start where the fit had them,
      # the parameter at the held value.
      start <- stats::setNames(.ctFitRawEstimate(fit), .ctFitRawParNames(fit))
      start <- start[!names(start) %in% .ctFitRawParNames(held$fit)]
      start[p] <- rawheld
      # As the fit was asked for: the held fit's call is `estonly`.
      .ctFitFlatRetryRefit(held$fit, base, data, cores, start = start,
        optimcontrol = list(estonly = estonly))
    }, error = function(e) NULL)
    after <- if (is.null(freed)) NA_real_ else .ctFitOptimisedLoglik(freed$fit)
    kept <- isTRUE(after > before)
    record[[length(record) + 1L]] <- data.frame(parameter = p, from = from,
      held = value, loglik_before = before,
      loglik_held = if (is.null(held)) NA_real_ else .ctFitOptimisedLoglik(held$fit),
      loglik_after = after, kept = kept, stringsAsFactors = FALSE)
    if (kept) {
      message(sprintf("Kept the refit: log likelihood %s -> %s.",
        format(before, nsmall = 2), format(after, nsmall = 2)))
      fit <- freed$fit
      warned <- freed$warnings
    } else {
      message(sprintf("The fit stands: the refit reached %s against %s.",
        if (is.finite(after)) format(after, nsmall = 2) else "no fit",
        format(before, nsmall = 2)))
    }
  }
  if (length(record)) {
    # The kept fit answers the caller's call, whichever refit reached it.
    fit$args <- call
    fit$flatretry <- do.call(rbind, record)
  }
  list(fit = fit, warnings = warned)
}
