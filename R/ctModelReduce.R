# Reducing a model whose fit leaves directions the data do not identify.
#
# A fit's identifiability report (`fit$identifiability$directions`) lists the
# directions along which the likelihood is flat, each with the parameters it
# is made of. A few kinds recur, and each has a model change that removes it:
#
#   unread        a parameter the likelihood does not read (a measurement
#                 error covariance between items scored one at a time, before
#                 that was refused): curvature exactly zero along it alone.
#                 Fixed at zero.
#   drift-intercept / drift-diffusion
#                 a latent's auto-effect against its intercept or its
#                 diffusion: the process settles faster than the data
#                 resolve, so only their ratio is determined. A latent
#                 without noise (a trend) is refitted as a linear trend, its
#                 auto-effect fixed at zero; any is refitted with its
#                 intercept or diffusion scaled by the drift (`driftScaled`),
#                 which states the same model along the ratio.
#   covariance    population sds, correlations or loadings: the covariance
#                 has more dimensions than the subjects identify. Refitted at
#                 lower ranks (`.ctImproveRankStep()`), from full rank with
#                 own variances (`popunique`).
#   no spread / weak
#                 random effects the data hold at zero, or barely determine
#                 (`fit$identifiability$effects`). Made fixed.
#
# `ctModelImprove(fit, reduceOnly = TRUE)` refits each change from the current
# estimate and keeps it unless the fuller model is better by the search's
# criterion; the full search does the same before it screens and after each
# change it keeps. On follhim/dis_ctsem (2026-11-21) the original fit's report
# held all four kinds at once: an ignored error covariance, a stress trend
# whose auto-effect had run to -37 an hour along its intercept (fixed at zero
# it fitted 21.6 nats better with three fewer parameters), and a 10 x 10
# random-effects covariance of about rank 6.

# A raw coordinate of a population covariance.
.ctReducePopCoordinate <- function(p) grepl('^(popsd_|rawcor_|poploading_|popunique_)', p)

# The cells of model parameters, by name: matrix, row and column of the cell
# that holds each (`.ctParamRows()`), NA where the name is not a parameter of
# the model's table (a predictor effect, a population coordinate).
.ctReduceCells <- function(model, params) {
  rows <- .ctParamRows(model, params)
  pars <- model$pars
  data.frame(param = params, matrix = pars$matrix[rows], row = pars$row[rows],
    col = pars$col[rows], stringsAsFactors = FALSE)
}

# The parameter name of a model cell, or NA where the cell is fixed or an
# expression.
.ctReduceCellParam <- function(model, matrix, row, col) {
  pars <- model$pars
  r <- which(pars$matrix == matrix & pars$row == row & pars$col == col)
  if (length(r) != 1L || is.na(pars$param[r]) || !is.na(pars$value[r]) ||
      !grepl('^[A-Za-z][A-Za-z0-9_.]*$', pars$param[r])) return(NA_character_)
  pars$param[r]
}

# Whether latent i has no noise of its own: every DIFFUSION cell in its row and
# column fixed at zero.
.ctReduceNoNoise <- function(model, i) {
  pars <- model$pars
  r <- which(pars$matrix == 'DIFFUSION' & (pars$row == i | pars$col == i))
  length(r) > 0L && all(is.na(pars$param[r]) &
    suppressWarnings(as.numeric(pars$value[r])) %in% 0)
}

# A reduction as `.ctImproveRebuild()` applies it: an edit of the model, never
# pruned (undoing it adds back what the data did not identify).
.ctReduceProposal <- function(key, label, edit) list(key = key, family = 'reduce',
  label = label, p = NA_real_, rank = 0, prunable = FALSE, edit = edit)

# `params` fixed at `value` in every cell that holds them, with no random or
# predictor effect.
.ctReduceFixEdit <- function(params, value) {
  force(params); force(value)
  function(model, data, covs, args = list()) {
    pars <- model$pars
    r <- which(pars$param %in% params)
    pars$value[r] <- value
    pars$param[r] <- NA_character_
    pars$transform[r] <- NA_character_
    model$pars <- pars
    model <- .ctReduceNoVariation(model, r)
    list(model = model, data = data, args = list())
  }
}

# Rows `r` of a model's table with no random or predictor effect.
.ctReduceNoVariation <- function(model, r) {
  pars <- model$pars
  pars$indvarying[r] <- if (is.character(pars$indvarying)) 'FALSE' else FALSE
  for (cc in intersect(paste0(model$TIpredNames, '_effect'), names(pars)))
    pars[[cc]][r] <- .ctDriftScaledFalse(pars[[cc]][r])
  model$pars <- pars
  .ctModelRawPopVarSync(model)
}

# Random effects removed: `params` vary over no subject.
.ctReduceFixedEffectsEdit <- function(params) {
  force(params)
  function(model, data, covs, args = list()) {
    list(model = .ctReduceNoVariation(model, which(model$pars$param %in% params)),
      data = data, args = list())
  }
}

# A matrix's cells written as what they accrue over one time unit
# (`.ctModelDriftScaled()`), beside any already so.
.ctReduceDriftScaledEdit <- function(which) {
  force(which)
  function(model, data, covs, args = list()) {
    model <- suppressMessages(.ctModelDriftScaled(model,
      unique(c(model$driftScaled, which))))
    list(model = model, data = data, args = list())
  }
}

# What a fit's identifiability report says, sorted into the kinds above: one
# finding per kind and latent, each with the parameters involved, a phrase
# saying what it is, and the changes that would remove it (`candidates`; none
# for a kind the rank step handles or for a direction of no known kind).
.ctReduceFindings <- function(fit) {
  model <- .ctFitBaseModel(fit)
  id <- fit$identifiability
  if (is.null(model) || is.null(id)) return(list())
  zero <- .ctEffectRows(id$effects, 'zero')
  weak <- .ctEffectRows(id$effects, 'weak')
  held <- c(if (!is.null(zero)) zero$coordinate, if (!is.null(weak)) weak$coordinate)
  out <- list()
  add <- function(f) if (!f$key %in% vapply(out, `[[`, '', 'key')) out[[length(out) + 1L]] <<- f
  latent <- function(i) model$latentNames[i]
  for (d in id$directions) {
    # A ridge of the augmented route between a population sd and its
    # correlations (`$partial`) is that route's, which the Laplace route
    # identifies (`.ctIdentifyAdvice()`): no change to the model.
    if (!is.null(d$partial)) next
    p <- setdiff(as.character(d$parameters), held)
    if (!length(p)) next
    pop <- .ctReducePopCoordinate(p)
    if (mean(pop) >= 0.5) {
      add(list(key = 'covariance', kind = 'covariance', parameters = p,
        text = 'the population covariance has more dimensions than the subjects identify',
        candidates = list()))
      next
    }
    mp <- p[!pop]
    cells <- .ctReduceCells(model, mp)
    found <- FALSE
    if (length(p) == 1L && isTRUE(abs(d$relative) < 1e-12) && !is.na(cells$matrix[1L])) {
      # Zero curvature, to roundoff, along one parameter alone: nothing the
      # likelihood computes reads it (`.ctOptimIdentifiedInverse()`'s 1e-12 for
      # a direction with no width at all). An off-diagonal or an effect is
      # fixed at zero; any other cell is named, not changed.
      offdiag <- cells$row[1L] != cells$col[1L] ||
        cells$matrix[1L] %in% c('TDPREDEFFECT', 'LAMBDA')
      add(list(key = paste0('unread:', mp), kind = 'unread', parameters = mp,
        text = paste0(mp, ' is not read by the likelihood'),
        candidates = if (offdiag) list(.ctReduceProposal(paste0('reduce:fix:', mp),
          paste0(mp, ' fixed at 0'), .ctReduceFixEdit(mp, 0))) else list()))
      next
    }
    for (i in unique(cells$row[cells$matrix %in% 'DRIFT' & cells$row == cells$col])) {
      drift <- cells$param[cells$matrix %in% 'DRIFT' & cells$row == i & cells$col == i][1L]
      cint <- .ctReduceCellParam(model, 'CINT', i, 1L)
      diffusion <- .ctReduceCellParam(model, 'DIFFUSION', i, i)
      if (!is.na(cint) && cint %in% mp) {
        found <- TRUE
        cand <- list()
        if (.ctReduceNoNoise(model, i)) cand[[length(cand) + 1L]] <- .ctReduceProposal(
          paste0('reduce:lineartrend:', drift), paste0(drift, ' fixed at 0 (a linear trend)'),
          .ctReduceFixEdit(drift, 0))
        if (!'CINT' %in% model$driftScaled) cand[[length(cand) + 1L]] <- .ctReduceProposal(
          'reduce:driftscaled:CINT', "intercepts scaled by the drift (driftScaled = 'CINT')",
          .ctReduceDriftScaledEdit('CINT'))
        add(list(key = paste0('driftcint:', i), kind = 'drift-intercept',
          parameters = c(drift, cint), text = paste0(drift, ' and ', cint,
            ' trade off: ', latent(i), ' settles faster than the data resolve'),
          candidates = cand))
      }
      if (!is.na(diffusion) && diffusion %in% mp) {
        found <- TRUE
        add(list(key = paste0('driftdiffusion:', i), kind = 'drift-diffusion',
          parameters = c(drift, diffusion), text = paste0(drift, ' and ', diffusion,
            ' trade off: ', latent(i), ' forgets its state faster than the data resolve'),
          candidates = if ('DIFFUSION' %in% model$driftScaled) list() else
            list(.ctReduceProposal('reduce:driftscaled:DIFFUSION',
              "diffusions scaled by the drift (driftScaled = 'DIFFUSION')",
              .ctReduceDriftScaledEdit('DIFFUSION')))))
      }
    }
    if (!found) add(list(key = paste0('other:', paste(sort(p), collapse = ',')),
      kind = 'other', parameters = p, text = paste(utils::head(p, 4), collapse = ', '),
      candidates = list()))
  }
  # Random effects the data hold at zero, or barely determine: on a model with
  # one level only, where `indvarying = FALSE` means what it says.
  single <- !length(model$groupIDnames)
  for (w in list(list(rows = zero, kind = 'no spread',
      text = ': the data hold their individual differences at zero'),
    list(rows = weak, kind = 'weak', text = ': the data barely determine their individual differences'))) {
    if (is.null(w$rows) || !nrow(w$rows) || !single) next
    eff <- unique(as.character(w$rows$effect))
    add(list(key = paste0(w$kind, ':', paste(sort(eff), collapse = ',')), kind = w$kind,
      parameters = eff, text = paste0(paste(eff, collapse = ', '), w$text),
      candidates = list(.ctReduceProposal(paste0('reduce:fixed:', paste(sort(eff), collapse = ',')),
        paste0(paste(eff, collapse = ', '), ' not varying over subjects'),
        .ctReduceFixedEffectsEdit(eff)))))
  }
  out
}

# Whether to keep a reduction: unless the model it reduces is better by the
# criterion. A reduction with fewer parameters is tested as the rank step
# tests a lower rank; one with as many (a reparameterisation) is kept unless it
# ends below the fit it restates by more than an optimiser's tolerance. A
# fuller fit below the reduced one was not at its optimum -- it had stopped on
# the ridge the reduction removes -- and the reduction is kept.
.ctReduceVerdict <- function(current, reduced, criterion, alpha, cores) {
  fail <- function(note) list(keep = FALSE, dloglik = NA_real_, dpar = NA_integer_,
    value = NA_real_, note = note)
  if (inherits(reduced, 'try-error')) return(fail(paste('failed:',
    trimws(conditionMessage(attr(reduced, 'condition'))))))
  ll0 <- .ctFitOptimisedLoglik(current); ll1 <- .ctFitOptimisedLoglik(reduced)
  k0 <- length(.ctFitRawEstimate(current)); k1 <- length(.ctFitRawEstimate(reduced))
  if (!is.finite(ll1)) return(fail('log likelihood not finite'))
  tol <- 0.01 * max(1, sqrt(abs(ll0)) / 10)
  out <- list(dloglik = ll1 - ll0, dpar = k1 - k0)
  if (k1 >= k0) return(c(out, list(keep = ll1 >= ll0 - tol, value = NA_real_,
    note = if (ll1 >= ll0 - tol) 'restates the model' else
      sprintf('restated, %.3g below', ll0 - ll1))))
  if (ll1 > ll0 + tol) return(c(out, list(keep = TRUE, value = NA_real_,
    note = sprintf('%.3g better with fewer parameters', ll1 - ll0))))
  v <- .ctImproveCompare(reduced, current, criterion, alpha, cores)
  if (!is.finite(v$dloglik)) return(fail(v$note))
  c(out, list(keep = !isTRUE(v$accept), value = v$value, note = v$note))
}

# The reduction stage: while the current fit's report has a finding with a
# change not yet refused, refit each of its changes, keep the best the
# criterion allows, and read the new fit's report; then lower the population
# covariance's rank while it is not identified (`.ctImproveRankStep()`).
# Returns as `.ctImproveRankStep()` does.
.ctImproveReduce <- function(base, current, kept, startdata, covariates, criterion,
    alpha, cores, say = function(...) invisible(), label = '', round = NA_integer_) {
  history <- list(); data <- NULL; args <- NULL; changed <- FALSE
  refused <- character()
  repeat {
    findings <- Filter(function(f) length(f$candidates) && !f$key %in% refused,
      .ctReduceFindings(current))
    if (!length(findings)) break
    f <- findings[[1L]]
    say(label, f$text)
    best <- NULL
    for (cand in f$candidates) {
      say(label, "refitting with ", cand$label)
      rebuilt <- try(.ctImproveRebuild(base, startdata, c(kept, list(cand)), covariates),
        silent = TRUE)
      refit <- if (inherits(rebuilt, 'try-error')) rebuilt else try(.ctImproveRefit(current,
        rebuilt$model, rebuilt$data, cores, NULL, rebuilt$args), silent = TRUE)
      verdict <- .ctReduceVerdict(current, refit, criterion, alpha, cores)
      history[[length(history) + 1L]] <- data.frame(round = round, family = 'reduce',
        change = cand$label, screen_p = NA_real_, dloglik = verdict$dloglik,
        dpar = verdict$dpar, criterion = verdict$value, kept = FALSE,
        note = verdict$note, stringsAsFactors = FALSE)
      if (isTRUE(verdict$keep) && (is.null(best) ||
          .ctFitOptimisedLoglik(refit) > .ctFitOptimisedLoglik(best$refit)))
        best <- list(cand = cand, refit = refit, rebuilt = rebuilt, row = length(history))
    }
    if (is.null(best)) { refused <- c(refused, f$key); next }
    history[[best$row]]$kept <- TRUE
    say(label, "kept ", best$cand$label)
    current <- best$refit; kept <- c(kept, list(best$cand)); changed <- TRUE
    data <- best$rebuilt$data; args <- best$rebuilt$args
  }
  ranked <- .ctImproveRankStep(base, current, kept, startdata, covariates, criterion,
    alpha, cores, say, label = label, round = round)
  if (ranked$changed) {
    current <- ranked$current; kept <- ranked$kept; data <- ranked$data
    args <- ranked$args; changed <- TRUE
  }
  list(current = current, kept = kept, data = data, args = args,
    history = c(history, ranked$history), changed = changed)
}

# The findings of the kinds in `describe`, in a sentence for a warning, the
# covariance naming two of its parameters, and the call that tries the changes
# for them. Random effects held at zero or barely determined have their own
# lines (`.ctIdentifyAdvice()`, `.ctEffectAdvice()`), so they are not
# described here. '' where no finding has a change to try.
.ctReduceWarnText <- function(findings,
    describe = c('unread', 'drift-intercept', 'drift-diffusion', 'covariance')) {
  shown <- Filter(function(f) f$kind %in% describe, findings)
  if (!length(shown)) return('')
  parts <- vapply(shown, function(f) if (identical(f$kind, 'covariance'))
    paste0(f$text, ' (', paste(utils::head(f$parameters, 2), collapse = ', '),
      if (length(f$parameters) > 2) ', ...' else '', ')') else f$text, '')
  text <- paste(parts, collapse = '; ')
  paste0(toupper(substring(text, 1, 1)), substring(text, 2), '. ',
    'ctModelImprove(fit, reduceOnly = TRUE) refits the changes that remove these ',
    'and keeps those the data allow.')
}
