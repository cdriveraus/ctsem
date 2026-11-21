# Reducing a model whose fit leaves directions the data do not identify.
#
# A fit's identifiability report (`fit$identifiability$directions`) lists the
# directions along which the likelihood is flat, each with the parameters it
# is made of. A few kinds recur, and each has a model change that removes it:
#
#   unread        a parameter the likelihood does not read (a measurement
#                 error covariance between items scored one at a time, before
#                 that was refused, or an auto-effect at its transform's bound
#                 near zero): curvature exactly zero along it alone. Fixed at
#                 zero where zero is a model (an off-diagonal, an effect, an
#                 auto-effect).
#   drift-intercept / drift-diffusion
#                 a latent's auto-effect against its intercept or its
#                 diffusion, with their spread over subjects: the process
#                 settles faster than the data resolve, so only their ratio is
#                 determined. A latent without noise (a trend) is refitted at
#                 both ends of its auto-effect, as a linear trend (auto-effect
#                 fixed at zero) and as a constant level (auto-effect and
#                 intercept fixed at zero); any is refitted with its intercept
#                 or diffusion scaled by the drift (`driftScaled`), which
#                 states the same model along the ratio.
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
# criterion, the one with fewest parameters where several are kept; the full
# search does the same before it screens and after each change it keeps. On follhim/dis_ctsem (2026-11-21) the original fit's report
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

# The effect a population coordinate is the spread of: its sd, own variance or
# loadings; '' for a correlation, which is the spread of no one effect. A level
# above the subject is a suffix (`popsd_<parameter>.<level>`).
.ctReduceCoordinateEffect <- function(p, model) {
  levels <- model$groupIDnames
  if (length(levels)) p <- sub(paste0('\\.(', paste(levels, collapse = '|'), ')$'), '', p)
  if (grepl('^poploading_.+_dim[0-9]+$', p)) return(sub('^poploading_(.+)_dim[0-9]+$', '\\1', p))
  if (grepl('^(popsd|popunique)_', p)) return(sub('^(popsd|popunique)_', '', p))
  ''
}

# A latent's auto-effect fixed at zero: a linear trend where the latent has no
# noise of its own, a random walk where it has. One key for both readings, so a
# change refused under one finding is not refitted under another.
.ctReduceDriftZero <- function(model, drift, i) .ctReduceProposal(
  paste0('reduce:drift0:', drift), paste0(drift, ' fixed at 0 (',
    if (.ctReduceNoNoise(model, i)) 'a linear trend' else 'a random walk', ')'),
  .ctReduceFixEdit(drift, 0))

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
  # The model as written; while ctFit() is still assembling the fit, which is
  # when its warnings are made, the one the engine ran (`.ctFitModelObject()`),
  # whose parameters and cells are the same.
  model <- .ctFitBaseModel(fit)
  if (is.null(model)) model <- tryCatch(.ctFitModelObject(fit), error = function(e) NULL)
  id <- fit$identifiability
  if (is.null(model) || is.null(id)) return(list())
  zero <- .ctEffectRows(id$effects, 'zero')
  weak <- .ctEffectRows(id$effects, 'weak')
  held <- c(if (!is.null(zero)) zero$coordinate, if (!is.null(weak)) weak$coordinate)
  out <- list()
  add <- function(f) if (!f$key %in% vapply(out, `[[`, '', 'key')) out[[length(out) + 1L]] <<- f
  latent <- function(i) model$latentNames[i]
  zeroown <- character(); zerosd <- character()
  spreadof <- character(); rest <- list()
  for (d in id$directions) {
    # A ridge of the augmented route between a population sd and its
    # correlations (`$partial`) is that route's, which the Laplace route
    # identifies (`.ctIdentifyAdvice()`): no change to the model.
    if (!is.null(d$partial)) next
    p <- setdiff(as.character(d$parameters), held)
    if (!length(p)) next
    # Own variances or sds alone at zero curvature: held at the floor of
    # their transforms, where every smaller value is as good. Several at once
    # span a flat space whose basis mixes them, so a direction of only such
    # coordinates names each. On follhim/dis_ctsem two own variances of a
    # rank-5 covariance sat there, and at rank 3 two more shared directions.
    if (isTRUE(abs(d$relative) < 1e-12) && all(grepl('^(popunique|popsd)_[^.]+$', p))) {
      zeroown <- unique(c(zeroown, sub('^popunique_', '', grep('^popunique_', p, value = TRUE))))
      zerosd <- unique(c(zerosd, sub('^popsd_', '', grep('^popsd_', p, value = TRUE))))
      next
    }
    pop <- .ctReducePopCoordinate(p)
    mp <- p[!pop]
    # The effect each population coordinate is the spread of; a correlation
    # names two, and is the spread of neither alone.
    own <- vapply(p[pop], .ctReduceCoordinateEffect, '', model = model)
    cells <- .ctReduceCells(model, mp)
    if (length(p) == 1L && isTRUE(abs(d$relative) < 1e-12) && !is.na(cells$matrix[1L])) {
      # Zero curvature, to roundoff, along one parameter alone: nothing the
      # likelihood computes reads it (`.ctOptimIdentifiedInverse()`'s 1e-12 for
      # a direction with no width at all), or it sits where its transform is
      # flat. An off-diagonal, an effect or an auto-effect is fixed at zero (an
      # auto-effect at its bound near zero is a linear trend, or a random walk
      # where the latent has noise); any other cell is named, not changed.
      cand <- list()
      if (cells$matrix[1L] %in% 'DRIFT' && cells$row[1L] == cells$col[1L]) {
        cand <- list(.ctReduceDriftZero(model, mp, cells$row[1L]))
      } else if (cells$row[1L] != cells$col[1L] ||
          cells$matrix[1L] %in% c('TDPREDEFFECT', 'LAMBDA')) {
        cand <- list(.ctReduceProposal(paste0('reduce:fix:', mp),
          paste0(mp, ' fixed at 0'), .ctReduceFixEdit(mp, 0)))
      }
      add(list(key = paste0('unread:', mp), kind = 'unread', parameters = mp,
        text = paste0(mp, ' is not read by the likelihood'), candidates = cand))
      next
    }
    # A latent's auto-effect with its intercept, or with its diffusion: the
    # spread of the intercept or diffusion over subjects scales with it too, so
    # its own population coordinates belong to the same finding. On the dis
    # model reduced to rank 3, the stress trend's one flat direction was its
    # auto-effect, its intercept and the intercept's three loadings and own
    # variance -- two thirds population coordinates, and no covariance problem.
    explained <- character()
    for (i in unique(cells$row[cells$matrix %in% 'DRIFT' & cells$row == cells$col])) {
      drift <- cells$param[cells$matrix %in% 'DRIFT' & cells$row == i & cells$col == i][1L]
      cint <- .ctReduceCellParam(model, 'CINT', i, 1L)
      diffusion <- .ctReduceCellParam(model, 'DIFFUSION', i, i)
      spread <- function(x) names(own)[own %in% x]
      if (!is.na(cint) && (cint %in% mp || length(spread(cint)))) {
        spreadof <- c(spreadof, cint)
        explained <- c(explained, drift, cint, spread(cint))
        cand <- list()
        if (.ctReduceNoNoise(model, i)) cand <- list(.ctReduceDriftZero(model, drift, i),
          .ctReduceProposal(paste0('reduce:level:', drift), paste0(latent(i),
            ' a constant level (', drift, ' and ', cint, ' fixed at 0)'),
            .ctReduceFixEdit(c(drift, cint), 0)))
        if (!'CINT' %in% model$driftScaled) cand[[length(cand) + 1L]] <- .ctReduceProposal(
          'reduce:driftscaled:CINT', "intercepts scaled by the drift (driftScaled = 'CINT')",
          .ctReduceDriftScaledEdit('CINT'))
        add(list(key = paste0('driftcint:', i), kind = 'drift-intercept',
          parameters = c(drift, cint), text = paste0(drift, ' and ', cint,
            ' trade off: ', latent(i), ' settles faster than the data resolve'),
          candidates = cand))
      }
      if (!is.na(diffusion) && (diffusion %in% mp || length(spread(diffusion)))) {
        spreadof <- c(spreadof, diffusion)
        explained <- c(explained, drift, diffusion, spread(diffusion))
        add(list(key = paste0('driftdiffusion:', i), kind = 'drift-diffusion',
          parameters = c(drift, diffusion), text = paste0(drift, ' and ', diffusion,
            ' trade off: ', latent(i), ' forgets its state faster than the data resolve'),
          candidates = if ('DIFFUSION' %in% model$driftScaled) list() else
            list(.ctReduceProposal('reduce:driftscaled:DIFFUSION',
              "diffusions scaled by the drift (driftScaled = 'DIFFUSION')",
              .ctReduceDriftScaledEdit('DIFFUSION')))))
      }
    }
    rest[[length(rest) + 1L]] <- list(p = p, rest = setdiff(p, explained),
      explained = length(explained) > 0L)
  }
  # What no pattern accounts for: mostly population coordinates is a
  # covariance with more dimensions than the subjects identify, which the rank
  # step reduces; anything else is named. The flat directions of one fit span a
  # space whose basis is arbitrary, so one direction can carry two problems (a
  # trend's trade-off beside a covariance's surplus dimensions), and the spread
  # of an intercept in a trade-off can lie in a direction of its own: on the dis
  # model at rank 3, the intercept's loadings and own variance alone. So the
  # spread of an effect in any trade-off is the trade-off's, in every direction.
  for (r in rest) {
    left <- r$rest[!vapply(r$rest, function(x) .ctReduceCoordinateEffect(x, model) %in%
      spreadof, logical(1))]
    if (length(left) && mean(.ctReducePopCoordinate(left)) >= 0.5) {
      add(list(key = 'covariance', kind = 'covariance',
        parameters = left[.ctReducePopCoordinate(left)],
        text = 'the population covariance has more dimensions than the subjects identify',
        candidates = list()))
    } else if (!r$explained && length(left)) add(list(key = paste0('other:',
      paste(sort(r$p), collapse = ',')), kind = 'other', parameters = r$p,
      text = paste(utils::head(r$p, 4), collapse = ', '), candidates = list()))
  }
  single <- !length(model$groupIDnames)
  # Own variances at zero leave the level's own variances, the rank kept: a
  # rank change (`.ctImproveRankChange()`) that replaces the level's.
  if (length(zeroown)) {
    own <- .ctReduceOwnVariances(fit, model)
    cand <- list()
    if (single && !is.null(own) && isTRUE(own$rank >= 1L) && all(zeroown %in% own$own)) {
      rest <- setdiff(own$own, zeroown)
      ch <- .ctImproveRankChange(own$rank, if (length(rest)) rest else FALSE)
      ch$label <- paste0('no own variance for ', paste(zeroown, collapse = ', '))
      cand <- list(ch)
    }
    add(list(key = paste0('zeroown:', paste(sort(zeroown), collapse = ',')),
      kind = 'own variance at zero', parameters = paste0('popunique_', zeroown),
      text = paste0('the own variance', if (length(zeroown) > 1L) 's', ' of ',
        paste(zeroown, collapse = ', '), if (length(zeroown) > 1L) ' are' else ' is',
        ' at zero'), candidates = cand))
  }
  if (length(zerosd) && single) add(list(key = paste0('zerosd:',
      paste(sort(zerosd), collapse = ',')), kind = 'sd at zero', parameters = zerosd,
    text = paste0(paste(zerosd, collapse = ', '), ': the data hold their individual ',
      'differences at zero'), candidates = list(.ctReduceProposal(paste0('reduce:fixed:',
      paste(sort(zerosd), collapse = ',')), paste0(paste(zerosd, collapse = ', '),
      ' not varying over subjects'), .ctReduceFixedEffectsEdit(zerosd)))))
  # Random effects the data hold at zero, or barely determine: on a model with
  # one level only, where `indvarying = FALSE` means what it says.
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

# Which effects have an own variance on a fit's level, and its rank, from the
# call the fit replays: `popunique` TRUE is every varying parameter, names are
# parameters or matrices. NULL where the call cannot be read or is by level.
.ctReduceOwnVariances <- function(fit, model) {
  replay <- tryCatch(.ctFitReplayArgs(fit), error = function(e) NULL)
  if (is.null(replay)) return(NULL)
  pu <- replay$popunique
  rank <- suppressWarnings(as.integer(replay$poprank))
  pars <- model$pars
  varying <- unique(pars$param[as.logical(pars$indvarying) %in% TRUE & !is.na(pars$param)])
  own <- if (isTRUE(pu)) varying else if (is.character(pu) && is.null(names(pu)))
    varying[varying %in% pu | pars$matrix[match(varying, pars$param)] %in% pu] else
      return(NULL)
  list(rank = if (length(rank) == 1L) rank else NA_integer_, own = own)
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
  # Said, not refused: a reduced fit that has not reached its maximum only
  # understates what the reduction gains.
  saddle <- if (isTRUE(reduced$identifiability$negative > 0L)) '; not at a maximum' else ''
  if (k1 >= k0) return(c(out, list(keep = ll1 >= ll0 - tol, value = NA_real_,
    note = paste0(if (ll1 >= ll0 - tol) 'restates the model' else
      sprintf('restated, %.3g below', ll0 - ll1), saddle))))
  if (ll1 > ll0 + tol) return(c(out, list(keep = TRUE, value = NA_real_,
    note = paste0(sprintf('%.3g better with fewer parameters', ll1 - ll0), saddle))))
  v <- .ctImproveCompare(reduced, current, criterion, alpha, cores)
  if (!is.finite(v$dloglik)) return(fail(v$note))
  c(out, list(keep = !isTRUE(v$accept), value = v$value, note = paste0(v$note, saddle)))
}

# The reduction stage: while the current fit's report has a finding with a
# change not yet refused, refit each of its changes, keep the best the
# criterion allows, and read the new fit's report (`.ctReduceFindingsStep()`);
# then lower the population covariance's rank while it is not identified
# (`.ctImproveRankStep()`), and where that keeps a rank, look again: a lower
# rank can leave own variances at zero, which the findings drop. Every change
# reduces, so this ends. Returns as `.ctImproveRankStep()` does.
.ctImproveReduce <- function(base, current, kept, startdata, covariates, criterion,
    alpha, cores, say = function(...) invisible(), label = '', round = NA_integer_) {
  st <- list(current = current, kept = kept, data = NULL, args = NULL, history = list(),
    changed = FALSE, tried = character())
  repeat {
    st <- .ctReduceFindingsStep(st, base, startdata, covariates, criterion, alpha,
      cores, say, label, round)
    ranked <- .ctImproveRankStep(base, st$current, st$kept, startdata, covariates,
      criterion, alpha, cores, say, label = label, round = round)
    st$history <- c(st$history, ranked$history)
    if (!ranked$changed) break
    st[c('current', 'kept', 'data', 'args')] <- ranked[c('current', 'kept', 'data', 'args')]
    st$changed <- TRUE
  }
  st[c('current', 'kept', 'data', 'args', 'history', 'changed')]
}

# The findings part of a reduction stage, on its state `st`. A change refused
# once is not refitted, under the same finding or another that proposes it
# (`st$tried`, by key and label).
.ctReduceFindingsStep <- function(st, base, startdata, covariates, criterion, alpha,
    cores, say, label, round) {
  id <- function(cand) paste0(cand$key, '|', cand$label)
  repeat {
    findings <- lapply(.ctReduceFindings(st$current), function(f) {
      f$candidates <- Filter(function(cand) !id(cand) %in% st$tried, f$candidates)
      f
    })
    findings <- Filter(function(f) length(f$candidates), findings)
    if (!length(findings)) return(st)
    f <- findings[[1L]]
    say(label, f$text)
    best <- NULL
    for (cand in f$candidates) {
      say(label, "refitting with ", cand$label)
      # A rank change replaces the level's (as `.ctImproveRankStep()` does).
      prior <- if (identical(cand$family, 'rank'))
        Filter(function(p) !identical(p$family, 'rank'), st$kept) else st$kept
      rebuilt <- try(.ctImproveRebuild(base, startdata, c(prior, list(cand)), covariates),
        silent = TRUE)
      refit <- if (inherits(rebuilt, 'try-error')) rebuilt else try(.ctImproveRefit(
        st$current, rebuilt$model, rebuilt$data, cores, NULL, rebuilt$args), silent = TRUE)
      verdict <- .ctReduceVerdict(st$current, refit, criterion, alpha, cores)
      st$history[[length(st$history) + 1L]] <- data.frame(round = round,
        family = 'reduce', change = cand$label, screen_p = NA_real_,
        dloglik = verdict$dloglik, dpar = verdict$dpar, criterion = verdict$value,
        kept = FALSE, note = verdict$note, stringsAsFactors = FALSE)
      st$tried <- c(st$tried, id(cand))
      # Of the changes kept, the one with fewest parameters, then the higher
      # likelihood: a restatement only where nothing smaller is allowed.
      better <- is.null(best) || verdict$dpar < best$dpar || (verdict$dpar == best$dpar &&
        .ctFitOptimisedLoglik(refit) > .ctFitOptimisedLoglik(best$refit))
      if (isTRUE(verdict$keep) && isTRUE(better))
        best <- list(cand = cand, refit = refit, rebuilt = rebuilt,
          row = length(st$history), prior = prior, dpar = verdict$dpar)
    }
    if (is.null(best)) next
    st$history[[best$row]]$kept <- TRUE
    say(label, "kept ", best$cand$label)
    st$current <- best$refit
    st$kept <- c(best$prior, list(best$cand))
    st$data <- best$rebuilt$data
    st$args <- best$rebuilt$args
    st$changed <- TRUE
  }
}

# The findings of the kinds in `describe`, in a sentence for a warning, the
# covariance naming two of its parameters, and the call that tries the changes
# for them. Random effects held at zero or barely determined have their own
# lines (`.ctIdentifyAdvice()`, `.ctEffectAdvice()`), so they are not
# described here. '' where no finding has a change to try.
.ctReduceWarnText <- function(findings,
    describe = c('unread', 'drift-intercept', 'drift-diffusion', 'covariance',
      'own variance at zero', 'sd at zero')) {
  shown <- Filter(function(f) f$kind %in% describe, findings)
  if (!length(shown)) return('')
  parts <- vapply(shown, function(f) if (identical(f$kind, 'covariance'))
    paste0(f$text, ' (', paste(utils::head(f$parameters, 2), collapse = ', '),
      if (length(f$parameters) > 2) ', ...' else '', ')') else f$text, '')
  # Capitalised only where it starts with a word: a part can start with a
  # parameter's name, which stays as it is.
  text <- paste(parts, collapse = '; ')
  if (startsWith(text, 'the ')) text <- paste0('T', substring(text, 2))
  paste0(text, '. ',
    'ctModelImprove(fit, reduceOnly = TRUE) refits the changes that remove these ',
    'and keeps those the data allow.')
}
