# Reducing a model whose fit leaves directions the data do not identify.
#
# Nothing here knows what a model's parameters mean. The fit's identifiability
# report (`fit$identifiability$directions`) gives the directions along which
# the likelihood is flat, and each parameter's transform gives its range. From
# those alone:
#
#   A coordinate is undetermined where most of its variation lies in the flat
#   subspace: its share, the squared length of its axis projected on an
#   orthonormal basis of the directions, is at least a half. An undetermined
#   model parameter is refitted at the simple values its range offers (each
#   finite end, and zero where zero is inside) and, where none of those is
#   kept, at its estimate, which states the same fit with one parameter fewer.
#   An undetermined population sd or own variance is refitted without it, zero
#   being the end of a variance's range.
#
#   A flat subspace that lies mostly in the other population coordinates is a
#   covariance with more dimensions than the subjects identify, even where no
#   one correlation or loading is undetermined: refitted at lower ranks
#   (`.ctImproveRankStep()`), from full rank with own variances (`popunique`).
#
#   Random effects the effects report holds at zero or barely determines are
#   made fixed.
#
# `ctModelImprove(fit, reduceOnly = TRUE)` refits each change from the current
# estimate, keeps it unless the fuller model is better by the search's
# criterion (of several kept, the one with fewest parameters), and reads the
# new fit's report again; the full search does the same before it screens and
# after each change it keeps.
#
# Rules written for the ridges one model showed do not carry to the next. On
# follhim/dis_ctsem (2026-11-21) a rule for an auto-effect against its
# intercept fired on the original fit by coincidence: that auto-effect had
# 0.003 of its variation in the flat subspace, which lay in the covariance's
# correlations (none with more than 0.28 of its own there). It became
# undetermined only once the rank was lowered, and the shares then said so.

# A raw coordinate of a population covariance.
.ctReducePopCoordinate <- function(p) grepl('^(popsd_|rawcor_|poploading_|popunique_)', p)

# A population sd or own variance at the subject level.
.ctReduceSpreadCoordinate <- function(p) grepl('^(popsd|popunique)_[^.]+$', p)

# The cells of model parameters, by name: matrix, row and column of the cell
# that holds each (`.ctParamRows()`), NA where the name is not a parameter of
# the model's table (a predictor effect, a population coordinate).
.ctReduceCells <- function(model, params) {
  rows <- .ctParamRows(model, params)
  pars <- model$pars
  data.frame(param = params, matrix = pars$matrix[rows], row = pars$row[rows],
    col = pars$col[rows], stringsAsFactors = FALSE)
}

# A reduction as `.ctImproveRebuild()` applies it: an edit of the model, never
# pruned (undoing it adds back what the data did not identify). A `fallback`
# is tried only where no other change for its finding was kept.
.ctReduceProposal <- function(key, label, edit, fallback = FALSE) list(key = key,
  family = 'reduce', label = label, p = NA_real_, rank = 0, prunable = FALSE,
  edit = edit, fallback = fallback)

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

# A varying parameter's population mean stated at `value`
# (`model$fixedpopmeans`), its random effect kept.
.ctReduceStateEdit <- function(param, value) {
  force(param); force(value)
  function(model, data, covs, args = list()) {
    means <- model$fixedpopmeans
    model$fixedpopmeans <- c(means[setdiff(names(means), param)],
      stats::setNames(value, param))
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

# The flat subspace of an identifiability report: each coordinate's share of
# it, and the coordinates with zero curvature alone, which nothing the
# likelihood computes reads. The share is taken on an orthonormal basis of the
# directions' loadings, since the report lists a direction twice where two
# eigenvalues are equal to roundoff (two own variances at zero on the dis model
# came as four directions) and keeps only the loadings above
# `.ctNullMassBar()`. Coordinates in `held` (random effects the effects report
# has a line for) are left out.
.ctReduceFlat <- function(id, held = character()) {
  dirs <- Filter(function(d) is.null(d$partial), id$directions)
  params <- setdiff(unique(unlist(lapply(dirs, function(d) as.character(d$parameters)))),
    held)
  out <- list(share = stats::setNames(numeric(), character()), unread = character())
  if (!length(params)) return(out)
  V <- matrix(0, length(params), length(dirs), dimnames = list(params, NULL))
  for (j in seq_along(dirs)) {
    p <- as.character(dirs[[j]]$parameters)
    l <- as.numeric(dirs[[j]]$loadings)
    if (length(l) != length(p) || !all(is.finite(l))) l <- rep(1, length(p))
    keep <- p %in% params
    V[p[keep], j] <- l[keep]
  }
  q <- qr(V, tol = 1e-6)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  out$share <- stats::setNames(pmin(rowSums(Q^2), 1), params)
  out$unread <- unique(unlist(lapply(dirs, function(d) {
    p <- setdiff(as.character(d$parameters), held)
    if (length(p) == 1L && isTRUE(abs(d$relative) < 1e-12)) p
  })))
  out
}

# How much of the flat subspace lies in population coordinates other than the
# sds and own variances undetermined on their own: at least half a dimension
# is a covariance with more dimensions than the subjects identify.
.ctReduceCovarianceMass <- function(share) {
  pop <- .ctReducePopCoordinate(names(share))
  alone <- .ctReduceSpreadCoordinate(names(share)) & share >= 0.5
  sum(share[pop & !alone])
}

# The values a parameter's range offers as simpler models, and its estimate,
# on the parameter's own scale: its cells evaluated at the estimate with its
# raw coordinate moved far out each way (`.ctBackendPopMeanSamples()`). An end
# is finite where the transform has stopped moving there; one within 1e-4 of
# zero is zero, since a transform keeps a variance or an auto-effect that far
# from it only to stay inside its range. Zero is offered too where it lies
# inside the range. A named list by parameter, each with `ends`, `zero` and
# `estimate`; parameters it cannot evaluate are left out.
.ctReduceValues <- function(fit, params) {
  raw <- .ctFitNameRaw(fit, fit$estimate$raw)
  index <- match(params, names(raw))
  params <- params[!is.na(index)]; index <- index[!is.na(index)]
  if (!length(params)) return(list())
  far <- c(-40, -30, 30, 40)
  rows <- do.call(rbind, lapply(index, function(k) {
    m <- matrix(raw, nrow = 1L + length(far), ncol = length(raw), byrow = TRUE)
    m[-1L, k] <- far
    m
  }))
  mapped <- tryCatch(.ctBackendPopMeanSamples(fit, samples = rows), error = function(e) NULL)
  if (is.null(mapped) || is.null(mapped$values)) return(list())
  values <- as.matrix(mapped$values)
  out <- list()
  for (i in seq_along(params)) {
    column <- match(index[i], as.integer(mapped$parnumber))
    if (is.na(column)) next
    v <- values[(i - 1L) * (1L + length(far)) + seq_len(1L + length(far)), column]
    still <- function(a, b) is.finite(a) && is.finite(b) &&
      abs(a - b) <= 1e-6 * max(1, abs(a))
    ends <- numeric()
    if (still(v[2], v[3])) ends <- c(ends, v[2])
    if (still(v[5], v[4])) ends <- c(ends, v[5])
    ends <- unique(ifelse(abs(ends) < 1e-4, 0, signif(ends, 6)))
    lo <- min(v[2], v[5]); hi <- max(v[2], v[5])
    out[[params[i]]] <- list(ends = ends,
      zero = is.finite(lo) && is.finite(hi) && lo < 0 && hi > 0 && !0 %in% ends,
      estimate = signif(v[1], 6))
  }
  out
}

# The changes for undetermined model parameters: each fixed at its range's
# simple values, then, as the fallback, the first (the largest share) at its
# estimate. A parameter varying over subjects keeps its random effect where
# its mean is stated inside its range; at an end, where no mean can be stated,
# it is fixed and not varying.
.ctReduceParameterChanges <- function(fit, model, params) {
  values <- .ctReduceValues(fit, params)
  varying <- .ctImproveVarying(model)
  num <- function(x) format(x, digits = 3)
  out <- list()
  for (p in names(values)) {
    v <- values[[p]]
    for (x in v$ends) out[[length(out) + 1L]] <- .ctReduceProposal(
      paste0('reduce:fix:', p, '=', x), paste0(p, ' fixed at ', num(x),
        if (p %in% varying) ', not varying'), .ctReduceFixEdit(p, x))
    if (isTRUE(v$zero)) out[[length(out) + 1L]] <- if (p %in% varying)
      .ctReduceProposal(paste0('reduce:mean:', p, '=0'), paste0(p, "'s mean fixed at 0"),
        .ctReduceStateEdit(p, 0)) else
      .ctReduceProposal(paste0('reduce:fix:', p, '=0'), paste0(p, ' fixed at 0'),
        .ctReduceFixEdit(p, 0))
  }
  first <- intersect(params, names(values))[1L]
  if (!is.na(first)) {
    x <- values[[first]]$estimate
    out[[length(out) + 1L]] <- if (first %in% varying)
      .ctReduceProposal(paste0('reduce:mean:', first, '=', x), paste0(first,
        "'s mean fixed at its estimate, ", num(x)), .ctReduceStateEdit(first, x),
        fallback = TRUE) else
      .ctReduceProposal(paste0('reduce:fix:', first, '=', x), paste0(first,
        ' fixed at its estimate, ', num(x)), .ctReduceFixEdit(first, x), fallback = TRUE)
  }
  out
}

# `x is` or `x, y are`, for a sentence.
.ctReduceAre <- function(x) paste(paste(x, collapse = ', '),
  if (length(x) > 1L) 'are' else 'is')

# What a fit's identifiability report says, sorted by the rules above into
# findings, each with the coordinates involved, a phrase saying what it is,
# whether `ctModelImprove(fit, reduceOnly = TRUE)` has a change for it
# (`reducible`), and, with `candidates`, those changes (which evaluates the
# parameters' ranges, so a warning leaves them out).
.ctReduceFindings <- function(fit, candidates = TRUE) {
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
  flat <- .ctReduceFlat(id, held)
  share <- flat$share
  coords <- names(share)
  single <- !length(model$groupIDnames)
  out <- list()
  finding <- function(key, kind, parameters, text, cand = list(), reducible = TRUE)
    out[[length(out) + 1L]] <<- list(key = key, kind = kind, parameters = parameters,
      text = text, candidates = cand, reducible = reducible)

  pop <- .ctReducePopCoordinate(coords)
  undetermined <- share >= 0.5 | coords %in% flat$unread
  incell <- !is.na(.ctReduceCells(model, coords)$matrix)
  mp <- coords[!pop & incell & undetermined]
  mp <- mp[order(-share[mp], match(mp, coords))]
  if (length(mp)) {
    unread <- intersect(mp, flat$unread)
    rest <- setdiff(mp, unread)
    finding('parameters', 'parameters', mp, paste(c(
      if (length(unread)) paste(.ctReduceAre(unread), 'not read by the likelihood'),
      if (length(rest)) paste(.ctReduceAre(rest), 'not determined by the data')),
      collapse = '; '), if (candidates) .ctReduceParameterChanges(fit, model, mp))
  }
  # Undetermined coordinates no change reaches: a predictor effect, a cell
  # written as an expression.
  other <- coords[!pop & !incell & undetermined]
  if (length(other)) finding(paste0('other:', paste(sort(other), collapse = ',')),
    'other', other, paste(.ctReduceAre(other), 'not determined by the data'),
    reducible = FALSE)

  # Undetermined own variances leave the level's own variances, the rank kept:
  # a rank change (`.ctImproveRankChange()`) that replaces the level's. All of
  # them at once, then, as the fallback, the one with the largest share.
  spread <- coords[.ctReduceSpreadCoordinate(coords) & undetermined]
  spread <- spread[order(-share[spread])]
  own <- sub('^popunique_', '', grep('^popunique_', spread, value = TRUE))
  if (length(own)) {
    info <- .ctReduceOwnVariances(fit, model)
    able <- single && !is.null(info) && isTRUE(info$rank >= 1L) && all(own %in% info$own)
    cand <- list()
    if (candidates && able) for (drop in unique(list(own, own[1L]))) {
      rest <- setdiff(info$own, drop)
      ch <- .ctImproveRankChange(info$rank, if (length(rest)) rest else FALSE)
      ch$label <- paste0('no own variance for ', paste(drop, collapse = ', '))
      ch$fallback <- length(cand) > 0L
      cand[[length(cand) + 1L]] <- ch
    }
    finding(paste0('own:', paste(sort(own), collapse = ',')), 'own variances',
      paste0('popunique_', own), paste0('the data do not determine the own variance',
        if (length(own) > 1L) 's', ' of ', paste(own, collapse = ', ')), cand,
      reducible = able)
  }
  sds <- sub('^popsd_', '', grep('^popsd_', spread, value = TRUE))
  if (length(sds)) {
    cand <- list()
    if (candidates && single) for (drop in unique(list(sds, sds[1L])))
      cand[[length(cand) + 1L]] <- .ctReduceProposal(paste0('reduce:fixed:',
        paste(sort(drop), collapse = ',')), paste0(paste(drop, collapse = ', '),
        ' not varying over subjects'), .ctReduceFixedEffectsEdit(drop),
        fallback = length(cand) > 0L)
    finding(paste0('sds:', paste(sort(sds), collapse = ',')), 'sds',
      paste0('popsd_', sds), paste0('the data do not determine how much ',
        paste(sds, collapse = ', '), if (length(sds) > 1L) ' vary' else ' varies',
        ' over subjects'), cand, reducible = single)
  }
  # The covariance: the rank step's (`.ctImproveWeakCovariance()`).
  if (.ctReduceCovarianceMass(share) >= 0.5) {
    covs <- coords[pop & !(coords %in% spread)]
    finding('covariance', 'covariance', covs[order(-share[covs])],
      'the population covariance has more dimensions than the subjects identify',
      reducible = single)
  }
  # Random effects the data hold at zero, or barely determine: on a model with
  # one level only, where `indvarying = FALSE` means what it says.
  for (w in list(list(rows = zero, kind = 'no spread',
      text = ': the data hold their individual differences at zero'),
    list(rows = weak, kind = 'weak', text = ': the data barely determine their individual differences'))) {
    if (is.null(w$rows) || !nrow(w$rows) || !single) next
    eff <- unique(as.character(w$rows$effect))
    finding(paste0(w$kind, ':', paste(sort(eff), collapse = ',')), w$kind, eff,
      paste0(paste(eff, collapse = ', '), w$text),
      if (candidates) list(.ctReduceProposal(paste0('reduce:fixed:', paste(sort(eff),
        collapse = ',')), paste0(paste(eff, collapse = ', '), ' not varying over subjects'),
        .ctReduceFixedEffectsEdit(eff))))
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
# then lower the population covariance's rank by one step if it is not
# identified (`.ctImproveRankStep()`), and where that keeps a rank, look again:
# a lower rank can leave a parameter or an own variance undetermined, which is
# removed before the next rank. On the dis model the stress trend's
# auto-effect became undetermined at rank 5 and walked to raw -760 through
# three more rank refits, each of which took several times as long. Every
# change reduces, so this ends. Returns as `.ctImproveRankStep()` does.
.ctImproveReduce <- function(base, current, kept, startdata, covariates, criterion,
    alpha, cores, say = function(...) invisible(), label = '', round = NA_integer_) {
  st <- list(current = current, kept = kept, data = NULL, args = NULL, history = list(),
    changed = FALSE, tried = character())
  repeat {
    st <- .ctReduceFindingsStep(st, base, startdata, covariates, criterion, alpha,
      cores, say, label, round)
    ranked <- .ctImproveRankStep(base, st$current, st$kept, startdata, covariates,
      criterion, alpha, cores, say, label = label, round = round, once = TRUE)
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
      if (isTRUE(cand$fallback) && !is.null(best)) next
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

# The findings of the kinds in `describe` that the reduction has changes for,
# in a sentence for a warning, the covariance naming two of its coordinates,
# and the call that tries the changes. Random effects held at zero or barely
# determined have their own lines (`.ctIdentifyAdvice()`,
# `.ctEffectAdvice()`), so they are not described here. '' where no finding
# has a change to try.
.ctReduceWarnText <- function(findings,
    describe = c('parameters', 'own variances', 'sds', 'covariance')) {
  shown <- Filter(function(f) f$kind %in% describe && isTRUE(f$reducible), findings)
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
