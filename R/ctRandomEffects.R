# Each group's conditional posterior of its random effects, stored on a julia
# fit as `fit$randomeffects`.
#
# Every person-level report needs the same thing: where a person's (or a
# group's) effects are, given their data, and how well their data determine
# them. On the augmented route the effects are filter states and the smoother
# gives both; on the Laplace route they are the inner modes and the inverse of
# the curvature there. Both come out of the one engine pass the fit already
# makes at its estimate to measure how much of each effect the data determine
# (`ctsem_effect_information`, R/ctBackendEffectInformation.R), so storing them
# costs nothing, and a later report reads them rather than rebuilding the
# objective and re-solving every mode.
#
# One list entry per level with random effects:
#   level  the id column the groups are of;
#   param  the raw parameters that vary there, and `parnumber` their
#          positions in the raw parameter vector;
#   ids    the groups, labelled as in the data, in the fit's order;
#   mean   groups x effects: the deviation of each group's raw parameter from
#          its population value (and TI-predictor effects) at the mode -- `L z`
#          on the Laplace route, the smoothed effect state's departure from
#          where the filter started it on the augmented one, both in raw units;
#   cov    groups x k(k+1)/2: the conditional covariance of those deviations,
#          packed column by column from the diagonal down
#          (`.ctRandomEffectsCov()` unpacks it);
#   route, and `evaluated_at`, the raw vector they were computed at, so a
#   reader can tell they belong to the estimate it reports.
#
# On the Laplace route the posterior is the Gaussian the approximation
# integrates: exact for effects that enter linearly with Gaussian measurement,
# an approximation otherwise. A level's own variances integrated in the filter
# (`popunique` with the hybrid) are coordinates of the layer this pass runs on,
# so they are here with the rest.

# Row and column of each packed entry, in packing order.
.ctPackedIndex <- function(k) {
  if (k < 1L) return(matrix(integer(), 0L, 2L))
  cols <- unlist(lapply(seq_len(k), function(j) rep(j, k - j + 1L)))
  rows <- unlist(lapply(seq_len(k), function(j) j:k))
  cbind(row = rows, col = cols)
}

# One group's covariance from a stored entry, or every group's as a
# groups x k x k array when `group` is NULL.
.ctRandomEffectsCov <- function(entry, group = NULL) {
  k <- length(entry$param)
  index <- .ctPackedIndex(k)
  unpack <- function(values) {
    m <- matrix(NA_real_, k, k, dimnames = list(entry$param, entry$param))
    m[index] <- values
    m[index[, 2:1, drop = FALSE]] <- values
    m
  }
  if (!is.null(group)) {
    if (is.character(group)) group <- match(group, as.character(entry$ids))
    return(unpack(entry$cov[group, ]))
  }
  out <- array(NA_real_, c(nrow(entry$cov), k, k),
    dimnames = list(as.character(entry$ids), entry$param, entry$param))
  for (g in seq_len(nrow(entry$cov))) out[g, , ] <- unpack(entry$cov[g, ])
  out
}

# The engine's per-level posterior, labelled for the fit: groups by their ids
# in the data rather than the engine's numbering.
.ctRandomEffectsStore <- function(fit, posterior) {
  if (is.null(posterior) || !length(posterior)) return(NULL)
  spec <- .ctBackendSpec(fit)
  subjectid <- .ctFitModelObject(fit)$subjectIDname
  ids <- as.character(.ctFitIdMap(fit)$original)
  data <- spec$data
  # Each varying parameter's position in the raw vector, by its label, the same
  # way on both routes: the augmented route's effect rows number their sds.
  cells <- .ctBackendFreeParameterCells(fit)
  cells <- cells[!cells$randomeffect, , drop = FALSE]
  labels_of <- .ctBackendParameterNames(cells)
  out <- lapply(posterior, function(entry) {
    parnumber <- as.integer(cells$parnumber[match(entry$param, labels_of)])
    if (anyNA(parnumber) && length(entry$parnumber) == length(parnumber))
      parnumber[is.na(parnumber)] <- entry$parnumber[is.na(parnumber)]
    k <- length(entry$param)
    ngroups <- nrow(entry$mean)
    level <- if (length(entry$level)) entry$level else subjectid
    if (!is.null(entry$first_member)) {
      # A sampler's blocks name their groups by their first member.
      member <- if (identical(entry$index, 1L) || is.null(data[[level]])) ids else
        as.character(data[[level]][match(ids, as.character(data[[subjectid]]))])
      labels <- member[entry$first_member]
    } else if (identical(entry$index, 1L) || !length(level) || is.null(data[[level]])) {
      labels <- if (length(ids) == ngroups) ids else as.character(seq_len(ngroups))
    } else {
      # The engine numbers a level's groups; each subject's group at that
      # level, in the fit's subject order, names them.
      member <- as.character(data[[level]][match(ids,
        as.character(data[[subjectid]]))])
      labels <- rep(NA_character_, ngroups)
      first <- !duplicated(entry$group)
      labels[entry$group[first]] <- member[first]
    }
    mean <- entry$mean
    dimnames(mean) <- list(labels, entry$param)
    index <- .ctPackedIndex(k)
    cov <- entry$cov
    dimnames(cov) <- list(labels, paste0(entry$param[index[, 1L]], ",",
      entry$param[index[, 2L]]))
    list(level = level, param = entry$param, parnumber = parnumber,
      ids = labels, mean = mean, cov = cov)
  })
  names(out) <- vapply(out, `[[`, "", "level")
  structure(out, route = if (!is.null(spec$laplace)) "laplace" else "augmented",
    evaluated_at = as.numeric(fit$estimate$raw))
}


# Each subject's own raw parameter vector at the estimate on the augmented
# route: the population values, its TI-predictor effects and its effects'
# posterior mean. What the Laplace route's engine builds as
# `ctsem_laplace_subject_values`; here the pieces are on the fit already.
# Missing predictor values count as the mean (zero, centred), and effects fixed
# to a value (`ti_fixed`) are left out.
.ctRandomEffectsSubjectRaw <- function(fit, entry) {
  spec <- .ctBackendSpec(fit)
  # With any stated population means after the free parameters
  # (R/ctFixedPopMeans.R), as a subject's vector carries them.
  estimate <- .ctBackendRawFull(fit, as.numeric(.ctFitRawEstimate(fit)))
  n <- nrow(entry$mean)
  raw <- matrix(estimate, n, length(estimate), byrow = TRUE)
  effects <- spec$ti_effects
  predictors <- spec$tipred_data
  if (!is.null(effects) && NROW(effects) && !is.null(predictors) &&
      NROW(predictors) == n) {
    effects <- as.data.frame(effects, stringsAsFactors = FALSE)
    predictors <- as.matrix(predictors)
    predictors[!is.finite(predictors)] <- 0
    for (e in seq_len(nrow(effects))) {
      j <- as.integer(effects$parameter[e])
      raw[, j] <- raw[, j] + estimate[as.integer(effects$coefficient[e])] *
        predictors[, as.integer(effects$predictor[e])]
    }
  }
  raw[, entry$parnumber] <- raw[, entry$parnumber] + entry$mean
  raw
}

# The augmented route's draws with each subject's own uncertainty added
# (`ctSubjectPars(uncertainty = 'total')`). There the draws come from the
# filter's subject matrices at each population draw, so the population part is
# exact and there is no raw vector per draw to perturb. Each subject's
# conditional draw `e ~ N(0, cov)` from the stored posterior is added on the
# parameter's own scale as `g(r + e) - g(r)`, `r` the subject's raw vector at
# the estimate: exact for a linear transform, and for a nonlinear one the
# spread is that of the transform at the estimate rather than at each draw.
# `values` is draws x subjects x the varying parameters `cells` names, in
# that order.
.ctRandomEffectsAugmentedDraws <- function(fit, values, cells) {
  entry <- fit$randomeffects[[1L]]
  if (is.null(entry) || is.null(entry$parnumber)) stop("This fit carries no ",
    "posterior of its subjects' effects (fits made before ctsem stored one); ",
    "refit it, or use uncertainty = 'population'.", call. = FALSE)
  position <- match(cells$parnumber, entry$parnumber)
  if (anyNA(position)) stop("Some varying parameters have no stored posterior ",
    "on this route; use uncertainty = 'population'.", call. = FALSE)
  ndraws <- dim(values)[1L]
  n <- dim(values)[2L]
  k <- length(entry$param)
  base <- .ctRandomEffectsSubjectRaw(fit, entry)
  covs <- .ctRandomEffectsCov(entry)
  noise <- array(0, c(ndraws, n, k))
  for (i in seq_len(n)) {
    factor <- tryCatch(chol(covs[i, , ]), error = function(e) NULL)
    if (is.null(factor)) next
    noise[, i, ] <- matrix(stats::rnorm(ndraws * k), ndraws, k) %*% factor
  }
  # Subject by subject, each subject's draws together, as `stacked` is.
  stacked <- base[rep(seq_len(n), each = ndraws), , drop = FALSE]
  shift <- do.call(rbind, lapply(seq_len(n), function(i)
    matrix(noise[, i, ], ndraws, k)))
  stacked[, entry$parnumber] <- stacked[, entry$parnumber] + shift
  layout <- .ctBackendSummaryLayout(fit)
  perturbed <- .ctBackendPopCellValues(fit, stacked, cells, layout)
  centre <- .ctBackendPopCellValues(fit, base, cells, layout)
  for (i in seq_len(n)) {
    rows <- (i - 1L) * ndraws + seq_len(ndraws)
    values[, i, ] <- values[, i, ] + sweep(perturbed[rows, , drop = FALSE], 2L,
      centre[i, ], "-")
  }
  values
}


# A sampled fit's per-group posterior (`fit$randomeffects`). Sampling the joint
# posterior draws every group's effects, and the samplers keep each group's
# mean and covariance over the draws (`effect_mean`, `effect_cov`, pooled over
# processes) whether or not the draws themselves are saved. A marginal sample
# integrated the effects out, so there the posterior is the one an optimised
# fit stores, at the posterior mean.
.ctRandomEffectsFromSample <- function(fit) {
  smp <- fit$sample
  spec <- .ctBackendSpec(fit)
  if (!length(smp$effect_mean) || !length(smp$effect_cov)) {
    check <- .ctEffectInformation(fit$model_spec, fit$estimate$raw,
      point = "at the posterior mean")
    if (is.null(check)) return(NULL)
    return(.ctRandomEffectsStore(fit, check$posterior))
  }
  layout <- .ctJuliaGet(.ctJuliaModule(spec$project)$ctsem_laplace_deviation_layout(
    .ctJuliaObjective(fit)))
  level <- as.integer(layout$level)
  within <- as.integer(layout$within)
  first <- as.integer(layout$first_member)
  starts <- which(within == 1L)
  structure <- .ctSpecRandomEffectLevels(fit$model_spec)
  width <- vapply(structure, function(s) length(s$params), integer(1))
  packed <- width[level[starts]] * (width[level[starts]] + 1L) / 2L
  covstart <- c(0L, cumsum(packed))[seq_along(starts)]
  mean <- as.numeric(smp$effect_mean)
  cov <- as.numeric(smp$effect_cov)
  if (length(mean) != length(level) || length(cov) != sum(packed)) return(NULL)
  posterior <- list()
  for (l in seq_along(structure)) {
    blocks <- which(level[starts] == l)
    k <- width[l]
    if (!length(blocks) || !k) next
    m <- matrix(vapply(blocks, function(b) mean[starts[b] + seq_len(k) - 1L],
      numeric(k)), ncol = k, byrow = TRUE)
    q <- k * (k + 1L) / 2L
    v <- matrix(vapply(blocks, function(b) cov[covstart[b] + seq_len(q)],
      numeric(q)), ncol = q, byrow = TRUE)
    posterior[[length(posterior) + 1L]] <- list(
      level = as.character(spec$laplace$levels[[l]]$name), index = l,
      param = as.character(structure[[l]]$params), mean = m, cov = v,
      first_member = first[starts[blocks]])
  }
  out <- .ctRandomEffectsStore(fit, posterior)
  if (!is.null(out)) attr(out, "route") <- "sampled"
  out
}
