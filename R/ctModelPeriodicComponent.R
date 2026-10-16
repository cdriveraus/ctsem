#' Add a cycle, such as a daily rhythm, that drives a model's latent processes
#'
#' Adds a pair of latents, \code{<name>_cos} and \code{<name>_sin}, that turn
#' at a fixed \code{period} with no noise of their own -- \code{DRIFT[cos, sin]
#' = -2 pi / period} and \code{DRIFT[sin, cos] = 2 pi / period} -- so they
#' trace the cosine and sine of the time of day (or week, or whatever the
#' period is). Each latent in \code{indices} is driven by both, through
#' estimated effects \code{DRIFT[k, <name>_cos]} and \code{DRIFT[k,
#' <name>_sin]}, so it is pushed through the cycle with an amplitude and phase
#' the data decide, and responds with its own dynamics: a fast process follows
#' the cycle closely, a slow one smooths and delays it.
#'
#' The cycle is set to the clock at each subject's first observation, from
#' two time independent predictors, \code{<name>_cos0} and \code{<name>_sin0},
#' the cosine and sine of that time over the period, which
#' \code{\link{ctPeriodicPhases}} adds to the data. So the data's time must
#' count from an origin shared by every subject -- hours since midnight of a
#' common day, say -- for the cycle to mean the same clock time for all.
#' With \code{individual = TRUE} each subject's cycle may also be shifted from
#' the clock and differ in strength, by an initial deviation of the pair with
#' an estimated standard deviation, which the filter learns from each
#' subject's data: rhythms that differ between people, or a time origin that
#' is not shared, are then still modelled.
#'
#' The effects start at zero, where they have a gradient: unlike a slow
#' process, the cycle is known, so a refit from the model without it needs no
#' particular start.
#'
#' Needs the julia backend, which reads time independent predictors in the
#' initial means.
#'
#' @param ctm A model from \code{ctModel(type = 'ct')}.
#' @param period The cycle's period, in the data's time units: 24 for a daily
#'   rhythm with time in hours.
#' @param indices Which latents the cycle drives, by position.
#' @param name Prefix for the added latents, predictors and parameters.
#' @param individual Let each subject's cycle deviate from the clock.
#' @return The extended model.
#' @seealso \code{\link{ctPeriodicPhases}} for the predictors it reads;
#'   \code{\link{ctFitAdequacy}}'s \code{periodic} family, which finds cycles a
#'   model lacks; \code{\link{ctModelImprove}}, which adds one where that
#'   family points.
#' @export
#' @examples
#' m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
#'   latentNames = c('eta1', 'eta2'), LAMBDA = diag(2))
#' m2 <- ctModelPeriodicComponent(m, period = 24)
#' m2$latentNames
#' m2$TIpredNames
#' p <- m2$pars
#' p[p$matrix %in% 'DRIFT' & (p$row > 2 | p$col > 2) &
#'   (!is.na(p$param) | p$value != 0), c('matrix', 'row', 'col', 'param', 'value')]
ctModelPeriodicComponent <- function(ctm, period, indices = seq_len(ctm$n.latent),
  name = 'cycle', individual = TRUE) {
  if (!inherits(ctm, 'ctStanModel')) stop("ctModelPeriodicComponent() takes a model ",
    "from ctModel(type = 'ct').", call. = FALSE)
  if (!isTRUE(ctm$continuoustime)) stop("A cycle is added to continuous time models.",
    call. = FALSE)
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0)
    stop("period is one positive number, in the data's time units.", call. = FALSE)
  nl <- ctm$n.latent
  indices <- as.integer(indices)
  if (!length(indices) || anyDuplicated(indices) || any(indices < 1L | indices > nl))
    stop("indices must be distinct latents, from 1 to ", nl, ".", call. = FALSE)
  newnames <- paste0(name, c('_cos', '_sin'))
  preds <- paste0(name, c('_cos0', '_sin0'))
  taken <- c(newnames[newnames %in% ctm$latentNames], preds[preds %in% ctm$TIpredNames])
  if (length(taken)) stop("The model already has ", paste(taken, collapse = ', '),
    "; give the cycle another name.", call. = FALSE)
  # The phase predictors are read by the initial means alone: no effects, on
  # these cells or on any a later change adds (`statedTIpreds`).
  for (p in preds) ctm$pars[[paste0(p, '_effect')]] <- 'FALSE'
  ctm$TIpredNames <- c(ctm$TIpredNames, preds)
  ctm$n.TIpred <- length(ctm$TIpredNames)
  ctm$statedTIpreds <- c(ctm$statedTIpreds, preds)
  ed <- .ctParsLatentEditor(ctm, newnames)
  allnames <- ed$names
  cs <- nl + 1L; sn <- nl + 2L
  w <- 2 * pi / period
  ed$fixed('DRIFT', cs, sn, -w)
  ed$fixed('DRIFT', sn, cs, w)
  for (i in indices) for (j in c(cs, sn)) ed$free('DRIFT', i, j,
    paste0('drift_', allnames[i], '_', allnames[j]), tipred = FALSE)
  # A cell naming a predictor in an expression takes the subject's value.
  ed$fixed('T0MEANS', cs, 1, 0); ed$fixed('T0MEANS', sn, 1, 0)
  t0 <- ed$finish()
  for (k in 1:2) {
    r <- which(t0$pars$matrix == 'T0MEANS' & t0$pars$row == nl + k)
    t0$pars$param[r] <- paste0(preds[k], ' * 1')
    t0$pars$value[r] <- NA
    t0$pars$transform[r] <- NA
  }
  if (isTRUE(individual)) {
    ed2 <- .ctParsLatentEditor(t0, character())
    for (j in c(cs, sn)) ed2$free('T0VAR', j, j, paste0(name, '_T0sd'), tipred = FALSE)
    t0 <- ed2$finish()
  }
  t0
}

#' The phase predictors a cycle added by ctModelPeriodicComponent reads
#'
#' Adds to long format data, for each subject, the cosine and sine of the
#' time of its first row over \code{period}: the columns \code{<name>_cos0}
#' and \code{<name>_sin0} that \code{\link{ctModelPeriodicComponent}} sets its
#' cycle from.
#'
#' @param datalong Long format data, one row per observation.
#' @param period The cycle's period, in the data's time units.
#' @param name The cycle's name, as given to
#'   \code{\link{ctModelPeriodicComponent}}.
#' @param id,time The subject and time columns.
#' @return \code{datalong} with the two columns added.
#' @export
#' @examples
#' d <- data.frame(id = rep(1:2, each = 3), time = c(8, 12, 20, 33, 37, 45))
#' ctPeriodicPhases(d, period = 24)
ctPeriodicPhases <- function(datalong, period, name = 'cycle', id = 'id', time = 'time') {
  datalong <- as.data.frame(datalong)
  if (!all(c(id, time) %in% names(datalong))) stop("datalong needs the '", id, "' and '",
    time, "' columns.", call. = FALSE)
  t0 <- stats::ave(as.numeric(datalong[[time]]), datalong[[id]], FUN = function(t)
    rep(min(t, na.rm = TRUE), length(t)))
  datalong[[paste0(name, '_cos0')]] <- cos(2 * pi * t0 / period)
  datalong[[paste0(name, '_sin0')]] <- sin(2 * pi * t0 / period)
  datalong
}

# The periods of the cycles a model already has: each clock pair
# ctModelPeriodicComponent() adds, <name>_cos and <name>_sin, turning at
# DRIFT[<name>_sin, <name>_cos] = 2 pi / period.
.ctModelCyclePeriods <- function(model) {
  lat <- model$latentNames
  out <- numeric()
  for (c in which(endsWith(lat, '_cos'))) {
    s <- match(paste0(sub('_cos$', '', lat[c]), '_sin'), lat)
    if (is.na(s)) next
    r <- model$pars[model$pars$matrix == 'DRIFT' & model$pars$row == s &
      model$pars$col == c, , drop = FALSE]
    w <- suppressWarnings(as.numeric(r$value))
    if (length(w) == 1L && is.finite(w) && w > 0) out <- c(out, 2 * pi / w)
  }
  out
}
