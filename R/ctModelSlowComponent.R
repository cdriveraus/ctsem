#' Add slow processes that drive a model's latent processes
#'
#' Each latent in \code{indices} gains a slow latent, named \code{slow_} and
#' the latent's name, that drives it through the drift: an estimated effect
#' \code{DRIFT[k, slow_k]} moves the level the latent returns to with the slow
#' process, so a process can wander over weeks while it fluctuates over hours.
#' Its equilibrium is \code{DRIFT[k, slow_k] / |DRIFT[k, k]|} times the slow
#' process, so a person with faster dynamics follows the slow process less
#' closely.
#'
#' A slow process is slower than the latent it drives by construction: its
#' drift is the latent's times \code{inv_logit(r)}, with \code{r} a
#' \code{PARS} parameter named \code{drift_} and the slow process's name, so
#' its rate is a fraction of the latent's. Left free, the two rates can trade
#' places: with only the latent measured, a latent of drift \code{a} driven by
#' a slow process of drift \code{b} has the spectrum of one with the two
#' exchanged and the effect adjusted, so the "slow" process can come out the
#' fast one, and refits from other starting values did. The latent's drift
#' parameter is declared in \code{PARS} for this, keeping its name,
#' transform, individual variation and predictor effects, so where the
#' latent's rate varies over people the slow process's follows it. Where the
#' latent's drift is fixed the fraction is of that value; where it is an
#' expression the slow process's persistence is left free, as before.
#'
#' A slow process has innovations of unit variance (the effect carries its
#' scale), correlated with those of the other slow processes, including ones
#' added earlier, and not with the fast ones.
#' It has no loadings and no intercept, and no individual variation or
#' predictor effects. Its initial mean and variance are its stationary ones
#' (\code{T0MEANS} and the \code{T0VAR} diagonal \code{'stationary'}, see
#' \code{\link{ctFit}}): zero and \code{1 / (2 |a|)} for persistence
#' \code{a}. Its initial covariance with the latent it drives is free -- that
#' latent follows it, so the two start related -- and with the other states
#' zero.
#'
#' The effect has to start away from zero. At zero neither the slow process nor
#' its effect has a gradient, so an optimiser started there stays. When
#' refitting from a model without the slow process, start the effect at one
#' (the slow process's sign is not identified, so either sign will do) and
#' the fraction small; \code{\link{ctModelImprove}} starts it at a tenth.
#'
#' @param ctm A model from \code{ctModel(type = 'ct')}.
#' @param indices Which latents get a slow process, by position.
#' @return The extended model.
#' @seealso \code{\link{ctModelHigherOrder}} for the other structural extension
#'   the search tries.
#' @export
#' @examples
#' m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
#'   latentNames = c('eta1', 'eta2'), LAMBDA = diag(2))
#' m2 <- ctModelSlowComponent(m, 1:2)
#' m2$latentNames
#' p <- m2$pars
#' p[p$matrix %in% c('DRIFT', 'DIFFUSION') & (p$row > 2 | p$col > 2) &
#'   (!is.na(p$param) | p$value != 0), c('matrix', 'row', 'col', 'param', 'value')]
ctModelSlowComponent <- function(ctm, indices = seq_len(ctm$n.latent)) {
  if (!inherits(ctm, 'ctStanModel')) stop("ctModelSlowComponent() takes a model ",
    "from ctModel(type = 'ct').", call. = FALSE)
  if (!isTRUE(ctm$continuoustime)) stop("Slow processes are for continuous time ",
    "models.", call. = FALSE)
  nl <- ctm$n.latent
  indices <- as.integer(indices)
  if (!length(indices) || anyDuplicated(indices) || any(indices < 1L | indices > nl))
    stop("indices must be distinct latents, from 1 to ", nl, ".", call. = FALSE)
  names0 <- ctm$latentNames
  newnames <- paste0('slow_', names0[indices])
  already <- newnames[newnames %in% names0]
  if (length(already)) stop("The model already has ", paste(already, collapse = ', '),
    ".", call. = FALSE)
  ed <- .ctParsLatentEditor(ctm, newnames)
  allnames <- ed$names
  slow <- c(grep('^slow_', names0), nl + seq_along(indices))
  for (k in seq_along(indices)) {
    i <- indices[k]; n <- nl + k
    .ctSlowDrift(ed, i, n, paste0('drift_', allnames[n]))
    ed$free('DRIFT', i, n, paste0('drift_', names0[i], '_', allnames[n]), tipred = FALSE)
    ed$fixed('DIFFUSION', n, n, 1)
    for (s in slow[slow < n]) ed$free('DIFFUSION', n, s,
      paste0('diff_', allnames[n], '_', allnames[s]), tipred = FALSE)
    ed$stationary('T0MEANS', n, 1)
    ed$free('T0VAR', n, i, paste0('T0var_', allnames[n], '_', names0[i]), tipred = FALSE)
    ed$stationary('T0VAR', n, n)
  }
  ed$finish()
}

# The slow process's drift as a fraction of latent i's: DRIFT[n, n] =
# <latent i's drift> * inv_logit(<name>), `name` a PARS parameter on the raw
# scale. A free drift of latent i is declared in PARS under its own name so
# the cell can read it; a fixed one enters as its value; an expression, or no
# negative value, leaves the slow drift free and unordered.
.ctSlowDrift <- function(ed, i, n, name) {
  rec <- ed$cell('DRIFT', i, i)
  latent <- NA_character_
  if (!is.null(rec)) {
    if (is.na(rec$param)) {
      v <- suppressWarnings(as.numeric(rec$value))
      if (is.finite(v) && v < 0) latent <- sprintf('(%s)', format(v, digits = 15))
    } else if (grepl('^[A-Za-z][A-Za-z0-9_.]*$', rec$param) && rec$param != 'stationary') {
      latent <- rec$param
      ed$declare('DRIFT', i, i)
    }
  }
  if (is.na(latent)) {
    ed$free('DRIFT', n, n, name, tipred = FALSE)
    return(invisible())
  }
  ed$free('PARS', ed$npars() + 1L, 1L, name, transform = 'param', tipred = FALSE)
  ed$expression('DRIFT', n, n, sprintf('%s * inv_logit(%s)', latent, name))
  invisible()
}
