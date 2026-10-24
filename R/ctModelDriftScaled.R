# `ctModel(driftScaled = c('DIFFUSION', 'CINT'))`: a process's diffusion sd and
# continuous intercept parameterised by what they accrue over one time unit,
# given the process's own drift diagonal `a`, rather than by their
# instantaneous rates.
#
#   DIFFUSION[i,i] = s * sqrt(2 a / expm1(2 a)),  s^2 the variance accrued
#   CINT[i]        = u * a / expm1(a),            u the level accrued
#
# `s^2` is the variance the process gains over one time unit from its own
# noise: the stationary variance where the process is fast (a -> -inf), the
# diffusion variance itself where it is slow (a -> 0), and defined for a > 0.
# `u` likewise moves from the stationary mean contribution to the intercept.
# The likelihood is the same as the default parameterisation's; what changes
# is that these parameters no longer have to move with the drift diagonal
# along the direction where faster mean reversion and more noise look alike
# -- the white-noise ridge, which the data at coarse intervals cannot close --
# so that direction lies along the drift alone. Only the diffusion's diagonal
# is rescaled: its correlations are still the continuous-time diffusion's, so
# fixing one to zero means what it always has. One time unit, because that is
# the unit the data and every other parameter are already in.
#
# Built from cells the model already supports, so both backends fit it as it
# stands: each free parameter moves, with its name, transform, random effect
# and predictor effects, into a new PARS cell, and its cell becomes an
# expression of that PARS cell and DRIFT[i,i]. Not where DRIFT[i,i] is fixed at
# zero (the factor is then one), nor for a cell that is already fixed or an
# expression. Exactly a = 0 is 0/0 in the expression; the default drift
# transform stops short of it.

# FALSE in a predictor-effect column's own type, which may be text.
.ctDriftScaledFalse <- function(x) if (is.character(x)) rep("FALSE", length(x)) else
  rep(FALSE, length(x))

# The rewrite of a continuous-time model's parameter table.
.ctModelDriftScaled <- function(model, which) {
  which <- unique(as.character(which))
  unknown <- setdiff(which, c("DIFFUSION", "CINT"))
  if (length(unknown)) stop("driftScaled names matrices to rescale by the drift ",
    "diagonal: 'DIFFUSION' and 'CINT'. Not ", paste(unknown, collapse = ", "), ".",
    call. = FALSE)
  if (!length(which)) return(model)
  if (!isTRUE(model$continuoustime)) stop("driftScaled applies to continuous-time ",
    "models (type = 'ct'): in discrete time the diffusion and intercept are ",
    "already per interval.", call. = FALSE)
  pars <- model$pars
  effectcols <- intersect(paste0(model$TIpredNames, "_effect"), names(pars))
  plain <- function(r) !is.na(pars$param[r]) && is.na(pars$value[r]) &&
    grepl("^[A-Za-z][A-Za-z0-9_.]*$", pars$param[r])
  driftzero <- function(i) {
    r <- which(pars$matrix == "DRIFT" & pars$row == i & pars$col == i)
    length(r) == 1L && is.na(pars$param[r]) && isTRUE(as.numeric(pars$value[r]) == 0)
  }
  parsrows <- which(pars$matrix == "PARS")
  nextrow <- if (length(parsrows)) max(pars$row[parsrows]) else 0L
  ncolpars <- if (length(parsrows)) max(pars$col[parsrows]) else 1L
  added <- list()
  for (mat in which) for (i in seq_len(model$n.latent)) {
    r <- which(pars$matrix == mat & pars$row == i & pars$col == (if (mat == "CINT") 1L else i))
    if (length(r) != 1L || !plain(r) || driftzero(i)) next
    nextrow <- nextrow + 1L
    moved <- pars[r, , drop = FALSE]
    moved$matrix <- "PARS"
    moved$row <- nextrow
    moved$col <- 1L
    # The rest of the new PARS row, fixed at zero, so the matrix has every cell.
    pad <- pars[rep(r, ncolpars - 1L), , drop = FALSE]
    if (nrow(pad)) {
      pad$matrix <- "PARS"
      pad$row <- nextrow
      pad$col <- seq_len(ncolpars)[-1L]
      pad$param <- NA_character_
      pad$transform <- NA_character_
      pad$value <- 0
      pad$indvarying <- FALSE
      for (cc in effectcols) pad[[cc]] <- .ctDriftScaledFalse(pad[[cc]])
    }
    factor <- if (mat == "DIFFUSION")
      sprintf("sqrt(2 * DRIFT[%d,%d] / expm1(2 * DRIFT[%d,%d]))", i, i, i, i) else
      sprintf("DRIFT[%d,%d] / expm1(DRIFT[%d,%d])", i, i, i, i)
    pars$param[r] <- sprintf("PARS[%d,1] * %s", nextrow, factor)
    pars$transform[r] <- NA_character_
    pars$value[r] <- NA_real_
    pars$indvarying[r] <- FALSE
    for (cc in effectcols) pars[[cc]][r] <- .ctDriftScaledFalse(pars[[cc]][r])
    added[[length(added) + 1L]] <- rbind(moved, pad)
  }
  if (length(added)) {
    pars <- rbind(pars, do.call(rbind, added))
    rownames(pars) <- NULL
  }
  model$pars <- pars
  model$driftScaled <- which
  .ctModelRawPopVarSync(model)
}
