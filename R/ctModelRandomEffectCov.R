# The population covariance, as ordinary T0VAR cells.
#
# NOT WIRED IN. `ctModel()` does not call this yet, and the two paragraphs at
# the end of this comment say what has to be true before it can. The code below
# works and is tested by hand; what does not work is the rest of the package's
# reaction to it.
#
# A varying parameter is a state under the augmented layout: its population mean
# becomes a T0MEANS row and its population covariance an entry of T0VAR. That
# has always been true, and it happened inside `ctStanModelIntOverPop()` on the
# way to a backend, so the model a user held never showed it. The consequence
# was that there was nowhere to *say* how large a population spread should be --
# the cells that hold the answer did not exist until after the point where
# anyone could edit them.
#
# So `ctModel()` extends T0VAR itself. The extra rows and columns are ordinary
# cells: free with a label by default, and a value or a different label may be
# written into them like any other matrix entry.
#
#   m$matrices$T0VAR[3, 3] <- 0.3   # population sd of the third dimension
#   m$matrices$T0VAR[3, 1] <- 0     # uncorrelated with the first
#
# The labelling follows what T0VAR already uses for latents -- `T0var_<name>` on
# the diagonal, `T0var_<row>_<col>` for a correlation below it, and a fixed zero
# above -- with the parameter's own name standing in for a latent's.
#
# This makes every model with individual differences visibly larger, which is
# the point: `indvarying` defaults to TRUE for T0MEANS, MANIFESTMEANS and CINT,
# so a great many models are multilevel without their author having thought
# about it, and the population covariance they are implicitly estimating was
# invisible until the fit had already made its assumptions.

# The parameters that add a dimension to T0VAR.
#
# T0MEANS parameters do not: their dimension is the latent itself, which T0VAR
# already covers. Everything else varying gets one, in the order
# `ctStanModelIntOverPop()` assigns states, so the two agree about which
# dimension is which.
#' @keywords internal
.ctModelRandomEffectNames <- function(pars) {
  if (is.null(pars$indvarying)) return(character())
  varying <- !is.na(pars$indvarying) & pars$indvarying &
    !is.na(pars$param)
  if (!any(varying)) return(character())
  t0means <- as.character(pars$param[varying & pars$matrix %in% "T0MEANS"])
  unique(setdiff(as.character(pars$param[varying]), t0means))
}

# Default contents for one new T0VAR cell, in T0VAR's own idiom.
#' @keywords internal
.ctModelRandomEffectCell <- function(rowname, colname, continuoustime) {
  defaults <- ctStanModelDefaultFreePar("T0VAR",
    row = if (identical(rowname, colname)) 1L else 2L,
    col = 1L, continuoustime = continuoustime)
  list(
    param = if (identical(rowname, colname)) paste0("T0var_", rowname) else
      paste0("T0var_", rowname, "_", colname),
    transform = defaults$transform)
}

# Extend T0VAR to cover the random effects, leaving anything already there
# alone.
#' @keywords internal
.ctModelExtendT0VAR <- function(pars, n.latent, latentNames, continuoustime,
  tieffects = character()) {
  extra <- .ctModelRandomEffectNames(pars)
  if (!length(extra)) return(pars)
  names <- c(latentNames, extra)
  size <- length(names)
  t0var <- which(pars$matrix %in% "T0VAR")
  if (!length(t0var)) return(pars)
  template <- pars[t0var[1L], , drop = FALSE]
  existing <- paste(pars$row[t0var], pars$col[t0var], sep = "_")

  additions <- list()
  for (row in seq_len(size)) {
    for (col in seq_len(size)) {
      # Only the cells the extension actually adds. A model round-tripped
      # through this function twice must not gain anything the second time, and
      # a user's own edits must survive it.
      if (row <= n.latent && col <= n.latent) next
      if (paste(row, col, sep = "_") %in% existing) next
      cell <- template
      cell$matrix <- "T0VAR"
      cell$row <- row
      cell$col <- col
      cell$indvarying <- FALSE
      if (length(tieffects)) cell[, tieffects] <- "FALSE"
      if (col > row) {
        # Above the diagonal T0VAR is a fixed zero, as it already is for the
        # latent block: the matrix is built from a lower factor.
        cell$param <- NA_character_
        cell$value <- 0
        cell$transform <- NA_character_
        cell$sdscale <- NA_real_
      } else {
        spec <- .ctModelRandomEffectCell(names[row], names[col], continuoustime)
        cell$param <- spec$param
        cell$value <- NA_real_
        cell$transform <- spec$transform
        cell$sdscale <- 1
      }
      additions[[length(additions) + 1L]] <- cell
    }
  }
  if (!length(additions)) return(pars)
  rbind(pars, do.call(rbind, additions))
}

# ---------------------------------------------------------------------------
# Why this is not called yet
#
# Extending T0VAR in `ctModel()` works -- a default two-latent two-manifest
# model grows a 4x4 T0VAR labelled `T0var_mm_Y1`, `T0var_mm_Y2_mm_Y1` and so
# on, values and labels can be written into it, and re-extending is a no-op.
# What breaks is everything downstream, for one reason:
#
#   **the population covariance is not an ordinary free parameter.**
#
# Stan parameterises it through `rawpopcovbase`/`rawpopsd` and the julia
# backend builds `julia_popcov_*` entries in its prepared table. Adding
# labelled T0VAR cells therefore counts it a second time: a one-latent model
# with one varying parameter went from 4 free parameters to 7, and the Stan
# round trip failed with "Number of unconstrained parameters does not match
# that of the model (8 vs 23)". Fourteen Laplace tests and eleven
# nonlinear-reporting tests failed with it.
#
# The second obstacle is subtler and rules out the obvious fix. Making the
# backends *honour a fixed cell instead of relabelling it* cannot be done by
# looking at the cell: `ctStanModelIntOverPop()` creates these placeholders
# with `value = 0`, and a user fixing a population covariance to zero writes
# exactly the same thing. Testing `!is.na(value)` pinned every population
# variance in every model to zero.
#
# So the cells need a marker of their own -- a `populationcov` column, say --
# saying "this is the population covariance, in its own parameterisation" so
# that:
#
#   * the Stan writer skips them when enumerating free parameters, since it
#     builds T0cov from rawpopcov and overwrites these cells anyway;
#   * the julia augmentation relabels an unmarked-as-set cell and honours one
#     the user actually set;
#   * a placeholder is distinguishable from a deliberate zero.
#
# Note also that the marker cannot be positional. For a varying T0MEANS,
# T0VAR[1,1] is *both* the initial-state variance and that parameter's
# population variance -- the same cell, two meanings -- so "row or column
# beyond n.latent" does not identify the population block.
