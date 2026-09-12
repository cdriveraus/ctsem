# The matsetup vocabulary.
#
# `matsetup` is built once by `.ctModelMatSetup()` (R/ctModelWriter.R) and read
# from R, from the two generated Stan programs, and from the julia spec builder.
# The Stan side reads it positionally, and CLAUDE.md's standing advice for that
# is to grep `matsetup[`. The R side reads it by *name*, so that grep finds none
# of it -- and the question "which row is a free population parameter" was
# written out inline in eleven places across five files, in seven spellings that
# did not agree. Some included `copyrow < 1` and some did not; one used
# `when %in% c(0,100)` where the others used `c(0,-1)`.
#
# `when` is an enum, assigned at R/ctModelWriter.R:577-629:
#
#   0     static parameter (the default)
#   1     t0 matrices, materialised in the second t0 pass
#   2     dynamic (DRIFT, CINT, DIFFUSION, JAx)
#   3     tdpred
#   4     measurement
#   100   PARS -- a wildcard meaning "materialise at every when"
#   -999  never used during calculations
#   -1    a population row carried alongside the above
#
# The 100 wildcard has already cost one silent wrong answer. The post-mortem is
# in the comment above `whenvecp` in `ctModelWriter.R`: it let a PARS carrier row
# answer a parameter lookup, and the parameters that lost their own transform
# lost a drift diagonal's negative definiteness and a diffusion diagonal's
# positivity and variance floor. Nothing errored.
#
# So this function does not hide the wildcard -- it makes a caller that wants it
# say so, by name, at the call site. `when` is an argument rather than a
# constant for exactly that reason.

#' Rows of matsetup that materialise a free population parameter
#'
#' The shared predicate behind every "which parameters are there" question asked
#' of `matsetup`. Returns a logical vector the length of `nrow(ms)`.
#'
#' @param ms matsetup, as a data.frame.
#' @param when which `when` values count. The default `c(0, -1)` is the
#'   population rows. Pass `c(0, 100)` to include the PARS wildcard -- and only
#'   where the wildcard is genuinely wanted, since a PARS row's number is not a
#'   parameter number in the same sense.
#' @param defining if TRUE, keep only the row that defines a parameter, not a
#'   row that copies it (`copyrow < 1`).
#' @param varying if TRUE, keep only parameters that vary over subjects, by a
#'   TI predictor or a random effect.
#' @param free if TRUE (the default), require `param > 0`. `nparams` is the one
#'   caller that wants the numbering range rather than the free rows, and it
#'   takes `free = FALSE` because a fixed row still carries a valid number.
#'
#' @return logical vector, `nrow(ms)` long.
#' @keywords internal
.ctMatsetupFreeRows <- function(ms, when = c(0, -1), defining = FALSE,
  varying = FALSE, free = TRUE) {
  keep <- ms$when %in% when
  if (free) keep <- keep & ms$param > 0
  if (defining) keep <- keep & ms$copyrow < 1
  if (varying) keep <- keep & (ms$tipred > 0 | ms$indvarying > 0)
  keep
}

#' The free population parameters of a model, one row each, in parameter order
#'
#' `.ctMatsetupFreeRows()` with the deduplication and ordering every caller of
#' this shape also wanted.
#'
#' @param m matsetup, as a data.frame.
#' @param ... passed to [.ctMatsetupFreeRows()].
#' @return `m`, subset to one row per parameter and ordered by parameter number.
#' @keywords internal
ctMatsetupFreePars <- function(m, ...) {
  m <- m[.ctMatsetupFreeRows(m, ...), , drop = FALSE]
  m <- m[match(unique(m$param), m$param), , drop = FALSE]
  m[order(m$param), , drop = FALSE]
}
