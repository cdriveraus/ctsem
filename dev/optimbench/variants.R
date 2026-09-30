# Variants: the factor under test, by name. A variant is a list with
#
#   optimcontrol  entries merged into ctFit(optimcontrol = ) for the cell
#   args          other ctFit() arguments, merged over the cell's defaults
#                 (priors, nlcontrol, poprank, ...)
#   objective     a tag for what is being maximised. Cells are scored against
#                 the best end point among cells with the same model, data,
#                 route and objective tag, so a variant that changes the
#                 objective itself (a prior, a floor) must say so here or its
#                 losses are measured against a different function.
#   fit           FALSE for a cell that evaluates at its start and does not fit
#                 (the reference check and the contamination control).
#
# Add a variant here rather than special-casing the harness. A variant is only
# a list of arguments; which code runs is the build under test's business.

BENCH_VARIANTS <- list(
  # The path a user takes: every default of the build under test.
  default = list(),

  # At juliaFit e2abf637 the prior warm-up (`carefulfit`) is skipped whenever
  # `priors` resolves TRUE, and the default 'randomCorr' does
  # (ctJuliaBackend.R warmspec(); plan section 2k item 1). priors = FALSE is
  # the one setting under which it runs. On a model without random-effect
  # correlations 'randomCorr' has no terms, so the objective is unchanged.
  priorsFALSE = list(args = list(priors = FALSE), objective = "noprior"),

  # The warm-up switched off explicitly, for the paired comparison once the
  # warm-up runs by default again.
  nocareful = list(optimcontrol = list(carefulfit = FALSE)),

  # The optimiser alone: no certification, uncertainty or Laplace correction.
  estonly = list(optimcontrol = list(estonly = TRUE)),

  # The overshoot probe off. Paired with `default` on the same cell, what the
  # probe costs; and where no run's probe found a gain, the two must end at
  # the same point, since the probe then changed no verdict.
  noprobe = list(optimcontrol = list(overshoot = "off")),

  # L-BFGS's initial inverse Hessian learned per coordinate (Gilbert &
  # Lemarechal's diagonal update) instead of one secant ratio. Paired with
  # `default` on the same build: whether it can be the default.
  diag = list(optimcontrol = list(lbfgs_diagonal = TRUE)),

  # No fit: evaluate the objective and the references at the start. Used with
  # a stored best-known point to check that a reference still reproduces.
  evalonly = list(fit = FALSE),

  # No fit, no reference: time three gradients at the start. The grid's
  # contamination control.
  control = list(fit = FALSE, control_only = TRUE)
)

bench_variant <- function(name) {
  v <- BENCH_VARIANTS[[name]]
  if (is.null(v)) stop("unknown variant '", name, "'; see variants.R")
  if (is.null(v$optimcontrol)) v$optimcontrol <- list()
  if (is.null(v$args)) v$args <- list()
  if (is.null(v$objective)) v$objective <- "default"
  if (is.null(v$fit)) v$fit <- TRUE
  v$name <- name
  v
}
