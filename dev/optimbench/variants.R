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
  # The scalar initial scale, the default before the diagonal became it.
  scalar = list(optimcontrol = list(lbfgs_diagonal = FALSE)),

  # How many curvature pairs L-BFGS keeps (default 20, chosen on a small
  # model): whether more pays where the curvature is broad and correlated --
  # the growth cells, bigp -- and costs nothing elsewhere.
  mem100 = list(optimcontrol = list(lbfgs_memory = 100L)),
  mem500 = list(optimcontrol = list(lbfgs_memory = 500L)),

  # A structured model of the curvature beside L-BFGS, from the per-unit score
  # outer products (optimcontrol$lbfgs_scores).
  scores = list(optimcontrol = list(lbfgs_scores = TRUE)),

  # The L-BFGS stage's coreset stages on the Laplace route forced on
  # (optimcontrol$coreset; the default, 'auto', leaves them off on every
  # bench cell, all under its 10,000 data cells), and off whatever the
  # build's default: paired on one build, what they save and that the
  # optimum is the same.
  coreset = list(optimcontrol = list(coreset = TRUE)),
  # The local curvature off (`lbfgs_local`, juliaFit localcurv): L-BFGS starts
  # from its learned diagonal, as before it. Paired with `default`, what the
  # exact columns of parameters few units read buy; where no parameter is
  # read by a tenth of the units or fewer, the two run the same code.
  nolocal = list(optimcontrol = list(lbfgs_local = FALSE)),
  nocoreset = list(optimcontrol = list(coreset = FALSE)),
  # The score outer product's model off (`lbfgs_scores`), the default since
  # it was measured to help resumes and refits: whether it still pays once the
  # anchored inner solve made Laplace iterations cheap.
  noscores = list(optimcontrol = list(lbfgs_scores = FALSE)),

  # The first L-BFGS step's length in the metric's norm (`initial_alpha`,
  # default 0.1, set on flat-transform starts): a third of it, three times it,
  # and the full unit step it replaced.
  alpha003 = list(optimcontrol = list(initial_alpha = 0.03)),
  alpha03 = list(optimcontrol = list(initial_alpha = 0.3)),
  alpha1 = list(optimcontrol = list(initial_alpha = 1)),

  # The prior warm-up forced on (the default only where an indicator is not
  # Gaussian), and its cap of 10 iterations at half and twice that.
  careful = list(optimcontrol = list(carefulfit = TRUE)),
  careful5 = list(optimcontrol = list(carefulfit = 5L)),
  careful20 = list(optimcontrol = list(carefulfit = 20L)),

  # An estimate, then what a value, a gradient and a Hessian cost there, and
  # single units against their weights (harness section 6): the measurements
  # one cost model for the pricers is set from. No references.
  costs = list(optimcontrol = list(estonly = TRUE), costs = TRUE),

  # The finish's steps on one curvature whatever the route would choose
  # (`_ctsem_finish_curvature`): a fresh exact Hessian whenever a step is slow,
  # or one Hessian walked with SR1 updates and refreshed when its walk stops
  # paying. Paired with `default`: what the route's own choice costs.
  exact = list(optimcontrol = list(newton = "exact")),
  chord = list(optimcontrol = list(newton = "chord")),

  # The sgd phase (`_ctsem_sgd`) before L-BFGS: optimcontrol$stochastic.
  stoch = list(optimcontrol = list(stochastic = TRUE)),

  # L-BFGS with the diagonal initial inverse Hessian and Grippo-Lampariello-
  # Lucidi acceptance over the last 20 values: the BFGS-side improvement alone.
  gll = list(optimcontrol = list(lbfgs_diagonal = TRUE, lbfgs_gll = 20L)),

  # The whole pipeline: the sgd phase, then that L-BFGS, then the usual finish.
  pipeline = list(optimcontrol = list(stochastic = TRUE, lbfgs_diagonal = TRUE,
    lbfgs_gll = 20L)),

  # SAEM before the optimiser (optimcontrol$saem): the random effects sampled,
  # Kesten's per-parameter steps, the drift rule; then the usual L-BFGS, finish
  # and certification.
  saem = list(optimcontrol = list(saem = TRUE)),
  # SAEM's own estimate: no L-BFGS iterations, no finish, no probe, nothing
  # certified -- what the sampling estimator reaches by itself, scored on the
  # exact reference like every other cell.
  saemonly = list(optimcontrol = list(saem = TRUE, maxiter = 0L, estonly = TRUE,
    overshoot = "off")),
  # The same two with f-SAEM's independence proposals for every subject.
  saemlap = list(optimcontrol = list(saem = TRUE, saem_proposal = "laplace")),
  saemlaponly = list(optimcontrol = list(saem = TRUE, saem_proposal = "laplace",
    maxiter = 0L, estonly = TRUE, overshoot = "off")),
  # Priors on every parameter, so every model is identified: the default
  # optimiser, SAEM before it, and SAEM alone, scored among themselves.
  allprior = list(args = list(priors = TRUE), objective = "allprior"),
  saem_allprior = list(optimcontrol = list(saem = TRUE), args = list(priors = TRUE),
    objective = "allprior"),
  saemonly_allprior = list(optimcontrol = list(saem = TRUE, maxiter = 0L,
    estonly = TRUE, overshoot = "off"), args = list(priors = TRUE),
    objective = "allprior"),

  # The refined fit: the Laplace fit, then the climb on the placed objective
  # with the check's flagged units on their rules (`intoverpop =
  # 'quadrature'`). Its climb is recorded in the lapcorr_* columns.
  refined = list(args = list(intoverpop = "quadrature"), objective = "quadrature"),

  # The certification's bar stated in standard errors: the gap is half the
  # squared distance from the maximum in them, so 0.01 se is 5e-5 nats and
  # 0.02 se 2e-4, against the default 1e-6 (0.0014 se). L-BFGS's own rule
  # follows at a hundredth of the bar (`.ctBackendInnerGapTol()`).
  gap5e5 = list(optimcontrol = list(gaptol = 5e-5)),
  gap2e4 = list(optimcontrol = list(gaptol = 2e-4)),
  refined_gap5e5 = list(args = list(intoverpop = "quadrature"),
    optimcontrol = list(gaptol = 5e-5), objective = "quadrature"),
  refined_gap2e4 = list(args = list(intoverpop = "quadrature"),
    optimcontrol = list(gaptol = 2e-4), objective = "quadrature"),

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
