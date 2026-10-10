# Bridges to the engine's per-row scores and conditional effect likelihoods
# (inst/julia/ContinuousTimeSEM/src/row_scores.jl), for the specification
# tests in ctScoreTest() and ctFitAdequacy().

# Each data row's score along each column of `directions` (raw parameters by
# K): the derivative of the row's log one-step-ahead predictive density, rows
# in the order ctKalmanArray reports them. On the Laplace route the rows are
# the sequential predictions' (each row given the subject's earlier rows), so
# the scores are a martingale difference sequence within a subject under the
# model, as on the filter routes. A subject's rows sum to its score.
.ctBackendRowScores <- function(fit, est, directions, cores = 2L) {
  directions <- as.matrix(directions)
  storage.mode(directions) <- 'double'
  if (nrow(directions) != length(est)) stop('directions needs one row per raw parameter.',
    call. = FALSE)
  module <- .ctJuliaModule(fit$model_spec$project)
  out <- .ctBackendWithMaxChunks(cores, .ctBackendJuliaValue(module$ctsem_row_scores(
    .ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(as.numeric(est)),
    .ctJuliaPut(directions))))
  matrix(as.numeric(out), ncol = ncol(directions),
    dimnames = list(NULL, colnames(directions)))
}

# Each unit's posterior over its standardised random effects, as the adaptive
# rule the quadrature integrates with (`nodes` a dimension): `points` (units by
# nodes^d by d), normalised `weights` (units by nodes^d), and `logmarginal`,
# each unit's quadrature term. Single-level Laplace fits.
.ctBackendPosteriorNodes <- function(fit, est, nodes = 5L, cores = 2L) {
  module <- .ctJuliaModule(fit$model_spec$project)
  out <- .ctBackendWithMaxChunks(cores, .ctBackendJuliaValue(
    module$ctsem_laplace_posterior_nodes(.ctJuliaObjective(fit, hybrid = TRUE),
      .ctJuliaNumericVector(as.numeric(est)), nodes = as.integer(nodes))))
  W <- as.matrix(out$weights)
  d <- length(out$points) / max(1L, length(W))
  list(points = array(as.numeric(out$points), c(nrow(W), ncol(W), d)), weights = W,
    logmarginal = as.numeric(out$logmarginal))
}

# Each unit's number of random-effect coordinates on the Laplace route.
.ctBackendUnitDims <- function(fit) {
  module <- .ctJuliaModule(fit$model_spec$project)
  as.integer(.ctBackendJuliaValue(module$ctsem_laplace_unit_dims(
    .ctJuliaObjective(fit, hybrid = TRUE))))
}
