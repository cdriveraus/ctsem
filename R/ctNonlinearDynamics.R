# Model implied dynamics for a context-dependent model -------------------------
#
# `ctDiscretePars()` reports expm(DRIFT * t): the response, over an interval, to
# a one unit impulse on each process. For a linear model that is exact. For a
# model whose DRIFT depends on the state it is the response of the model
# *linearised at one point*, which is a different and weaker claim -- the real
# system's DRIFT changes as the trajectory moves, and over a long interval it
# may change a lot.
#
# Two honest alternatives, both of which need the nonlinear system integrated
# rather than exponentiated:
#
#   'simulate'   Integrate from a starting state, and again from that state plus
#                one unit on process j, and difference the two trajectories.
#                This is the actual model implied regression: it starts at 1 by
#                construction and, for a stable process, decays to 0 -- so it
#                reads exactly like the linear panel it replaces, and reduces to
#                it when no cell depends on the state.
#
#   'empirical'  The interval Jacobians the filter already computed along each
#                subject's realised trajectory. Not implemented here yet; the
#                engine returns them from ctsem_kalman as `transition` and the R
#                side does not surface them.
#
# The integrator mirrors the engine's own: over each substep the system is
# treated as affine, dx/dt = A x + b with A the analytic Jacobian JAx at the
# current state and b whatever makes that consistent with the true derivative
# there. That is what the extended Kalman filter does between observations, so
# a trajectory from here matches one the filter would have produced, and it
# costs one engine call per substep rather than the four an RK4 would.

# Derivative and local Jacobian of the deterministic system at one state.
.ctNonlinearDerivative <- function(fit, state, tipreds = NULL, nlatent, raw = NULL) {
  mats <- suppressMessages(ctBackendParMatrices(fit, raw = raw, state = state,
    tipreds = tipreds, trim = FALSE))
  index <- seq_len(nlatent)
  drift <- mats$DRIFT[index, index, drop = FALSE]
  cint <- as.numeric(mats$CINT[index, 1])
  jacobian <- mats$JAx[index, index, drop = FALSE]
  # A model with no state dependence has no analytic JAx to speak of; the drift
  # matrix is then exactly the Jacobian.
  if (!length(jacobian) || all(!is.finite(jacobian)) || all(jacobian == 0)) {
    jacobian <- drift
  }
  list(derivative = as.numeric(drift %*% state[index]) + cint, jacobian = jacobian)
}

# One affine step of length dt, exact for the affine system it linearises.
.ctNonlinearStep <- function(state, derivative, jacobian, dt, nlatent) {
  index <- seq_len(nlatent)
  transition <- as.matrix(Matrix::expm(jacobian * dt))
  # b is the affine offset consistent with the true derivative at this state.
  offset <- derivative - as.numeric(jacobian %*% state[index])
  shift <- try(solve(jacobian, (transition - diag(nlatent)) %*% offset), silent = TRUE)
  # A singular Jacobian means no closed form for the offset's contribution;
  # over one substep the first order term is a good enough stand-in.
  if (inherits(shift, 'try-error')) shift <- offset * dt
  state[index] <- as.numeric(transition %*% state[index]) + as.numeric(shift)
  state
}

#' Deterministic trajectory of the (possibly nonlinear) system
#'
#' @param fit A ctJuliaFit.
#' @param state Starting state, at engine (augmented) length.
#' @param times Times at which to report the state. Must start at 0.
#' @param tipreds Time independent predictor values, or NULL for population.
#' @param maxstep Largest integration substep.
#' @return length(times) by nlatent matrix.
#' @noRd
.ctNonlinearTrajectory <- function(fit, state, times, tipreds = NULL, maxstep = 0.1,
  raw = NULL) {
  nlatent <- .ctFitNlatent(fit)
  times <- sort(unique(c(0, as.numeric(times))))
  out <- matrix(NA_real_, length(times), nlatent)
  out[1L, ] <- state[seq_len(nlatent)]
  current <- state
  for (k in seq_along(times)[-1L]) {
    span <- times[k] - times[k - 1L]
    nsteps <- max(1L, ceiling(span / maxstep))
    dt <- span / nsteps
    for (step in seq_len(nsteps)) {
      local <- .ctNonlinearDerivative(fit, current, tipreds, nlatent, raw)
      current <- .ctNonlinearStep(current, local$derivative, local$jacobian, dt, nlatent)
      if (!all(is.finite(current[seq_len(nlatent)]))) return(out)
    }
    out[k, ] <- current[seq_len(nlatent)]
  }
  out
}

#' Model implied regression by simulation, for a nonlinear system
#'
#' The response to a one unit impulse on each process, as a difference from the
#' unperturbed trajectory: integrate from \code{state}, integrate again from
#' \code{state} with one unit added to process j, and subtract. Column j of the
#' result at time t is what a one unit change in process j is worth, t later.
#'
#' Taking the difference rather than the perturbed trajectory itself is what
#' makes this readable: the unperturbed run carries the system's own drift
#' towards its attractor, which is not an effect of the impulse. Differenced,
#' the response starts at exactly 1 on the diagonal and decays to 0 for a stable
#' process, the same shape the linear panel has.
#'
#' @return length(times) by nlatent by nlatent array, [time, response, impulse].
#' @noRd
.ctNonlinearImpulseResponse <- function(fit, state, times, tipreds = NULL,
  maxstep = 0.1, raw = NULL) {

  nlatent <- .ctFitNlatent(fit)
  times <- sort(unique(c(0, as.numeric(times))))
  baseline <- .ctNonlinearTrajectory(fit, state, times, tipreds, maxstep, raw)
  out <- array(NA_real_, dim = c(length(times), nlatent, nlatent))
  for (impulse in seq_len(nlatent)) {
    perturbed <- state
    perturbed[impulse] <- perturbed[impulse] + 1
    out[, , impulse] <- .ctNonlinearTrajectory(fit, perturbed, times, tipreds,
      maxstep, raw) - baseline
  }
  dimnames(out) <- list(NULL, .ctFitModelObject(fit)$latentNames,
    .ctFitModelObject(fit)$latentNames)
  out
}

# ctDiscretePars(method='simulate') ------------------------------------------
#
# Same output shape as the linearised path -- [sample, subject, time, row, col]
# -- so ctDiscreteParsPlot and everything downstream is unchanged.
#
# One trajectory pair per posterior draw per process, each of many substeps,
# each substep an engine call. That is thousands of calls where the linearised
# path does one matrix exponential, so the default draw count is small and
# stated rather than silently inherited from a plotting argument sized for a
# cheap computation.
.ctDiscreteParsSimulate <- function(fit, times, state = 'asymptotic',
  nsamples = 10, maxstep = 0.1, quiet = FALSE) {

  if (!inherits(fit, 'ctJuliaFit')) stop(call. = FALSE,
    paste0("method='simulate' requires backend='julia': integrating the ",
      "nonlinear system needs the matrices re-materialised at each step, which ",
      "the stan backend cannot do after fitting."))

  resolved <- .ctResolveState(fit, state)
  start <- if (is.null(resolved$state)) .ctContextBaseState(fit) else resolved$state

  draws <- .ctBackendRawSamples(fit)
  nsamples <- max(1L, min(as.integer(nsamples), nrow(draws)))
  if (nsamples < nrow(draws)) {
    draws <- draws[round(seq(1, nrow(draws), length.out = nsamples)), , drop = FALSE]
  }
  if (!quiet) message('Simulating the nonlinear impulse response from ',
    resolved$label, ', over ', nsamples, ' draw', if (nsamples > 1) 's' else '',
    '. This integrates the system rather than exponentiating a frozen DRIFT, ',
    'so it is much slower than method="linearise".')

  times <- sort(unique(c(0, as.numeric(times))))
  nlatent <- .ctFitNlatent(fit)
  out <- array(NA_real_, dim = c(nsamples, 1L, length(times), nlatent, nlatent))
  for (draw in seq_len(nsamples)) {
    out[draw, 1L, , , ] <- .ctNonlinearImpulseResponse(fit, start, times,
      maxstep = maxstep, raw = draws[draw, ])
  }
  attr(out, 'times') <- times
  attr(out, 'stateLabel') <- resolved$label
  out
}
