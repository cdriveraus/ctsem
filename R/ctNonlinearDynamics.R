# Model implied dynamics for a context-dependent model -------------------------
#
# `ctDiscretePars()` reports expm(J * t), J the drift's Jacobian at one point
# (JAx; DRIFT itself for a linear model): the response, over an interval, to a
# one unit impulse on each process. For a linear model that is exact. For a
# model whose dynamics depend on the state it is the response of the model
# *linearised at one point*, which is a different and weaker claim -- the real
# system's Jacobian changes as the trajectory moves, and over a long interval it
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
# The trajectories are the engine's (`ctsem_state_trajectories`): the filter's
# own mean step at the state, over substeps of at most `maxstep`, every draw and
# every trajectory in one call. Integrated here with an engine call per substep
# of each trajectory, 27 lags of the ctsem book's nonlinear model at 20 draws
# took over seven minutes.

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

#' Deterministic trajectories of the (possibly nonlinear) system
#'
#' @param fit A ctJuliaFit.
#' @param raw Raw parameters, npar by K: one column per trajectory.
#' @param starts Starting states, nlatent by K.
#' @param times Times at which to report the states, sorted, the first 0.
#' @param tipreds Time independent predictor values, or NULL for population.
#' @param maxstep Largest integration substep.
#' @return length(times) by nlatent by K array; NaN from where a trajectory
#'   leaves the finite reals.
#' @noRd
.ctNonlinearTrajectories <- function(fit, raw, starts, times, tipreds = NULL,
  maxstep = 0.1) {
  raw <- as.matrix(raw)
  starts <- as.matrix(starts)
  storage.mode(raw) <- storage.mode(starts) <- "double"
  module <- .ctJuliaModule(.ctBackendSpec(fit)$project)
  arguments <- list(.ctJuliaObjective(fit), .ctJuliaPut(raw), .ctJuliaPut(starts),
    .ctJuliaVector(as.numeric(times)), maxstep = as.numeric(maxstep)[1L])
  # Passed only when non-empty: JuliaConnectoR hangs on a zero-length vector.
  if (length(tipreds)) arguments$tipreds <- .ctJuliaVector(as.numeric(tipreds))
  out <- .ctBackendJuliaValue(do.call(module$ctsem_state_trajectories, arguments))
  array(as.numeric(out), c(length(times), nrow(starts), ncol(raw)))
}

#' Model implied regressions by simulation, for many draws at once
#'
#' For each draw, the response to a shock on each process as a difference from
#' that draw's unperturbed trajectory (see \code{.ctNonlinearImpulseResponse}),
#' every trajectory of every draw integrated in one engine call.
#'
#' @param raw Raw parameters, npar by D.
#' @param shocks nlatent by nlatent by D: column j of slice d the perturbation
#'   applied for impulse j under draw d.
#' @return D by length(times) by nlatent by nlatent array,
#'   [draw, time, response, impulse].
#' @noRd
.ctNonlinearImpulseResponses <- function(fit, state, times, raw, shocks,
  tipreds = NULL, maxstep = 0.1) {
  nlatent <- .ctFitNlatent(fit)
  times <- sort(unique(c(0, as.numeric(times))))
  raw <- as.matrix(raw)
  ndraws <- ncol(raw)
  start <- as.numeric(state[seq_len(nlatent)])
  # Per draw, its baseline then one trajectory per impulse.
  width <- nlatent + 1L
  starts <- matrix(start, nlatent, ndraws * width)
  for (d in seq_len(ndraws)) {
    starts[, (d - 1L) * width + 1L + seq_len(nlatent)] <- start + shocks[, , d]
  }
  paths <- .ctNonlinearTrajectories(fit, raw[, rep(seq_len(ndraws), each = width),
    drop = FALSE], starts, times, tipreds, maxstep)
  out <- array(NA_real_, c(ndraws, length(times), nlatent, nlatent))
  for (d in seq_len(ndraws)) {
    baseline <- paths[, , (d - 1L) * width + 1L]
    for (impulse in seq_len(nlatent)) {
      out[d, , , impulse] <- paths[, , (d - 1L) * width + 1L + impulse] - baseline
    }
  }
  out
}

#' Model implied regression by simulation, for a nonlinear system
#'
#' The response to a shock on each process, as a difference from the unperturbed
#' trajectory: integrate from \code{state}, integrate again from \code{state}
#' plus \code{shock[, j]}, and subtract.
#'
#' Taking the difference rather than the perturbed trajectory itself is what
#' makes this readable: the unperturbed run carries the system's own drift
#' towards its attractor, which is not an effect of the shock. Differenced, the
#' response decays to 0 for a stable process, the same shape the linear panel
#' has.
#'
#' @param shock nlatent by nlatent matrix, column j the perturbation applied for
#'   impulse j, in the processes' own units. Defaults to the identity, one unit
#'   on process j alone. **The magnitude matters**: a nonlinear system's
#'   response to a shock of size 2 is not twice its response to a shock of size
#'   1, so the caller has to choose a magnitude that means something rather than
#'   inherit an arbitrary 1.
#' @return length(times) by nlatent by nlatent array, [time, response, impulse].
#' @noRd
.ctNonlinearImpulseResponse <- function(fit, state, times, tipreds = NULL,
  maxstep = 0.1, raw = NULL, shock = NULL) {

  nlatent <- .ctFitNlatent(fit)
  times <- sort(unique(c(0, as.numeric(times))))
  if (is.null(shock)) shock <- diag(nlatent)
  if (is.null(raw)) raw <- fit$estimate$raw
  out <- .ctNonlinearImpulseResponses(fit, state, times, raw = matrix(raw),
    shocks = array(shock, c(nlatent, nlatent, 1L)), tipreds = tipreds,
    maxstep = maxstep)
  out <- array(out[1L, , , ], c(length(times), nlatent, nlatent))
  latentNames <- .ctFitModelObject(fit)$latentNames[seq_len(nlatent)]
  dimnames(out) <- list(NULL, latentNames, latentNames)
  out
}

# The shock applied for each process, and what to divide the response by.
#
# A nonlinear response does not scale with the shock, so "one unit" is not a
# neutral default -- it is a choice, and a bad one for a process whose natural
# spread is 0.02 or 200. The shock is therefore one *standard deviation* of the
# process, from the stationary within-subject covariance at the evaluation
# point: a magnitude the model itself supplies, in the processes' own units.
#
# What comes along with that shock is the companion matrix, and which one is
# the caller's choice -- see R/ctCompanionShock.R. Note that the magnitude and
# the pattern come from different places on purpose: the magnitude from the
# stationary covariance, because that is what puts a shock on the scale of the
# process, and the pattern from whichever covariance the chosen interpretation
# is about.
#
# The response is divided so the panel reads as the linear one does: by sd_r
# when standardising, by sd_c otherwise. Both reduce exactly to the linearised
# answer for a linear model, whichever companion matrix is in use -- the
# scalings cancel in the same way.
.ctNonlinearShockSpec <- function(mats, nlatent, impulseType, standardise,
  magnitude = 1) {

  index <- seq_len(nlatent)
  variance <- diag(as.matrix(mats$asymDIFFUSIONcov[index, index, drop = FALSE]))
  if (any(!is.finite(variance)) || any(variance < 0)) return(NULL)
  sdv <- sqrt(variance + 1e-10)

  companion <- .ctCompanionMatrix(.ctCompanionType(impulseType),
    mats$DIFFUSIONcov, mats$asymDIFFUSIONcov, nlatent)
  if (is.null(companion)) return(NULL)
  # Column c: a shock of one sd_c in process c, times the companion ratios, so
  # shock[c,c] = sd_c and shock[r,c] = sd_c * C[r,c]. The magnitude comes from
  # the stationary covariance because that is what puts a shock on the scale of
  # the process; the pattern comes from whichever covariance the chosen
  # interpretation asks about. Different questions, different matrices.
  shock <- companion %*% diag(sdv, nlatent)

  divisor <- if (isTRUE(standardise)) sdv %o% rep(1, nlatent) else rep(1, nlatent) %o% sdv
  list(shock = shock * magnitude, divisor = divisor * magnitude, sd = sdv)
}

# Each draw's impulse responses from `state`, each divided as its shock says
# (`.ctNonlinearShockSpec()`): draws x length(times) x nlatent x nlatent, NA for
# a draw with no stationary variance at the state, which attribute
# 'nonstationary' counts. Two engine calls whatever the number of draws: one
# for every draw's matrices at the state, one for every trajectory.
.ctNonlinearDrawResponses <- function(fit, draws, state, times, tipreds = NULL,
  impulseType = 'unit', standardise = FALSE, magnitude = 1, maxstep = 0.1) {
  nlatent <- .ctFitNlatent(fit)
  ndraws <- nrow(draws)
  out <- array(NA_real_, c(ndraws, length(times), nlatent, nlatent))
  # Per draw: the shock's own scale is a function of the parameters, so it
  # moves with them rather than being fixed from the point estimate.
  mats <- .ctBackendPopArrays(fit, samples = draws, tipreds = tipreds,
    state = state)
  slice <- function(a, d) matrix(a[d, , ], dim(a)[2L], dim(a)[3L])
  specs <- lapply(seq_len(ndraws), function(d) .ctNonlinearShockSpec(
    list(asymDIFFUSIONcov = slice(mats$pop_asymDIFFUSIONcov, d),
      DIFFUSIONcov = slice(mats$pop_DIFFUSIONcov, d)),
    nlatent, impulseType, standardise, magnitude))
  ok <- which(!vapply(specs, is.null, logical(1)))
  if (length(ok)) {
    shocks <- vapply(specs[ok], function(s) as.matrix(s$shock),
      matrix(0, nlatent, nlatent))
    response <- .ctNonlinearImpulseResponses(fit, state, times,
      raw = t(draws[ok, , drop = FALSE]),
      shocks = array(shocks, c(nlatent, nlatent, length(ok))),
      tipreds = tipreds, maxstep = maxstep)
    for (j in seq_along(ok)) for (k in seq_along(times)) {
      out[ok[j], k, , ] <- response[j, k, , ] / specs[[ok[j]]]$divisor
    }
  }
  attr(out, 'nonstationary') <- ndraws - length(ok)
  out
}

# ctDiscretePars(method='simulate') ------------------------------------------
#
# Same output shape as the linearised path -- [sample, subject, time, row, col]
# -- so ctDiscreteParsPlot and everything downstream is unchanged.
#
# One trajectory pair per posterior draw per process, each of many substeps.
# Far more work than the linearised path's one matrix exponential a lag, so the
# default draw count is small and stated rather than silently inherited from a
# plotting argument sized for a cheap computation.
.ctDiscreteParsSimulate <- function(fit, times, state = 'asymptotic',
  nsamples = 10, maxstep = 0.1, quiet = FALSE, impulseType = 'unit',
  standardise = FALSE, magnitude = 1) {

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
    '. This integrates the system rather than exponentiating its Jacobian, ',
    'so it is slower than method="linearise".')

  times <- sort(unique(c(0, as.numeric(times))))
  nlatent <- .ctFitNlatent(fit)
  responses <- .ctNonlinearDrawResponses(fit, draws, start, times,
    impulseType = impulseType, standardise = standardise, magnitude = magnitude,
    maxstep = maxstep)
  # [draw, time, row, col] to [draw, subject, time, row, col]: a unit dimension
  # in second place leaves the order of the values as it is.
  out <- array(responses, c(nsamples, 1L, length(times), nlatent, nlatent))
  nonstationary <- attr(responses, 'nonstationary')
  if (nonstationary > 0 && !quiet) message(nonstationary, ' of ', nsamples,
    ' draws had no stationary variance at this state, so no shock magnitude ',
    'could be derived from the model; those draws are NA.')

  attr(out, 'times') <- times
  attr(out, 'stateLabel') <- resolved$label
  out
}

# ctPredictTIP's dynamics panel, for a context-dependent model ----------------
#
# The linear panel is built by fabricating a dataset with every manifest set to
# missing, one pseudo-subject per covariate level, refitting the data and
# reading each pseudo-subject's matrices. For a linear model that is the right
# design: it is the simplest and most interpretable way to isolate the
# covariate's effect, and the dynamics do not depend on where the trajectory
# went, so the fabrication costs nothing.
#
# For a model whose matrices depend on the state it costs everything: what gets
# read is the DRIFT at the last row of a dataset with no observations at all,
# i.e. wherever that pseudo-subject's unconstrained prior trajectory happened to
# drift to over the full time range. The covariate effect is then confounded
# with a fabrication artefact.
#
# The analogue that does work: for each covariate level, evaluate at the state
# that level implies -- its own asymptote -- and take the impulse response as a
# difference from it. So a one unit impulse still starts at 1 and decays to zero
# for a stable process, reading exactly as the linear panel does, while the
# level-to-level differences are the covariate's real effect on the dynamics
# rather than an effect of where each fabricated trajectory ended up.
#
# It also does away with the fabricated data: covariate values go to the engine
# directly as `tipreds`, so nothing depends on a pseudo-subject at all.
.ctPredictTIPDynamics <- function(fit, tipredIndex, values, times, ntipred,
  nsamples = 5, latentNames, quiet = FALSE, impulseType = 'unit',
  standardise = FALSE) {

  draws <- .ctBackendRawSamples(fit)
  nsamples <- max(1L, min(as.integer(nsamples), nrow(draws)))
  if (nsamples < nrow(draws)) {
    draws <- draws[round(seq(1, nrow(draws), length.out = nsamples)), , drop = FALSE]
  }
  times <- sort(unique(c(0, as.numeric(times))))
  # A discrete time model steps whole time units, as ctDiscretePars() reads it.
  if (!isTRUE(.ctFitModelObject(fit)$continuoustime)) times <- unique(round(times))
  nlatent <- length(latentNames)
  out <- array(NA_real_, dim = c(nsamples, length(values), length(times),
    nlatent, nlatent))

  if (!quiet) message('Model matrices depend on the latent state, so the ',
    'dynamics for each covariate level are simulated from that level\'s own ',
    'asymptotic state rather than read off a frozen DRIFT.')

  for (level in seq_along(values)) {
    tipreds <- rep(0, ntipred)
    tipreds[tipredIndex] <- values[level]
    state <- try(.ctContextAsymptoticState(fit, tipreds = tipreds), silent = TRUE)
    if (inherits(state, 'try-error')) {
      # No fixed point at this covariate level: report nothing for it rather
      # than silently substituting a different level's state.
      warning(call. = FALSE, 'No asymptotic state at covariate value ',
        signif(values[level], 3), '; its dynamics are not shown.')
      next
    }
    out[, level, , , ] <- .ctNonlinearDrawResponses(fit, draws, state, times,
      tipreds = tipreds, impulseType = impulseType, standardise = standardise)
  }

  dimnames(out) <- list(Sample = seq_len(nsamples), Subject = seq_along(values),
    `Time interval` = times, row = latentNames, col = latentNames)
  attributes(out)$impulseType <- .ctCompanionType(impulseType)
  attributes(out)$cov <- FALSE
  attributes(out)$standardise <- standardise
  attributes(out)$method <- 'simulate'
  out
}
