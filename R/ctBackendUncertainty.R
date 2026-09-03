# Uncertainty for optimisation-backend fits (backend='julia').
#
# `ctOptimUncertainty()` was written against `ctStanFit`, but almost none of its
# machinery is Stan-specific: `ctOptimComputeUncertainty()` needs a raw parameter
# vector, a function returning the log probability and its gradient, and the
# subject / datapoint counts. Everything else it reaches for -- the Stan model
# object, `standata` -- is used only by the score and bootstrap methods.
#
# So this file does not reimplement any of it. It builds those three things from
# a `ctJuliaFit` and calls the same `ctOptimComputeUncertainty()`,
# `ctOptimNormalDraws()` and `imis_is()` the Stan path uses. A second
# implementation of a Hessian or an importance sampler is exactly the kind of
# duplicated numerical surface this project has already been bitten by.
#
# One thing worth stating because it would otherwise need re-deriving: the Stan
# path evaluates `log_prob(..., adjust_transform = TRUE)`, and these engines
# compute the `adjust_transform = FALSE` quantity. Here that is the same number.
# ctsem declares every Stan parameter as an unconstrained `vector` (see the
# `parameters` block of inst/stan/ctsm.stan) -- the parameter transforms live
# inside the model, not in Stan's constraint syntax -- so there is no Jacobian
# adjustment to include. The Stan/Julia parity suite compares
# against `adjust_transform = FALSE` for the same reason.

# Wrap a backend's evaluate function as the `lpgFunc` contract
# `ctOptimComputeUncertainty` expects: a numeric log probability carrying its
# gradient as an attribute, and a finite fallback rather than an error at a
# point the Hessian's finite differences happen to wander into.
.ctBackendLpgFunc <- function(fit) {
  if (!inherits(fit, "ctJuliaFit")) {
    stop("Unsupported fit class for backend uncertainty.", call. = FALSE)
  }
  evaluate <- function(parm) ctJuliaEvaluate(fit, parm, gradient = TRUE)
  function(parm) {
    result <- try(evaluate(as.numeric(parm)), silent = TRUE)
    value <- if (inherits(result, "try-error")) NaN else as.numeric(result$value)[1L]
    gradient <- if (inherits(result, "try-error")) NULL else as.numeric(result$gradient)
    if (!is.finite(value) || is.null(gradient) || length(gradient) != length(parm) ||
        any(!is.finite(gradient))) {
      # Matches the Stan path's own guard: a large finite penalty with a zero
      # gradient, so a finite-difference step into an invalid region degrades
      # the local approximation rather than aborting the whole calculation.
      value <- -1e100
      gradient <- rep(0, length(parm))
    }
    attributes(value) <- list(gradient = gradient)
    value
  }
}

# The two counts `ctOptimCheckUncertaintyData` needs. The Julia spec has no
# `standata`, but it does carry the subject starts and the observation times,
# which is the same information.
.ctBackendDataShape <- function(fit) {
  spec <- fit$model_spec
  list(nsubjects = length(spec$subject_starts), ndatapoints = length(spec$times),
    # `.ctFitRowSubject()` rather than a second spelling of the same expansion.
    # This one wrote `seq_len(max(1L, nsubjects))` where that writes
    # `seq_along(starts)`; they agree on every input either can be given, but
    # two spellings of one index map is how they come to disagree.
    subject = .ctFitRowSubject(fit))
}

# `fullbootstrap` is the one method still out of reach: it resamples subjects
# and re-optimises each sample, which needs the model rebuilt per resample
# rather than just re-evaluated. Everything else these engines can serve, now
# that they produce per-subject scores directly (see .ctBackendScoreMatrix).
.ctBackendUncertaintySupported <- c("hessian", "surrogate", "is", "opg",
  "sandwich", "bootstrap")

.ctBackendUncertainty <- function(fit, uncertainty, draws, finishsamples,
  cores, control, verbose) {

  if (!uncertainty %in% .ctBackendUncertaintySupported) {
    stop("uncertainty='", uncertainty, "' is not available for backend='",
      fit$backend, "' fits. ",
      "fullbootstrap resamples subjects and re-optimises each sample, which ",
      "needs the model rebuilt per resample rather than re-evaluated. ",
      "Available: ", paste(.ctBackendUncertaintySupported, collapse = ", "),
      ".", call. = FALSE)
  }

  # A state-explicit fit has one uncertainty route and not five. The others
  # all reach for `lpgFunc`, which evaluates the *marginal* likelihood over
  # the parameters -- a different objective from the one this fit maximised,
  # and using it here would produce standard errors for a model that was not
  # fitted. `hessian` needs only the curvature, which the engine supplies
  # correctly for this target.
  if (isFALSE(fit$args$intoverstates) && !identical(uncertainty, "hessian")) {
    stop("uncertainty='", uncertainty, "' is not available for an ",
      "intoverstates=FALSE fit: it would score the marginal likelihood, ",
      "which is not what this fit maximised. Use uncertainty='hessian', ",
      "which profiles the states out, or fit with optimize=FALSE to sample ",
      "the states and read the posterior directly.", call. = FALSE)
  }

  est <- as.numeric(fit$estimate$raw)
  if (!length(est) || any(!is.finite(est))) {
    stop("The fit has no finite raw parameter estimate to work from.", call. = FALSE)
  }
  shape <- .ctBackendDataShape(fit)
  lpgFunc <- .ctBackendLpgFunc(fit)

  # `cores` here means engine threads, not R processes: the Julia engine splits
  # its own subject loop, so each log-probability evaluation is parallel and
  # there is nothing for an R cluster to do.
  #
  # The ceiling the caller asked for.
  #
  # This forced one chunk, to mitigate a penalty that does not reproduce. The
  # record it carried, swept in one direction on one machine, was:
  #
  #     chunks              1       2       4       8       1 again
  #     4 latents, 120 subj  31.2 s  68.9 s  43.0 s  30.8 s  31.5 s
  #     2 latents, 300 subj   5.0 s  11.3 s   7.3 s   5.3 s   5.4 s
  #
  # Re-measured on dev1 -- 23 homogeneous cores, no SMT -- with the sweep run
  # up and back down rather than in one direction, on the same two shapes.
  # Whole-phase seconds, minimum of three:
  #
  #     chunks               1     2     4     8     8     4     2     1
  #     4 latents, 120 subj 3.68  3.73  3.58  3.69  3.68  3.64  3.65  3.63
  #     2 latents, 300 subj 2.94  2.98  2.96  2.95  2.98  2.97  2.99  2.90
  #
  # The chunk count moves the phase by 4% and 3%, against a baseline that
  # repeats to 2.9% and 0.1%. There is no penalty here to mitigate.
  #
  # The mechanism agrees. This phase spends almost all of its wall time in
  # JuliaConnectoR round trips rather than in the subject loop: one gradient
  # costs 0.417 s called from R and 0.005 s measured inside the engine, and a
  # round trip costs 40-80 ms whatever it carries. Chunking the subject loop
  # cannot move a number that is one percent subject loop, in either
  # direction. Whatever produced the first table, it was not the chunking.
  #
  # What the chunk count does change is the exact Hessian, the one component
  # here with enough work in it to divide: 0.40 s to 0.22 s at eight chunks on
  # the first shape, 0.21 s to 0.16 s on the second. Forcing one chunk gave
  # that up and bought nothing.
  #
  # Set explicitly rather than left alone, because the session may carry a
  # ceiling from elsewhere; the engine still caps it at its own thread count.
  chunks <- suppressWarnings(as.integer(cores)[1L])
  if (is.na(chunks) || chunks < 1L) chunks <- 1L
  previous <- .ctBackendSetMaxChunks(chunks)
  # Restored on exit, like the optimiser's. This runs *after* the optimiser has
  # put its own ceiling back, so without a restore here the session still ended
  # a fit reconfigured -- measured at 4 after a `cores=4` Laplace fit, where the
  # optimiser alone had correctly left it at whatever it found.
  on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)

  # The engines produce per-subject scores from one traced pass, so they are
  # computed here and handed in rather than reconstructed a subject at a time.
  scores <- if (uncertainty %in% c("opg", "sandwich", "bootstrap")) {
    .ctBackendScoreMatrix(fit, est)
  } else NULL

  # Likewise the Hessian: exact rather than finite-differenced, because the
  # engine can differentiate its own gradient. `control$analyticHessian=FALSE`
  # falls back to the shared finite-difference path, which is worth keeping
  # reachable -- it is the only way to compare the two on a real model.
  needsHessian <- uncertainty %in% c("hessian", "sandwich", "bootstrap", "is") ||
    (uncertainty == "surrogate" && is.null(control$initialCov))
  hessian <- if (needsHessian && !identical(control$analyticHessian, FALSE)) {
    .ctBackendHessian(fit, est, verbose = verbose)
  } else NULL
  control$analyticHessian <- NULL

  uncertaintyfit <- ctOptimComputeUncertainty(est = est, standata = shape,
    sm = NULL, lpgFunc = lpgFunc, uncertainty = uncertainty,
    finishsamples = finishsamples, cores = cores, matsetup = NA,
    control = control, verbose = verbose, scores = scores, hessian = hessian)

  if (draws == "imis") {
    if (is.null(control$imisMaxIter)) control$imisMaxIter <- 50
    # Wider than the curvature says. The
    # proposal starts from the Hessian covariance, which on a modest sample is
    # *narrower* than the posterior -- the very thing importance sampling is
    # being asked to correct -- and a proposal narrower than its target cannot
    # correct it, because the region carrying the missing mass is never
    # visited. This previously defaulted to 1.1, and returned standard errors
    # within 10% of the Hessian's where the posterior was up to twice as wide.
    if (is.null(control$imisScaleInit)) control$imisScaleInit <- 1.5
    if (is.null(control$imisTailScale)) control$imisTailScale <- 1.2
    # Normal, not t. See `imis_is`: the heavier-tailed proposal was measured
    # and was worse at equal scale, and only competitive at a scale that
    # collapsed the effective sample size.
    if (is.null(control$imisDf)) control$imisDf <- Inf
    if (is.null(control$isESS)) control$isESS <- 100
    if (is.null(control$isitersize)) control$isitersize <- 1000
    is_res <- imis_is(lpgFunc, mu_hat = est, Sigma_hat = uncertaintyfit$cov,
      max_iter = control$imisMaxIter, scale_init = control$imisScaleInit,
      tail_scale = control$imisTailScale, df = control$imisDf,
      target_ess = control$isESS,
      n_batch = control$isitersize, cl = NA, finishsamples = finishsamples,
      verbose = verbose > 0)
    samples <- is_res$theta
    uncertaintyfit$proposal_cov <- uncertaintyfit$cov
    if (!is.null(is_res$covariance) && all(is.finite(is_res$covariance))) {
      uncertaintyfit$cov <- ctOptimSafeCov(is_res$covariance)
    } else if (nrow(samples) > 1) {
      uncertaintyfit$cov <- ctOptimSafeCov(stats::cov(samples))
    }
    uncertaintyfit$imis <- is_res
    uncertaintyfit$details$importance_sampling <- list(ess = is_res$ess,
      df_used = is_res$df_used,
      covariance = "weighted importance-sampling covariance")
  } else {
    samples <- ctOptimNormalDraws(est, uncertaintyfit$cov, finishsamples)
  }

  storedControl <- control
  storedControl$initialCov <- NULL
  uncertaintyfit$draws <- draws
  uncertaintyfit$settings <- list(method = uncertainty, draws = draws,
    finishsamples = finishsamples, cores = cores, control = storedControl)

  fit$estimate$cov <- uncertaintyfit$cov
  fit$estimate$se <- sqrt(diag(uncertaintyfit$cov))
  fit$estimate$rawposterior <- samples
  fit$uncertainty <- uncertaintyfit
  # New draws mean the fit's constrained draws describe the previous ones, so
  # they are refreshed here rather than left to be noticed downstream.
  fit$transformedpars <- .ctBackendConstrain(fit)
  # Deliberately no transformed-parameter summary: `ctOptimUpdateTransformed()`
  # goes through `rstan::constrain_pars` and the Stan model object, and neither
  # backend has a parameter-matrix reconstruction API yet (the same gap that
  # makes `ctSummaryMatrices()` unavailable for these fits). The raw-scale
  # covariance and draws are complete and usable; anything on the transformed
  # scale would have to be reconstructed, not merely relabelled.
  fit
}

# --- priors ----------------------------------------------------------------
#
# The generated Stan model's prior block is a sum of `normal_lpdf(x/scale|0,1)`
# terms over the raw parameter vector, and the Julia engine uses that
# same vector in that same order (which is why the parity tests can hand the
# identical `raw` to all three). So the whole of ctsem's prior semantics reduces
# to a list of (index, scale) pairs, decided here where the semantics live, and
# evaluated in the engine as a normal log-density.
#
# The layout, from `ctModelWriter.R`'s model block and `ctData.R`:
#   1..nparams                     rawpopmeans          scale 1
#   next nindvarying               rawpopsdbase         scale 1
#   next nindvaryingoffdiagonals   sqrtpcov             scale 1
#   TI effect coefficients         tipredeffectparams   scale tipredeffectscale
#
# `.ctJuliaAugmentRandomEffects` appends the population SD and correlation
# parameters in exactly that order, and `.ctJuliaTIEffects` appends the TI
# coefficients after them, so the indices line up without a separate mapping.
# Priors for the Laplace route, built from *its* raw layout rather than Stan's.
#
# `.ctBackendPriorSpec` below encodes the generated Stan model's ordering:
# population means, one block of population scales, one block of correlations,
# then TI-predictor effects. That is right for a single level and wrong for a
# hierarchy, where there is a scale block and a correlation block *per level* --
# it refuses rather than mis-assigning, which is how the gap surfaced.
#
# Every level gets the prior shape the subject level has: `normal(0,1)` on the
# untransformed scales and on the unconstrained correlation coordinates, which
# is what `rawpopsdbase` and `sqrtpcov` carry in the generated model. That is a
# deliberate choice rather than an inherited one -- the parameterisation is
# identical at every level, so the prior should be too, while noting that the
# same nominal prior does far more work over six studies than over three
# hundred subjects.
#
# `.ctBackendLaplacePriorSpec` reproduces `.ctBackendPriorSpec` exactly when
# there is one level, and `test-julia-laplace.R` asserts that rather than
# trusting it.
.ctBackendLaplacePriorSpec <- function(standata, laplace, npar) {
  if (is.null(standata)) {
    stop("priors=TRUE needs the prepared model data; this fit was built without it.",
      call. = FALSE)
  }
  .ctBackendRejectLaplacePriors(standata)

  nparams <- as.integer(standata$nparams)[1L]
  if (is.na(nparams)) nparams <- 0L
  index <- seq_len(nparams)
  scale <- rep(1, nparams)

  for (level in laplace$levels) {
    if (length(level$sd_index)) {
      index <- c(index, as.integer(level$sd_index))
      scale <- c(scale, rep(1, length(level$sd_index)))
    }
    if (length(level$cor_index)) {
      index <- c(index, as.integer(level$cor_index))
      scale <- c(scale, rep(1, length(level$cor_index)))
    }
  }

  ntipredeffects <- as.integer(standata$ntipredeffects)[1L]
  if (!is.na(ntipredeffects) && ntipredeffects > 0L) {
    tipredscale <- as.numeric(standata$tipredeffectscale)[1L]
    if (!is.finite(tipredscale) || tipredscale <= 0) tipredscale <- 1
    index <- c(index, laplace$npar + seq_len(ntipredeffects))
    scale <- c(scale, rep(tipredscale, ntipredeffects))
  }

  if (length(index) != npar) {
    stop("Cannot map ctsem's priors onto this model's raw parameters: the ",
      "Laplace layout accounts for ", length(index), " of ", npar,
      " free parameters. Please report this model shape.", call. = FALSE)
  }
  priormod <- as.numeric(standata$priormod)[1L]
  if (!is.finite(priormod)) priormod <- 1
  nsubsets <- as.numeric(standata$nsubsets)[1L]
  if (!is.finite(nsubsets) || nsubsets <= 0) nsubsets <- 1
  list(index = as.integer(index), scale = as.numeric(scale),
    weight = priormod / nsubsets)
}

# The prior families these engines cannot evaluate, refused identically
# whichever layout is in use.
.ctBackendRejectLaplacePriors <- function(standata) {
  laplaceprior <- as.integer(standata$laplaceprior)
  if (length(laplaceprior) && any(laplaceprior == 1L)) {
    stop("Laplace priors are not implemented for backend='julia'. ",
      "The generated Stan model uses a smoothed double-exponential density for ",
      "these, which these engines do not evaluate; use backend='stan', or drop ",
      "laplaceprior for the affected matrices.", call. = FALSE)
  }
  if (isTRUE(as.integer(standata$laplacetipreds)[1L] == 1L)) {
    stop("Laplace priors on TI predictor effects are not implemented for ",
      "backend='julia'; use backend='stan'.", call. = FALSE)
  }
  if (isTRUE(as.integer(standata$laplaceprioronly)[1L] == 1L)) {
    stop("laplaceprioronly is not implemented for backend='julia'; ",
      "use backend='stan'.", call. = FALSE)
  }
  invisible(TRUE)
}

.ctBackendPriorSpec <- function(standata, npar) {
  if (is.null(standata)) {
    stop("priors=TRUE needs the prepared model data; this fit was built without it.",
      call. = FALSE)
  }
  .ctBackendRejectLaplacePriors(standata)

  nparams <- as.integer(standata$nparams)[1L]
  nindvarying <- as.integer(standata$nindvarying)[1L]
  noffdiagonals <- as.integer(standata$nindvaryingoffdiagonals)[1L]
  if (is.na(nparams)) nparams <- 0L
  if (is.na(nindvarying)) nindvarying <- 0L
  if (is.na(noffdiagonals)) noffdiagonals <- 0L

  index <- seq_len(nparams)
  scale <- rep(1, nparams)
  position <- nparams
  if (nindvarying > 0L) {
    index <- c(index, position + seq_len(nindvarying))
    scale <- c(scale, rep(1, nindvarying))
    position <- position + nindvarying
    if (nindvarying > 1L && noffdiagonals > 0L) {
      index <- c(index, position + seq_len(noffdiagonals))
      scale <- c(scale, rep(1, noffdiagonals))
      position <- position + noffdiagonals
    }
  }
  ntipredeffects <- as.integer(standata$ntipredeffects)[1L]
  if (!is.na(ntipredeffects) && ntipredeffects > 0L) {
    tipredscale <- as.numeric(standata$tipredeffectscale)[1L]
    if (!is.finite(tipredscale) || tipredscale <= 0) tipredscale <- 1
    index <- c(index, position + seq_len(ntipredeffects))
    scale <- c(scale, rep(tipredscale, ntipredeffects))
    position <- position + ntipredeffects
  }

  # If the engine's free-parameter count and the Stan-side layout disagree, the
  # indices are meaningless and a silently mis-scaled posterior is far worse
  # than a refusal.
  if (position != npar) {
    stop("Cannot map ctsem's priors onto this model's raw parameters: the Stan ",
      "layout accounts for ", position, " of ", npar, " free parameters. ",
      "Please report this model shape.", call. = FALSE)
  }
  priormod <- as.numeric(standata$priormod)[1L]
  if (!is.finite(priormod)) priormod <- 1
  nsubsets <- as.numeric(standata$nsubsets)[1L]
  if (!is.finite(nsubsets) || nsubsets <= 0) nsubsets <- 1
  list(index = as.integer(index), scale = as.numeric(scale),
    weight = priormod / nsubsets)
}

# --- the Hessian ------------------------------------------------------------

# The exact Hessian of the log posterior, from the engine.
#
# The shared path finite-differences the gradient: `2 * npar` reverse sweeps,
# accurate to about the square root of machine precision, and only if the step
# suits the parameter's scale -- which one global step cannot do for a vector
# mixing log standard deviations with unconstrained correlations. The engine
# instead differentiates its own reverse-mode gradient in forward mode, which
# costs `ceil(npar / chunksize)` sweeps and is exact. Measured against
# `ForwardDiff.hessian` of the primal on a state-dependent 13-parameter model,
# the exact route agrees to 2e-16 relative and the finite difference to 2e-6.
#
# Returns NULL rather than erroring if anything goes wrong, so the caller falls
# back to the finite-difference Hessian: a worse covariance is a better outcome
# than no fit at all, and the two are interchangeable at this seam.
.ctBackendHessian <- function(fit, est, verbose = 0, gradient = "adjoint") {
  # A state-explicit fit maximised a different object, so its curvature is a
  # different object too. See `.ctBackendJointHessian`.
  if (isFALSE(fit$args$intoverstates)) return(.ctBackendJointHessian(fit, est))
  module <- .ctJuliaModule(.ctBackendSpec(fit)$project)
  # `gradient='forward'` selects `ctsem_hessian_forward` instead of
  # `ctsem_hessian`: a model with a sampled (missing) TI predictor value
  # forces exactly this (see `.ctFitJuliaBackendImpl`), because
  # `ctsem_hessian` nests over the adjoint gradient, which refuses such a
  # model outright rather than silently omitting its cotangent.
  hessian_fn_name <- if (identical(gradient, "forward")) "ctsem_hessian_forward" else "ctsem_hessian"
  # A user whose cached engine environment predates this function has no such
  # function, and that is a silent fallback rather than an error: the engine
  # environment is keyed on a hash of the engine's source, so it refreshes
  # itself the moment that source changes.
  available <- isTRUE(tryCatch(is.function(module[[hessian_fn_name]]),
    error = function(e) FALSE))
  if (!available) return(NULL)
  # Said out loud, as the finite-difference path says "Estimating Hessian",
  # because the first call on a given model shape spends about ten seconds in
  # Julia compiling the reverse pass for dual numbers. That is once per model
  # per session -- the second call on this 28-parameter model takes 0.08 s,
  # against 0.78 s for the 2*npar gradient evaluations a finite difference
  # needs -- but an unexplained ten-second pause is worth a line of output.
  message("Computing exact Hessian")
  result <- try(.ctBackendJuliaValue(module[[hessian_fn_name]](
    .ctJuliaObjective(fit), .ctJuliaVector(as.numeric(est)))), silent = TRUE)
  if (inherits(result, "try-error")) {
    warning("The engine could not differentiate its gradient here; ",
      "falling back to the finite-difference Hessian.", call. = FALSE)
    return(NULL)
  }
  hessian <- matrix(as.numeric(result), nrow = length(est), ncol = length(est))
  if (any(!is.finite(hessian))) {
    warning("The exact Hessian was not finite at the estimate; ",
      "falling back to the finite-difference Hessian.", call. = FALSE)
    return(NULL)
  }
  hessian
}

# --- per-subject scores -----------------------------------------------------

# The score matrix `ctOptimScoreMatrix()`/`bootstrapHessian()` consume, from the
# engines' own per-subject adjoint contributions rather than by re-initialising
# a model per subject the way `scorecalc()` must for Stan.
.ctBackendScoreMatrix <- function(fit, est) {
  module <- .ctJuliaModule(fit$model_spec$project)
  result <- JuliaConnectoR::juliaGet(module$ctsem_subject_gradients(
    .ctJuliaObjective(fit), .ctJuliaVector(as.numeric(est))))
  scores <- as.matrix(result$scores)
  if (any(!is.finite(scores))) {
    stop("The engine returned non-finite per-subject scores at the estimate; ",
      "score-based uncertainty cannot be computed here.", call. = FALSE)
  }
  scores
}

# The Hessian of an `intoverstates=FALSE` fit.
#
# Not the second derivative of the joint density in the parameter block: that
# is the curvature at *fixed* states, which treats a trajectory that was
# estimated as though it had been observed and gives intervals far too
# narrow. What a standard error needs here is the curvature of the profile,
# which is the Schur complement of the joint Hessian -- equivalently the
# Laplace approximation to the observed information with the states
# integrated back out. The engine forms it per subject, so no matrix of the
# full state dimension is ever built; see `ctsem_joint_hessian`.
.ctBackendJointHessian <- function(fit, est) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  innovations <- as.numeric(fit$estimate$innovations)
  if (!length(innovations)) return(NULL)
  joint <- try(.ctJuliaJointObjective(fit, length(est)), silent = TRUE)
  if (inherits(joint, "try-error")) return(NULL)
  message("Computing exact Hessian, with the states profiled out")
  result <- try(.ctBackendJuliaValue(module$ctsem_joint_hessian(joint,
    .ctJuliaNumericVector(c(as.numeric(est), innovations)),
    profile = TRUE)), silent = TRUE)
  if (inherits(result, "try-error") || is.null(result)) {
    warning("The profiled Hessian could not be formed at the estimate -- ",
      "usually a subject whose trajectory the data does not determine. ",
      "Falling back to the finite-difference Hessian, which differentiates ",
      "the marginal likelihood and so does not describe this fit.",
      call. = FALSE)
    return(NULL)
  }
  hessian <- matrix(as.numeric(result), nrow = length(est), ncol = length(est))
  if (any(!is.finite(hessian))) return(NULL)
  hessian
}
