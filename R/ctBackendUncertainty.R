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
  nsubjects <- length(spec$subject_starts)
  ndatapoints <- length(spec$times)
  list(nsubjects = nsubjects, ndatapoints = ndatapoints,
    subject = rep(seq_len(max(1L, nsubjects)),
      times = diff(c(spec$subject_starts, ndatapoints + 1L))))
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

  est <- as.numeric(fit$estimate$raw)
  if (!length(est) || any(!is.finite(est))) {
    stop("The fit has no finite raw parameter estimate to work from.", call. = FALSE)
  }
  shape <- .ctBackendDataShape(fit)
  lpgFunc <- .ctBackendLpgFunc(fit)

  # `cores` here means engine threads, not R processes: the Julia engine splits
  # its own subject loop, so each log-probability evaluation is parallel and
  # there is nothing for an R cluster to do.
  if (inherits(fit, "ctJuliaFit") && cores > 1L) {
    try(JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!",
      as.integer(cores)), silent = TRUE)
  }

  # The engines produce per-subject scores from one traced pass, so they are
  # computed here and handed in rather than reconstructed a subject at a time.
  scores <- if (uncertainty %in% c("opg", "sandwich", "bootstrap")) {
    .ctBackendScoreMatrix(fit, est)
  } else NULL

  uncertaintyfit <- ctOptimComputeUncertainty(est = est, standata = shape,
    sm = NULL, lpgFunc = lpgFunc, uncertainty = uncertainty,
    finishsamples = finishsamples, cores = cores, matsetup = NA,
    control = control, verbose = verbose, scores = scores)

  if (draws == "imis") {
    if (is.null(control$imisMaxIter)) control$imisMaxIter <- 50
    if (is.null(control$imisScaleInit)) control$imisScaleInit <- 1.1
    if (is.null(control$imisTailScale)) control$imisTailScale <- 1.1
    if (is.null(control$isESS)) control$isESS <- 100
    if (is.null(control$isitersize)) control$isitersize <- 1000
    is_res <- imis_is(lpgFunc, mu_hat = est, Sigma_hat = uncertaintyfit$cov,
      max_iter = control$imisMaxIter, scale_init = control$imisScaleInit,
      tail_scale = control$imisTailScale, target_ess = control$isESS,
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
.ctBackendPriorSpec <- function(standata, npar) {
  if (is.null(standata)) {
    stop("priors=TRUE needs the prepared model data; this fit was built without it.",
      call. = FALSE)
  }
  laplace <- as.integer(standata$laplaceprior)
  if (length(laplace) && any(laplace == 1L)) {
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
