# Uncertainty for optimisation-backend fits (backend='julia', backend='cpp').
#
# `ctOptimUncertainty()` was written against `ctStanFit`, but almost none of its
# machinery is Stan-specific: `ctOptimComputeUncertainty()` needs a raw parameter
# vector, a function returning the log probability and its gradient, and the
# subject / datapoint counts. Everything else it reaches for -- the Stan model
# object, `standata` -- is used only by the score and bootstrap methods.
#
# So this file does not reimplement any of it. It builds those three things from
# a `ctJuliaFit` or `ctCppFit` and calls the same `ctOptimComputeUncertainty()`,
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
# adjustment to include. The Stan/Julia and Stan/C++ parity suites compare
# against `adjust_transform = FALSE` for the same reason.

# Wrap a backend's evaluate function as the `lpgFunc` contract
# `ctOptimComputeUncertainty` expects: a numeric log probability carrying its
# gradient as an attribute, and a finite fallback rather than an error at a
# point the Hessian's finite differences happen to wander into.
.ctBackendLpgFunc <- function(fit) {
  evaluate <- if (inherits(fit, "ctJuliaFit")) {
    function(parm) ctJuliaEvaluate(fit, parm, gradient = TRUE)
  } else if (inherits(fit, "ctCppFit")) {
    function(parm) ctCppEvaluate(fit, parm, gradient = TRUE)
  } else {
    stop("Unsupported fit class for backend uncertainty.", call. = FALSE)
  }
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

# The two counts `ctOptimCheckUncertaintyData` needs. The Julia/C++ spec has no
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

.ctBackendUncertaintySupported <- c("hessian", "surrogate", "is")

.ctBackendUncertainty <- function(fit, uncertainty, draws, finishsamples,
  cores, control, verbose) {

  if (!uncertainty %in% .ctBackendUncertaintySupported) {
    stop("uncertainty='", uncertainty, "' is not available for backend='",
      fit$backend, "' fits. ",
      "The score-based methods (opg, sandwich, bootstrap) need per-subject ",
      "score contributions and fullbootstrap needs subject resampling with ",
      "refits, neither of which these engines expose yet. Available: ",
      paste(.ctBackendUncertaintySupported, collapse = ", "), ".", call. = FALSE)
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

  uncertaintyfit <- ctOptimComputeUncertainty(est = est, standata = shape,
    sm = NULL, lpgFunc = lpgFunc, uncertainty = uncertainty,
    finishsamples = finishsamples, cores = cores, matsetup = NA,
    control = control, verbose = verbose)

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
