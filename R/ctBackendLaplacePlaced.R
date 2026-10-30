# The refined fit with its flagged units' rules placed at every evaluation
# (`laplace_placed.jl`), in place of the rounds (`.ctLaplaceContinueRun()`):
# `ctsem_laplace_placed_optimize` climbs the placed objective from the Laplace
# estimate, L-BFGS and then the Newton finish on a Hessian of differenced exact
# gradients, and the uncertainty is drawn on that same objective, whose gradient
# is its value's. Units on frozen draws are placed again only where their
# weights' effective size has collapsed, the run's L-BFGS memory dropped. Option
# B of review/SURVEY-quadrature-2026-10-30.md; `.ctLaplaceContinueDefaults$
# placed = FALSE` sends a fit to the rounds instead.
#
# Called by `.ctLaplaceContinue()` once the continuation `cont` is placed at the
# estimate `est` and `record` describes that placement; `report` writes the
# quadrature log likelihood into the fit and `timed` the wall seconds by stage.
# `stay` moves nothing: the rules placed at the estimate, its Hessian and the
# uncertainty drawn there. Returns the fit.
.ctLaplacePlacedRefine <- function(fit, module, cont, est, record, control, tol,
  stay, chunks, verbose, sink, method, uncertaintycontrol, report, timed) {
  get <- .ctJuliaGet
  num <- .ctJuliaNumericVector
  npar <- length(est)
  started <- proc.time()[["elapsed"]]
  since <- function(from) proc.time()[["elapsed"]] - from
  stages <- c()
  placed <- module$ctsem_laplace_placed(cont)
  climb_sink <- if (is.null(sink)) NULL else function(text, kind = "update")
    sink(text, if (identical(kind, "done")) "update" else kind)
  x <- est
  hc <- NULL
  if (!stay) {
    res <- get(module$ctsem_laplace_placed_optimize(placed, num(est),
      maxiter = as.integer(.ctJuliaOr(control$placed_maxiter, 1000L)),
      gap_tol = as.numeric(tol), converge_tol = as.numeric(tol), newton = TRUE,
      certify = TRUE, lbfgs_pairs = "objective", progress = !is.null(climb_sink),
      progress_sink = climb_sink, progress_overwrite = .ctProgressOverwrite(verbose),
      progress_label = "Refined fit, rules placed"))
    stages["climb"] <- since(started)
    x <- as.numeric(res$minimizer)
    reuse <- .ctBackendHessianReuse(.ctBackendFitOptimcontrol(fit))
    hc <- as.numeric(res$hessian)
    hc <- if (length(hc) == npar^2 && all(is.finite(hc)) &&
        isTRUE(as.numeric(res$hessian_distance) <= reuse)) matrix(hc, npar, npar)
    record$optimiser <- list(iterations = as.integer(res$iterations),
      f_calls = as.integer(res$f_calls), g_calls = as.integer(res$g_calls),
      converged = isTRUE(res$converged), hessians = as.integer(res$newton_hessians),
      replacements = as.integer(res$replacements))
    record$continuation <- if (isTRUE(res$converged)) "converged" else "stopped"
  }
  if (is.null(hc)) {
    if (!is.null(sink)) sink(sprintf("Refined fit, rules placed: hessian"), "update")
    hessian_started <- proc.time()[["elapsed"]]
    hc <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_laplace_placed_hessian(placed, num(x)))), npar, npar), silent = TRUE)
    if (inherits(hc, "try-error")) hc <- NULL
    stages["hessian"] <- since(hessian_started)
  }
  final <- get(module$ctsem_laplace_placed_evaluate(placed, num(x), gradient = TRUE))
  units <- as.numeric(final$unit_loglik)
  record$route <- "placed"
  # A run asked to stay placed its rules at the estimate and moved nothing.
  record$status <- if (stay) "placed" else "continued"
  record$applied <- !stay
  se <- as.numeric(fit$estimate$se)
  delta <- x - est
  delta_se <- if (length(se) == npar) ifelse(se > 0, delta / se, NA_real_) else
    rep(NA_real_, npar)
  names(delta) <- names(delta_se) <- names(fit$estimate$se)
  record$delta <- delta
  record$delta_se <- delta_se
  record$material <- any(is.finite(delta_se)) &&
    isTRUE(max(abs(delta_se), na.rm = TRUE) >= control$material)
  record$loglik_quadrature <- sum(units)
  record$logposterior_quadrature <- as.numeric(final$value)
  record$objective <- as.numeric(final$value)
  fit$estimate$raw <- x
  fit <- report(fit, list(quadrature_units = units, quadrature = final$value,
    quadrature_subjects = as.numeric(final$subject_loglik)))
  usable <- is.matrix(hc) && all(is.finite(hc))
  record["hessian"] <- list(NULL)
  certification <- NULL
  if (usable) {
    certification_started <- proc.time()[["elapsed"]]
    gradient <- as.numeric(final$gradient)
    gap <- .ctBackendOptimGap(hc, gradient)
    saddle <- NULL
    if (isTRUE(gap$ok) && isTRUE(gap$nnegative > 0L)) {
      saddle <- .ctBackendSaddleProbe(module, placed, x, hc, ascent = gradient)
    }
    probe <- NULL
    if (isTRUE(gap$ok) && isTRUE(gap$nflat > 0L || gap$nnegative > 0L) &&
        isTRUE(gap$residual_norm > 0)) {
      out <- try(get(module$ctsem_flat_probe(placed, num(x), num(gap$residual))),
        silent = TRUE)
      if (!inherits(out, "try-error")) probe <- .ctBackendProbeFields(out, npar)
    }
    parnames <- try(.ctBackendRawParameterNames(fit, npar), silent = TRUE)
    if (inherits(parnames, "try-error")) parnames <- NULL
    certification <- .ctBackendCertificationRecord(gap, probe, tolerance = tol,
      saturated = isTRUE(fit$optim$saturated), overshot = FALSE,
      parnames = parnames, saddle = saddle)
    record$certification <- certification
    record$hessian <- hc
    stages["certification"] <- since(certification_started)
  }
  method <- .ctJuliaOr(method, .ctJuliaOr(fit$uncertainty$settings$method, "hessian"))
  own <- usable && method %in% c("hessian", "checked", "is")
  record$hessian_at <- if (own) "estimate" else "laplace_estimate"
  if (own) {
    uncertainty_started <- proc.time()[["elapsed"]]
    settings <- fit$uncertainty$settings
    finishsamples <- .ctJuliaOr(settings$finishsamples,
      if (!is.null(fit$estimate$rawposterior)) nrow(fit$estimate$rawposterior) else 1000L)
    fit$uncertainty <- list(hessian = hc, evaluated_at = x,
      certification = certification)
    fit <- .ctBackendUncertainty(fit, uncertainty = method,
      draws = .ctResolveDraws(method, "auto"),
      finishsamples = finishsamples, cores = chunks,
      control = .ctJuliaOr(uncertaintycontrol, .ctJuliaOr(settings$control, list())),
      verbose = verbose,
      lpg = .ctLaplaceContinueLpg(module, cont, objective = placed, unitsum = FALSE))
    record$draws <- "redrawn"
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("draws about the refined estimate, by '", method,
        "' on its objective with the rules placed at each point"),
      nodes = record$nodes)
    stages["uncertainty"] <- since(uncertainty_started)
  } else {
    post <- fit$estimate$rawposterior
    if (!is.null(post) && ncol(post) == npar) {
      fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
    }
    record$draws <- "recentred"
  }
  # What the run cost: the placed objective's calls and seconds, and the
  # continuation's (the screen's placement), with the hybrid's own error at
  # the screen.
  after <- get(module$ctsem_laplace_continuation_info(cont))
  record$error <- .ctLaplaceIntegrationError(after, unscored = record$wide)
  info <- get(module$ctsem_laplace_placed_info(placed))
  record$placed <- info
  record$evaluations <- c(values = as.integer(info$value_calls),
    gradients = as.integer(info$gradient_calls),
    placements = as.integer(after$recentres))
  record$engine_seconds <- c(values = as.numeric(info$seconds_values),
    gradients = as.numeric(info$seconds_gradients),
    placements = as.numeric(after$seconds_placements),
    hessian = as.numeric(info$seconds_hessian))
  record$seconds <- c(timed(), stages)
  fit$laplace$correction <- record
  fit
}
