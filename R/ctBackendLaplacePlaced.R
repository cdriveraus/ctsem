# The refined fit with its flagged units' rules placed at every evaluation
# (`laplace_placed.jl`): `ctsem_laplace_placed_optimize` climbs the placed
# objective from the Laplace estimate, L-BFGS and then the Newton finish on a
# Hessian of differenced exact gradients, and the uncertainty is drawn on that
# same objective, whose gradient is its value's. Units on frozen draws are
# placed again only where their weights' effective size has collapsed, the
# run's L-BFGS memory dropped. Option B of
# review/SURVEY-quadrature-2026-10-30.md, which replaced the rounds
# (review/PLACED-rules-race-2026-10-30.md).

# The climb alone, shared by the fit and `ctLaplaceCheck(refine = TRUE)`: the
# placed objective `placed` from `est`, to the fit's certification tolerance
# `tol`. Returns the engine's result (`ctsem_optimize`'s, with `replacements`).
.ctLaplacePlacedClimb <- function(module, placed, est, tol, control, verbose = 0L,
  sink = NULL, certify = TRUE) {
  climb_sink <- if (is.null(sink)) NULL else function(text, kind = "update")
    sink(text, if (identical(kind, "done")) "update" else kind)
  .ctJuliaGet(module$ctsem_laplace_placed_optimize(placed,
    .ctJuliaNumericVector(est),
    maxiter = as.integer(.ctJuliaOr(control$maxiter, 1000L)),
    gap_tol = as.numeric(tol), converge_tol = as.numeric(tol), newton = TRUE,
    certify = isTRUE(certify), lbfgs_pairs = "objective",
    progress = !is.null(climb_sink), progress_sink = climb_sink,
    progress_overwrite = .ctProgressOverwrite(verbose),
    progress_label = "Refined fit, rules placed"))
}

# Called by `.ctLaplaceContinue()` once the continuation `cont` is placed at the
# estimate `est` and `record` describes that placement; `report` writes the
# quadrature log likelihood into the fit and `timed` the wall seconds by stage.
# `stay` moves nothing: the rules placed at the estimate, its Hessian and the
# uncertainty drawn there. The guard reads the one objective at the start and
# the end: a change past `guard` nats,
# or `guard_per_subject` a subject, reverts the run, and an end lower than the
# start by `value_tol` keeps the Laplace optimum. Returns the fit.
.ctLaplacePlacedRefine <- function(fit, module, cont, est, record, control, tol,
  stay, chunks, verbose, sink, method, uncertaintycontrol, report, timed) {
  get <- .ctJuliaGet
  num <- .ctJuliaNumericVector
  npar <- length(est)
  nsubjects <- length(fit$model_spec$subject_starts)
  started <- proc.time()[["elapsed"]]
  since <- function(from) proc.time()[["elapsed"]] - from
  stages <- c()
  placed <- module$ctsem_laplace_placed(cont)
  start <- get(module$ctsem_laplace_placed_evaluate(placed, num(est), gradient = FALSE))
  # What a fit reports at a point: every unit by quadrature, the flagged units
  # by their rules and the rest by the rule the screen would place for each
  # (`ctsem_laplace_placed_quadrature`), as the rounds reported theirs.
  quadrature_at <- function(at) {
    e <- get(module$ctsem_laplace_placed_quadrature(placed, num(at)))
    list(quadrature_units = as.numeric(e$unit_loglik), quadrature = as.numeric(e$value),
      quadrature_subjects = as.numeric(e$subject_loglik))
  }
  # What the screen left on Laplace and what it flagged, from the placement.
  after <- get(module$ctsem_laplace_continuation_info(cont))
  record$route <- "placed"
  record$flagged_units <- as.integer(after$flagged)[seq_len(after$nflagged)]
  record$bar <- as.numeric(after$bar)
  record$residual <- as.numeric(after$residual)
  # Units on importance draws, frozen for the run, and their draws.
  record$importance <- c(units = as.integer(after$nimportance),
    draws = as.integer(after$importance_draws))
  accounted <- function(record) {
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
    record
  }
  x <- est
  hc <- NULL
  # Where the Hessian was formed: the finish's own point, within the reuse
  # distance of the estimate, or the estimate.
  hessian_at <- est
  if (!stay) {
    res <- .ctLaplacePlacedClimb(module, placed, est, tol, control, verbose, sink)
    stages["climb"] <- since(started)
    x <- as.numeric(res$minimizer)
    reuse <- .ctBackendHessianReuse(.ctBackendFitOptimcontrol(fit))
    hc <- as.numeric(res$hessian)
    hc <- if (length(hc) == npar^2 && all(is.finite(hc)) &&
        isTRUE(as.numeric(res$hessian_distance) <= reuse)) matrix(hc, npar, npar)
    hessian_at <- if (is.null(hc)) x else as.numeric(res$hessian_evaluated_at)
    record$optimiser <- list(iterations = as.integer(res$iterations),
      f_calls = as.integer(res$f_calls), g_calls = as.integer(res$g_calls),
      converged = isTRUE(res$converged), hessians = as.integer(res$newton_hessians),
      replacements = as.integer(res$replacements))
    record$continuation <- if (isTRUE(res$converged)) "converged" else "stopped"
  }
  final <- get(module$ctsem_laplace_placed_evaluate(placed, num(x), gradient = TRUE))
  change <- as.numeric(final$value) - as.numeric(start$value)
  limit <- max(control$guard, control$guard_per_subject * nsubjects)
  record$guard <- list(limit = limit, change = change,
    fired = !stay && (!is.finite(change) || abs(change) > limit))
  worse <- !stay && is.finite(change) && change < -as.numeric(control$value_tol)
  record$ended_lower <- worse
  if (isTRUE(record$guard$fired) || worse) {
    if (isTRUE(record$guard$fired)) {
      if (!is.null(sink)) sink("", "break")
      warning(sprintf(paste0("The refined fit moved the quadrature objective by ",
        "%.3g nats, which no approximation error of this size explains; the fit ",
        "keeps the Laplace optimum, and the rejected point is ",
        "fit$laplace$correction$rejected_estimate."), change), call. = FALSE)
      record$rejected_estimate <- x
      record$status <- "reverted"
    } else {
      record$status <- "no_gain"
    }
    record$applied <- FALSE
    q <- quadrature_at(est)
    record$loglik_quadrature <- sum(q$quadrature_units)
    record$logposterior_quadrature <- q$quadrature
    record$gap_reported <- q$quadrature - as.numeric(after$laplace)
    fit <- report(fit, q)
    fit$laplace$correction <- accounted(record)
    return(fit)
  }
  if (is.null(hc)) {
    if (!is.null(sink)) sink("Refined fit, rules placed: hessian", "update")
    hessian_started <- proc.time()[["elapsed"]]
    hc <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_laplace_placed_hessian(placed, num(x)))), npar, npar), silent = TRUE)
    if (inherits(hc, "try-error")) hc <- NULL
    stages["hessian"] <- since(hessian_started)
  }
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
  q <- quadrature_at(x)
  record$loglik_quadrature <- sum(q$quadrature_units)
  record$logposterior_quadrature <- q$quadrature
  record$gap_reported <- q$quadrature - as.numeric(after$laplace)
  record$objective <- as.numeric(final$value)
  fit$estimate$raw <- x
  fit <- report(fit, q)
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
    fit$uncertainty <- list(hessian = hc, evaluated_at = hessian_at,
      certification = certification)
    fit <- .ctBackendUncertainty(fit, uncertainty = method,
      draws = .ctResolveDraws(method, "auto"),
      finishsamples = finishsamples, cores = chunks,
      control = .ctJuliaOr(uncertaintycontrol, .ctJuliaOr(settings$control, list())),
      verbose = verbose,
      lpg = .ctLaplacePlacedLpg(module, placed))
    record$draws <- "redrawn"
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("draws about the refined estimate, by '", method,
        "' on its objective with the rules placed at each point"),
      nodes = record$nodes,
      # What the density the draws came from is known to, in nats.
      error = .ctLaplaceIntegrationError(after, unscored = record$wide))
    stages["uncertainty"] <- since(uncertainty_started)
  } else {
    post <- fit$estimate$rawposterior
    if (!is.null(post) && ncol(post) == npar) {
      fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
    }
    record$draws <- "recentred"
  }
  fit$laplace$correction <- accounted(record)
  fit
}

# A log-posterior function on the placed objective `placed`, in the `lpgFunc`
# contract `.ctBackendUncertainty()` takes, with the engine, batch and
# gradient-batch attributes `.ctBackendLpgFunc()`'s density carries, so that
# 'checked' runs DAIS and its Stein check here too and draws cross the bridge
# a batch at a time. Every point is evaluated with the flagged units' rules
# placed there: read where they were placed at the estimate, the rules said
# less the further a draw moved, and gD1's draws put a population sd 2.9
# posterior sds high. It carries no 'unitsum' mark, so no coreset stands in
# for its units (`.ctCoresetLpg()`): the engine's coresets are the Laplace
# objective's.
.ctLaplacePlacedLpg <- function(module, placed) {
  engine <- function() list(module = module, objective = placed)
  batches <- .ctBackendEngineBatches(engine)
  function(fit, gradient = TRUE) {
    wantgrad <- isTRUE(gradient)
    fn <- function(parm) {
      at <- matrix(as.numeric(parm), nrow = 1L)
      if (!wantgrad) return(batches$batch(at))
      result <- batches$gradbatch(at)
      value <- result$value[1L]
      grad <- as.numeric(result$gradient)
      # `.ctBackendLpgFunc()`'s guard: a point the objective cannot use is a
      # large finite penalty with no gradient.
      if (!is.finite(value) || value <= -1e100 || any(!is.finite(grad))) {
        value <- -1e100
        grad <- rep(0, length(parm))
      }
      attributes(value) <- list(gradient = grad)
      value
    }
    attr(fn, "engine") <- engine
    attr(fn, "batch") <- batches$batch
    attr(fn, "gradbatch") <- batches$gradbatch
    fn
  }
}
