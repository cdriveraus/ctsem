# Chains in separate processes.
#
# Chains share nothing: no state, no workspace, no barrier, and they only meet
# when the draws are pooled at the end. That makes them the axis worth putting
# in separate processes, where the subject split is not -- units have to meet at
# a barrier on every gradient, and threads do that for free.
#
# What a worker costs is a fresh Julia and a fresh compile for the model's
# dimensions, measured at 26-43 s and unavoidable per process (see
# `inst/julia/ContinuousTimeSEM/src/workspace_buffers.jl` for why). That is
# affordable only because it need not be paid when sampling starts:
# `.ctBackendWarmWorkers` sets the workers compiling during the optimisation
# that has to run first anyway. Measured, workers warmed alongside a 39.8 s
# optimisation were ready with 0.0 s of waiting.
#
# Each worker runs one chain through the ordinary `ctSample` path, so there is
# no second sampler implementation to keep in step with the first. The parent
# pools the draws and recomputes R-hat and effective size over all of them,
# through the same Julia routine the single-process path uses -- those are
# properties of the whole run and cannot be averaged from per-chain values.

#' Run each chain in its own process
#'
#' @param fit A `ctJuliaFit`, already optimised.
#' @param chains Number of chains, one per worker.
#' @param warmup,draws Per chain.
#' @param cores Total threads to divide among the workers.
#' @param handles Optional warmed pool from [.ctBackendWarmWorkers()]. When
#'   absent the workers are started here and the compile is paid in full.
#' @param control,saveEffects,seed,verbose As for [ctSample()].
#' @return A fit with pooled draws and diagnostics, or `NULL` if the workers
#'   could not be used, in which case the caller samples in-process.
#' @keywords internal
.ctBackendSampleProcesses <- function(fit, chains, warmup, draws, cores = 1L,
  handles = NULL, control = list(), saveEffects = FALSE, seed = 1L,
  verbose = FALSE) {

  if (!.ctBackendCanWarm() || chains < 2L) return(NULL)
  if (!inherits(fit, "ctJuliaFit")) return(NULL)

  # Threads left over after one process per chain. A worker's own subject split
  # then uses them, which is the same nesting the in-process path does, except
  # that the chains no longer share an allocator.
  per_worker <- max(1L, as.integer(cores) %/% as.integer(chains))

  if (is.null(handles)) {
    handles <- .ctBackendWarmWorkers(fit, workers = chains,
      values = fit$estimate$raw)
    if (is.null(handles)) return(NULL)
  }
  .ctBackendWarmWait(handles, verbose = verbose)

  # `seed + k - 1`, which makes this path reproduce the in-process one exactly.
  #
  # The engine gives chain `c` the stream `Xoshiro(seed + c)`. A worker runs a
  # single chain, so its chain is `c = 1` and it draws `Xoshiro(S + 1)` from
  # whatever seed `S` it was handed. Setting `S = seed + k - 1` makes worker `k`
  # draw `Xoshiro(seed + k)` -- the stream in-process chain `k` would have used.
  # Same data, same start, same metric, same stream, so the draws come back
  # identical element for element.
  #
  # That is worth more than tidiness. Whether a user sampled in one process or
  # four should not change their numbers, and it makes the pooling *testable*:
  # pool both ways and compare. A chain-major layout error is otherwise
  # invisible, because mixing the chains still yields plausible means and a
  # reassuring R-hat -- it is the one bug here that hides rather than shouts.
  #
  # Distinct seeds per chain remain essential either way: chains sharing a
  # stream are the same chain, and R-hat over copies of one chain is 1 however
  # wrong they are.
  #
  # **Reproduction is close, not exact, and that is expected.** Measured against
  # the in-process path on the same fit and seed, with `warmup = 0` so the first
  # kept draw sits as near the shared start as the sampler ever gets: the first
  # draw differs by 6.2e-10 and the twelfth by 6.7e-09. With a normal warmup the
  # difference reaches order 1, which is what an early run of this comparison
  # reported and misread as a failure.
  #
  # The residue is process-local numerical state, not a fault in the pooling.
  # Each unit's mode comes from an inner Newton solve warm-started from whatever
  # the objective last held, and the parent carries an optimisation's worth of
  # history where a fresh worker carries one warm-up evaluation. Both land on
  # the same mode to solver tolerance rather than to the last bit, and NUTS
  # amplifies the difference: it is chaotic, so 1e-10 becomes order 1 within a
  # few dozen transitions.
  #
  # What the measurement does establish is the part that could have been wrong
  # and would not have announced itself. A mismatched stream or a chain-major
  # layout error would show at the *first* draw, at the scale of the posterior's
  # own width -- order 1, not 1e-10. Neither does.
  results <- lapply(seq_len(chains), function(k) {
    tryCatch(
      future::future(ctsem:::.ctBackendSampleOneChain(fit, warmup, draws,
        per_worker, control, saveEffects, as.integer(seed) + k - 1L),
        seed = TRUE),
      error = function(e) NULL)
  })
  drawn <- lapply(results, function(h) {
    if (is.null(h)) return(NULL)
    tryCatch(future::value(h), error = function(e) NULL)
  })
  # A worker that failed returns an empty list carrying an `error` attribute,
  # not NULL, so testing for NULL alone would let it through and the failure
  # would surface later as an empty matrix in the pooling.
  ok <- vapply(drawn, function(d)
    !is.null(d) && !is.null(d$posterior) && length(d$posterior) > 0,
    logical(1))
  if (sum(ok) < chains) {
    warning(sum(!ok), " of ", chains, " chains failed in their worker ",
      "process, so sampling fell back to this session. The first error was: ",
      .ctBackendFirstError(drawn), call. = FALSE)
    return(NULL)
  }
  .ctBackendPoolChains(fit, drawn, chains, warmup, draws)
}

# One chain, in a worker. Deliberately the ordinary entry point: a separate
# sampler for the process path would be a second thing to keep correct.
#' @keywords internal
.ctBackendSampleOneChain <- function(fit, warmup, draws, threads, control,
  saveEffects, seed) {
  tryCatch({
    if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = threads)
    s <- suppressWarnings(suppressMessages(ctsem::ctSample(fit, chains = 1L,
      warmup = warmup, draws = draws, cores = threads, seed = seed,
      saveEffects = saveEffects, control = control)))
    # Only what pooling needs. Returning the whole fit would send the data and
    # the model back across for every chain, having already sent them out.
    list(posterior = s$estimate$rawposterior,
      divergent = s$sample$divergent,
      warmup_divergent = s$sample$warmup_divergent,
      saturated = s$sample$saturated,
      max_depth = s$sample$max_depth,
      stepsize = s$sample$stepsize,
      ebfmi = s$sample$ebfmi,
      accept = s$sample$accept)
  }, error = function(e) structure(list(), error = conditionMessage(e)))
}

# Pool the chains and recompute what only the pool can say.
#' @keywords internal
.ctBackendPoolChains <- function(fit, drawn, chains, warmup, draws) {
  # `saveEffects` does not reach here: `rawposterior` holds the population part
  # only, whether or not the effects were kept, so pooling is the same either
  # way.
  posteriors <- lapply(drawn, function(d) as.matrix(d$posterior))
  npar <- ncol(posteriors[[1]])
  if (!all(vapply(posteriors, ncol, integer(1)) == npar)) return(NULL)

  # `rbind` stacks chain 1's draws, then chain 2's, so the transpose below is
  # `npar x (nchains * ndraws)` with column `(c-1)*ndraws + t` holding chain
  # `c`'s draw `t` -- the chain-major layout `ctsem_sample_diagnostics` indexes.
  # Every chain must have contributed the same number of draws for that to hold.
  if (length(unique(vapply(posteriors, nrow, integer(1)))) != 1L) return(NULL)
  pooled <- do.call(rbind, posteriors)
  colnames(pooled) <- colnames(posteriors[[1]])

  # Chain-major, `ndim x (nchains * ndraws)`, which is the layout
  # `ctsem_sample_diagnostics` documents and the one the single-process path
  # hands it. Getting this wrong would not error -- it would silently mix the
  # chains and report R-hat over the mixture, which is always reassuring.
  module <- .ctJuliaModule(fit$model_spec$project)
  # Not `diag`: that shadows `base::diag`, which is called a few lines below to
  # take standard errors off the covariance, and the failure would read as
  # "could not find function" rather than as anything about names.
  pooldiag <- JuliaConnectoR::juliaGet(module$ctsem_sample_diagnostics(
    JuliaConnectoR::juliaPut(t(pooled)), as.integer(chains)))

  out <- fit
  out$estimate$laplace_raw <- as.numeric(fit$estimate$raw)
  out$estimate$rawposterior <- pooled
  out$estimate$raw <- as.numeric(colMeans(pooled))
  out$estimate$cov <- stats::cov(pooled)
  out$estimate$se <- sqrt(diag(out$estimate$cov))
  out$uncertainty <- list(method = "sampling", hessian = fit$uncertainty$hessian,
    settings = list(chains = chains, warmup = warmup, draws = draws,
      processes = TRUE))
  out$sample <- list(
    chains = chains, warmup = warmup, draws = draws,
    rhat = stats::setNames(as.numeric(pooldiag$rhat)[seq_len(npar)], colnames(pooled)),
    ess = stats::setNames(as.numeric(pooldiag$ess)[seq_len(npar)], colnames(pooled)),
    # Summed across chains, because they count events; the step size and E-BFMI
    # are per chain and stay per chain.
    divergent = sum(vapply(drawn, function(d) as.integer(d$divergent), integer(1))),
    warmup_divergent = sum(vapply(drawn,
      function(d) as.integer(d$warmup_divergent), integer(1))),
    saturated = sum(vapply(drawn, function(d) as.integer(d$saturated), integer(1))),
    max_depth = max(vapply(drawn, function(d) as.integer(d$max_depth), integer(1))),
    stepsize = unlist(lapply(drawn, function(d) as.numeric(d$stepsize))),
    ebfmi = unlist(lapply(drawn, function(d) as.numeric(d$ebfmi))),
    accept = unlist(lapply(drawn, function(d) as.numeric(d$accept))),
    processes = TRUE)
  out$sample$unidentified <- .ctBackendFlatParameters(out, npar)
  .ctSampleWarn(out$sample)
  out
}

# The first thing that actually went wrong, for the fallback warning.
#' @keywords internal
.ctBackendFirstError <- function(drawn) {
  for (d in drawn) {
    msg <- attr(d, "error")
    if (!is.null(msg)) return(msg)
  }
  "not reported"
}
