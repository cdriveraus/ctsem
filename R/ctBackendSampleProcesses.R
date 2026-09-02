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
  .ctBackendPoolChains(fit, drawn, chains, warmup, draws, saveEffects)
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
      # The effect draws only when they were asked for -- they are the one
      # field here big enough for the bridge to notice -- but their summary
      # always, because a pooled fit that reported no random effects at all is
      # what dropping it produced.
      effects = s$sample$effects,
      effect_mean = s$sample$effect_mean,
      effect_sd = s$sample$effect_sd,
      divergent = s$sample$divergent,
      warmup_divergent = s$sample$warmup_divergent,
      saturated = s$sample$saturated,
      max_depth = s$sample$max_depth,
      stepsize = s$sample$stepsize,
      ebfmi = s$sample$ebfmi,
      accept = s$sample$accept,
      depth = s$sample$depth,
      energy = s$sample$energy)
  }, error = function(e) structure(list(), error = conditionMessage(e)))
}

# Pool the chains, then hand them to the ordinary assembler.
#
# What only the pool can say is R-hat and effective sample size, which are
# properties of the whole run rather than averages of per-chain values, and the
# effect summaries, which have to be recombined rather than concatenated.
# Everything after that -- where the draws go, what becomes the point estimate,
# which diagnostics warn, what class the object carries -- is what the
# single-process path already does, so this builds the engine's own result shape
# and calls `.ctBackendSampleAssemble()` rather than filling the fit in again by
# hand.
#
# It was written the other way first, and four things went missing in the copy:
# the diagnostics carried no class, so `print()` fell back to printing a list;
# the Laplace start was absent; the constrained draws still described the
# optimised fit's; and the effect summaries were dropped entirely, so every
# multi-chain fit -- which is the default -- reported no random effects at all.
#' @keywords internal
.ctBackendPoolChains <- function(fit, drawn, chains, warmup, draws,
  saveEffects = FALSE) {
  posteriors <- lapply(drawn, function(d) as.matrix(d$posterior))
  npar <- ncol(posteriors[[1]])
  if (!all(vapply(posteriors, ncol, integer(1)) == npar)) return(NULL)

  # Every chain must have contributed the same number of draws for the layout
  # below to hold, and an effective-sample-size target can break that: a chain
  # that reached `minEss` early stops before one that did not. Truncating to the
  # shortest is the honest repair -- these are all post-warmup draws from the
  # same stationary distribution, so dropping the tail of the longer chains
  # costs a little precision and nothing else, where refusing to pool would
  # discard every chain and sample the whole run again in this session.
  counts <- vapply(posteriors, nrow, integer(1))
  ndraws <- min(counts)
  if (ndraws < 1L) return(NULL)
  if (any(counts != ndraws)) {
    message("Chains returned ", paste(counts, collapse = ", "),
      " draws, so the first ", ndraws, " of each were pooled.")
  }
  rows <- function(x) {
    if (is.null(x) || !length(x)) return(NULL)
    x <- as.matrix(x)
    x[seq_len(ndraws), , drop = FALSE]
  }
  perdraw <- function(field) unlist(lapply(drawn, function(d) {
    v <- as.numeric(d[[field]])
    if (length(v) >= ndraws) v[seq_len(ndraws)] else v
  }))

  # `rbind` stacks chain 1's draws, then chain 2's, so the transpose below is
  # `ndim x (nchains * ndraws)` with column `(c-1)*ndraws + t` holding chain
  # `c`'s draw `t` -- the chain-major layout `ctsem_sample_diagnostics` indexes
  # and the assembler reshapes. Getting this wrong would not error -- it would
  # silently mix the chains and report R-hat over the mixture, which is always
  # reassuring.
  pooled <- do.call(rbind, lapply(posteriors, rows))
  colnames(pooled) <- colnames(posteriors[[1]])
  effectdraws <- lapply(drawn, function(d) rows(d$effects))
  pooledeffects <- if (all(vapply(effectdraws, is.matrix, logical(1))))
    do.call(rbind, effectdraws) else NULL

  # The effect summaries recombine rather than concatenate. The pooled mean is
  # the mean of the chains' means; the pooled variance is the within-chain sum
  # of squares plus the between-chain one, over the pooled degrees of freedom.
  # Averaging the chains' standard deviations instead would understate the
  # spread by exactly the part between chains -- the part R-hat is about.
  means <- do.call(rbind, lapply(drawn, function(d) as.numeric(d$effect_mean)))
  sds <- do.call(rbind, lapply(drawn, function(d) as.numeric(d$effect_sd)))
  effectmean <- numeric(0)
  effectsd <- numeric(0)
  if (!is.null(means) && ncol(means) > 0L && identical(dim(means), dim(sds))) {
    effectmean <- colMeans(means)
    total <- (ndraws - 1) * colSums(sds^2) +
      ndraws * colSums(sweep(means, 2, effectmean)^2)
    effectsd <- sqrt(total / max(1L, ndraws * nrow(means) - 1L))
  }
  neffects <- max(length(effectmean),
    if (is.null(pooledeffects)) 0L else ncol(pooledeffects))
  # The draws themselves can only be returned if every chain actually sent
  # them; asking the assembler to reshape to a width the matrix does not have
  # would interleave parameters and effects into plausible-looking nonsense.
  keepeffects <- isTRUE(saveEffects) && !is.null(pooledeffects) &&
    ncol(pooledeffects) == neffects

  module <- .ctJuliaModule(fit$model_spec$project)
  # Not `diag`, which would shadow `base::diag` for the rest of the function.
  pooldiag <- JuliaConnectoR::juliaGet(module$ctsem_sample_diagnostics(
    JuliaConnectoR::juliaPut(t(pooled)), as.integer(chains)))

  # The engine's own result shape, so that one assembler serves both paths.
  result <- list(
    draws = if (keepeffects) t(cbind(pooled, pooledeffects)) else t(pooled),
    npar = as.integer(npar),
    # Above `npar` whenever the chains summarised any effects, which is what
    # tells the assembler this was the joint sampler and not a marginal one.
    ndim = as.integer(npar + neffects),
    ndraws = as.integer(ndraws),
    rhat = as.numeric(pooldiag$rhat), ess = as.numeric(pooldiag$ess),
    # Summed across chains, because they count events; the step size and E-BFMI
    # are per chain and stay per chain.
    ndivergent = sum(vapply(drawn, function(d) as.integer(d$divergent), integer(1))),
    warmup_divergent = sum(vapply(drawn,
      function(d) as.integer(d$warmup_divergent), integer(1))),
    nsaturated = sum(vapply(drawn, function(d) as.integer(d$saturated), integer(1))),
    max_depth = max(vapply(drawn, function(d) as.integer(d$max_depth), integer(1))),
    stepsize = unlist(lapply(drawn, function(d) as.numeric(d$stepsize))),
    ebfmi = unlist(lapply(drawn, function(d) as.numeric(d$ebfmi))),
    accept = perdraw("accept"), depth = perdraw("depth"),
    energy = perdraw("energy"),
    effect_mean = effectmean, effect_sd = effectsd)

  out <- .ctBackendSampleAssemble(fit, result, npar, keepeffects,
    as.integer(chains), warmup, ndraws, fit$uncertainty$hessian,
    as.numeric(fit$estimate$raw))
  # Recorded after the fact because it changes nothing about the draws and
  # everything about how they were produced.
  out$uncertainty$settings$processes <- TRUE
  out$sample$processes <- TRUE
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
