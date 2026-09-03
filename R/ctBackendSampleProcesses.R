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
#' @param fit A `ctJuliaFit`, or the shell `ctFit(optimize = FALSE)` is about to
#'   fill. Carries the model and data a worker rebuilds its objective from.
#' @param target From [.ctBackendSampleTarget()]: what to sample, in a form that
#'   survives serialisation.
#' @param chains Number of chains, one per worker.
#' @param warmup,draws Per chain.
#' @param cores Total threads to divide among the workers.
#' @param handles Optional warmed pool from [.ctBackendWarmWorkers()]. When
#'   absent the workers are started here and the compile is paid in full.
#' @param control,saveEffects,seed,verbose As for [ctSample()].
#' @param progress Report chain progress from the parent while the workers
#'   run. Separate from `verbose` for the same reason `.ctBackendSampleEngine`
#'   keeps the two apart -- see there.
#' @return A fit with pooled draws and diagnostics, or `NULL` if the workers
#'   could not be used, in which case the caller samples in-process.
#' @keywords internal
.ctBackendSampleProcesses <- function(fit, target, chains, warmup, draws,
  cores = 1L, handles = NULL, control = list(), saveEffects = FALSE,
  seed = 1L, verbose = FALSE, progress = .ctVerboseOn(verbose)) {

  if (!.ctBackendCanWarm() || chains < 2L) return(NULL)
  if (!inherits(fit, "ctJuliaFit")) return(NULL)

  # Threads left over after one process per chain. A worker's own subject split
  # then uses them, which is the same nesting the in-process path does, except
  # that the chains no longer share an allocator.
  per_worker <- max(1L, as.integer(cores) %/% as.integer(chains))

  if (is.null(handles)) {
    handles <- .ctBackendWarmWorkers(fit, workers = chains,
      values = target$estimate)
    if (is.null(handles)) return(NULL)
  }
  .ctBackendWarmWait(handles, verbose = verbose)

  # Each chain runs in its own process, so its printed output sits in that
  # process's own stdout buffer and only reaches the parent -- all at once,
  # after the fact -- when `future::value()` collects it. `ctSample(verbose =
  # TRUE)` under `processes = TRUE` used to print nothing at all for exactly
  # this reason: the reporting existed, in the worker, and had nowhere to go
  # until the run was already over.
  #
  # The fix is to report from the parent instead of hoping a worker's console
  # output arrives. Each worker's engine call is given a callback -- the same
  # `progress_callback` a GUI would use, see `.ctBackendSampleEngine` -- that
  # writes a one-line snapshot to a small file rather than printing, and the
  # parent polls those files while it waits and prints one line per chain.
  # That is the same information the single-process path prints live, just
  # relayed through a file because a process boundary is in the way.
  report <- isTRUE(progress)
  progress_files <- if (report) vapply(seq_len(chains), function(i)
    tempfile(pattern = sprintf("ctsem_sample_chain%d_", i), fileext = ".progress"),
    character(1)) else NULL
  if (report) on.exit(unlink(progress_files, force = TRUE), add = TRUE)

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
    chain_file <- if (report) progress_files[k] else NULL
    tryCatch(
      future::future(ctsem:::.ctBackendSampleOneChain(fit, target, warmup,
        draws, per_worker, control, saveEffects, as.integer(seed) + k - 1L,
        progress_file = chain_file),
        seed = TRUE),
      error = function(e) NULL)
  })
  if (report) {
    .ctBackendReportProcesses(results, progress_files, chains = chains)
  }
  drawn <- lapply(results, function(h) {
    if (is.null(h)) return(NULL)
    tryCatch(future::value(h), error = function(e) NULL)
  })
  # A worker that failed returns an empty list carrying an `error` attribute,
  # not NULL, so testing for NULL alone would let it through and the failure
  # would surface later as an empty matrix in the pooling.
  ok <- vapply(drawn, function(d)
    !is.null(d) && !is.null(d$draws) && length(d$draws) > 0,
    logical(1))
  if (sum(ok) < chains) {
    warning(sum(!ok), " of ", chains, " chains failed in their worker ",
      "process, so sampling fell back to this session. The first error was: ",
      .ctBackendFirstError(drawn), call. = FALSE)
    return(NULL)
  }
  .ctBackendPoolChains(fit, target, drawn, chains, warmup, draws, saveEffects)
}

# One chain, in a worker.
#
# Deliberately the shared runner's own engine call: a second sampler for the
# process path would be a second thing to keep correct, and it would have to be
# told the same things anyway. What the worker does not run is the assembly --
# constraining a single chain's draws only to throw them away when the pool is
# assembled is work nobody reads.
#
# It used to call `ctSample()`, which fixed the target as well as the code: the
# joint entry, and a refusal for any fit without a Laplace spec. That is right
# for `ctSample()`'s own callers and wrong for two of the three routes
# `ctFit(optimize = FALSE)` can take, which is why the target now travels
# explicitly.
#
# `progress_file`, when given, replaces `control$callback` for this call: a
# user's own callback is an R closure over the parent session (a plot device,
# a Shiny reactive) and calling it from here would try to reach across a
# process boundary that does not carry it. What can cross is a path, and the
# parent polls what gets written there -- see `.ctBackendReportProcesses`.
#' @keywords internal
.ctBackendSampleOneChain <- function(fit, target, warmup, draws, threads,
  control, saveEffects, seed, progress_file = NULL) {
  tryCatch({
    if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = threads)
    callback <- if (is.null(progress_file)) NULL else
      .ctBackendProgressFileWriter(progress_file)
    result <- suppressWarnings(suppressMessages(
      .ctBackendSampleEngine(fit, target, chains = 1L, warmup = warmup,
        draws = draws, cores = threads, saveEffects = saveEffects, seed = seed,
        control = control, verbose = FALSE, progress = FALSE,
        callback = callback)))
    .ctBackendChainResult(result)
  }, error = function(e) structure(list(), error = conditionMessage(e)))
}

# A callback that writes one line rather than printing one -- the worker's
# stdout is not read live, so a callback that `cat()`ed here would be exactly
# as invisible as the printed progress line this whole mechanism exists to
# work around.
#
# Written to a temp path and renamed into place, so the parent, reading
# concurrently, never sees a half-written line: `file.rename` within one
# filesystem is atomic, a plain `writeLines` to the final path is not.
# Wrapped in `tryCatch` because a reporting write must never be the thing that
# fails a chain -- a full disk or a deleted temp directory should cost a
# missed update, not the sample.
#' @keywords internal
.ctBackendProgressFileWriter <- function(path) {
  tmp <- paste0(path, ".tmp")
  function(phase, iteration, total, logp, divergent) {
    tryCatch({
      writeLines(paste(as.character(phase), as.integer(iteration),
        as.integer(total), sprintf("%.6f", as.numeric(logp)),
        as.integer(divergent)), tmp)
      file.rename(tmp, path)
    }, error = function(e) NULL)
    NULL
  }
}

# The counterpart read: one line, back into its fields, or `NULL` for
# anything that does not parse -- a file not yet written, or caught mid-write
# despite the rename (a stale reader on a slow network share, say).
#' @keywords internal
.ctBackendReadProgressFile <- function(path) {
  if (!file.exists(path)) return(NULL)
  line <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) NULL)
  if (is.null(line) || !length(line) || !nzchar(line)) return(NULL)
  parts <- strsplit(line, "\\s+")[[1]]
  if (length(parts) < 5L) return(NULL)
  iteration <- suppressWarnings(as.integer(parts[2]))
  total <- suppressWarnings(as.integer(parts[3]))
  logp <- suppressWarnings(as.numeric(parts[4]))
  divergent <- suppressWarnings(as.integer(parts[5]))
  if (anyNA(c(iteration, total, logp, divergent))) return(NULL)
  list(phase = parts[1], iteration = iteration, total = total, logp = logp,
    divergent = divergent)
}

#' Report per-chain progress from the parent while chain processes run
#'
#' Polls each chain's progress file on a short interval and prints one line
#' per chain still running -- the shape asked for when chains are processes:
#' each worker's own printed progress sits in output that never reaches the
#' parent until the chain is already done, so the parent reports instead,
#' from what the workers wrote rather than from what they printed.
#'
#' Not overwritten in place. `CTSEMProgress` overwrites a single line because
#' it owns the whole of what is on it; here several chains share the console
#' and a later one finishing does not mean an earlier one's last line should
#' vanish. A short block of chains, printed occasionally, is simple, will not
#' garble on any terminal, and reads fine at the couple-of-seconds cadence
#' this polls at -- faster would not show anything a chain-level report needs.
#'
#' @param results Future handles from [.ctBackendSampleProcesses()], one per
#'   chain, possibly containing `NULL` for a chain that never started.
#' @param progress_files One path per chain, written by
#'   [.ctBackendProgressFileWriter()].
#' @param chains Number of chains.
#' @param interval Seconds between polls.
#' @return `NULL`, invisibly. Called for its printing.
#' @keywords internal
.ctBackendReportProcesses <- function(results, progress_files, chains,
  interval = 2) {
  now <- Sys.time()
  phase_started <- rep(now, chains)
  phase_seen <- rep(NA_character_, chains)
  repeat {
    resolved <- vapply(results, function(h) is.null(h) || future::resolved(h),
      logical(1))
    lines <- character(0)
    for (k in seq_len(chains)) {
      if (resolved[k]) next
      info <- .ctBackendReadProgressFile(progress_files[k])
      if (is.null(info)) next
      if (is.na(phase_seen[k]) || !identical(phase_seen[k], info$phase)) {
        phase_started[k] <- Sys.time()
        phase_seen[k] <- info$phase
      }
      elapsed <- as.numeric(difftime(Sys.time(), phase_started[k], units = "secs"))
      rate <- if (info$iteration > 0 && elapsed > 0) info$iteration / elapsed else 0
      eta <- if (rate > 0 && info$total > info$iteration)
        paste0(" | ", .ctDuration((info$total - info$iteration) / rate),
          " at this rate") else ""
      lines <- c(lines, sprintf("  chain %d/%d: %-8s %5d/%-5d%s | logp %.2f",
        k, chains, info$phase, info$iteration, info$total, eta, info$logp))
    }
    if (length(lines)) cat(paste(lines, collapse = "\n"), "\n", sep = "")
    if (all(resolved)) break
    Sys.sleep(interval)
  }
  invisible(NULL)
}

# What a chain sends home.
#
# The engine's own result, minus the two things the pool recomputes -- R-hat and
# effective sample size are properties of the whole run and cannot be averaged
# from per-chain values -- and minus anything that would be a second copy of the
# model. Returning the whole fit would send the data and the model back across
# for every chain, having already sent them out.
#' @keywords internal
.ctBackendChainResult <- function(result) {
  ndraws <- as.integer(result$ndraws)
  list(
    # `kept x ndraws`, which is the engine's own layout, so pooling the chains
    # is a `cbind` and the assembler reshapes the pool exactly as it reshapes a
    # single-process result.
    draws = matrix(as.numeric(result$draws), ncol = ndraws),
    npar = as.integer(result$npar), ndim = as.integer(result$ndim),
    ndraws = ndraws,
    ndivergent = as.integer(result$ndivergent),
    warmup_divergent = as.integer(result$warmup_divergent),
    nsaturated = as.integer(result$nsaturated),
    max_depth = as.integer(result$max_depth),
    stepsize = as.numeric(result$stepsize),
    ebfmi = as.numeric(result$ebfmi),
    accept = as.numeric(result$accept),
    depth = as.numeric(result$depth),
    energy = as.numeric(result$energy),
    # Always, even unsaved: a pooled fit that reported no random effects at all
    # is what dropping these produced.
    effect_mean = as.numeric(result$effect_mean),
    effect_sd = as.numeric(result$effect_sd))
}

# Pool the chains into one engine result, then hand it to the ordinary assembler.
#
# What only the pool can say is R-hat and effective sample size, and the effect
# summaries, which have to be recombined rather than concatenated. Everything
# after that -- where the draws go, what becomes the point estimate, which
# diagnostics warn, what class the object carries -- is what the single-process
# path already does, so this produces the engine's own shape and calls
# `.ctBackendSampleAssemble()` rather than filling the fit in again by hand.
#
# It was written the other way first, and four things went missing in the copy:
# the diagnostics carried no class, so `print()` fell back to printing a list;
# the Laplace start was absent; the constrained draws still described the
# optimised fit's; and the effect summaries were dropped entirely, so every
# multi-chain fit -- which is the default -- reported no random effects at all.
#' @keywords internal
.ctBackendPoolChains <- function(fit, target, drawn, chains, warmup, draws,
  saveEffects = FALSE) {
  mats <- lapply(drawn, function(d) as.matrix(d$draws))
  npar <- as.integer(target$npar)
  kept <- nrow(mats[[1]])
  if (!all(vapply(mats, nrow, integer(1)) == kept) || kept < npar) return(NULL)

  # Every chain must have contributed the same number of draws for the layout
  # below to hold, and an effective-sample-size target can break that: a chain
  # that reached `minEss` early stops before one that did not. Truncating to the
  # shortest is the honest repair -- these are all post-warmup draws from the
  # same stationary distribution, so dropping the tail of the longer chains
  # costs a little precision and nothing else, where refusing to pool would
  # discard every chain and sample the whole run again in this session.
  counts <- vapply(mats, ncol, integer(1))
  ndraws <- min(counts)
  if (ndraws < 1L) return(NULL)
  if (any(counts != ndraws)) {
    message("Chains returned ", paste(counts, collapse = ", "),
      " draws, so the first ", ndraws, " of each were pooled.")
  }
  perdraw <- function(field) unlist(lapply(drawn, function(d) {
    v <- as.numeric(d[[field]])
    if (length(v) >= ndraws) v[seq_len(ndraws)] else v
  }))

  # `cbind` puts chain 1's draws first, then chain 2's, which is the chain-major
  # `ndim x (nchains * ndraws)` layout `ctsem_sample_diagnostics` indexes and the
  # assembler reshapes: column `(c-1)*ndraws + t` holds chain `c`'s draw `t`.
  # Getting this wrong would not error -- it would silently mix the chains and
  # report R-hat over the mixture, which is always reassuring.
  pooled <- do.call(cbind,
    lapply(mats, function(m) m[, seq_len(ndraws), drop = FALSE]))
  ndim <- as.integer(drawn[[1]]$ndim)
  if (!isTRUE(is.finite(ndim))) ndim <- kept
  keepeffects <- isTRUE(saveEffects) && kept > npar
  # A chain that sent effect draws nobody asked for: the assembler reshapes to
  # `npar` rows in that case, so the extra rows have to go rather than be
  # interleaved into plausible-looking nonsense.
  if (!keepeffects && kept > npar) pooled <- pooled[seq_len(npar), , drop = FALSE]

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

  module <- .ctJuliaModule(fit$model_spec$project)
  # The population block alone: R-hat over every random effect as well would
  # cost more than it says, and the assembler reads only the first `npar`.
  # Not `diag`, which would shadow `base::diag` for the rest of the function.
  pooldiag <- JuliaConnectoR::juliaGet(module$ctsem_sample_diagnostics(
    JuliaConnectoR::juliaPut(pooled[seq_len(npar), , drop = FALSE]),
    as.integer(chains)))

  result <- list(
    draws = pooled, npar = npar, ndim = ndim, ndraws = ndraws,
    rhat = as.numeric(pooldiag$rhat), ess = as.numeric(pooldiag$ess),
    # Summed across chains, because they count events; the step size and E-BFMI
    # are per chain and stay per chain.
    ndivergent = sum(vapply(drawn, function(d) as.integer(d$ndivergent), integer(1))),
    warmup_divergent = sum(vapply(drawn,
      function(d) as.integer(d$warmup_divergent), integer(1))),
    nsaturated = sum(vapply(drawn, function(d) as.integer(d$nsaturated), integer(1))),
    max_depth = max(vapply(drawn, function(d) as.integer(d$max_depth), integer(1))),
    stepsize = unlist(lapply(drawn, function(d) as.numeric(d$stepsize))),
    ebfmi = unlist(lapply(drawn, function(d) as.numeric(d$ebfmi))),
    accept = perdraw("accept"), depth = perdraw("depth"),
    energy = perdraw("energy"),
    effect_mean = effectmean, effect_sd = effectsd)

  out <- .ctBackendSampleAssemble(fit, result, npar, keepeffects,
    as.integer(chains), warmup, ndraws, target$hessian,
    target$estimate[seq_len(npar)])
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
