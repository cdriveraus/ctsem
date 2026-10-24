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
# There are at most `cores` workers -- the call's own ceiling, 2 unless asked
# -- so four chains at the default run two to a worker, as one engine call
# whose chains keep the seeds they would have in this session. One worker per
# chain used to start four processes, and four cores, whatever `cores` said.
#
# Each worker runs the ordinary sampling engine call, so there is no second
# sampler implementation to keep in step with the first. The run stops by one
# rule over every chain, as it does in one session: after each batch a worker
# hands its draws to the parent through a file and waits; the parent pools
# every worker's batch, asks the engine's `ctsem_sample_verdict` how many more
# draws the chains need, and shares them out in proportion to each worker's
# speed, so no worker waits on another. The workers' shares of `cores` need
# not be equal -- 6 over 4 chains is 2, 2, 1, 1 -- and their chains end at
# different lengths, which the engine's diagnostics weight by length. Workers
# that stopped on their own share of the target, as they once did, returned
# chains the pool then cut to the shortest: 230 of 930, 1365, 230 and 696
# draws on one run, below the target the run was asked for.

#' Run each chain in its own process
#'
#' @param fit A `ctJuliaFit`, or the shell `ctFit(optimize = FALSE)` is about to
#'   fill. Carries the model and data a worker rebuilds its objective from.
#' @param target From [.ctBackendSampleTarget()]: what to sample, in a form that
#'   survives serialisation.
#' @param chains Number of chains.
#' @param warmup,draws Per chain.
#' @param cores The most processes to run at once, and the threads they share.
#' @param handles Optional warmed pool from [.ctBackendWarmWorkers()]. When
#'   absent the workers are started here and the compile is paid in full.
#' @param control,saveEffects,seed,verbose As for [ctFitUncertainty()].
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
  # `.ctFitIsJulia()`, not the class literal: the named predicate is the one
  # place this question is spelled, and `test-duplication-ratchet.R` counts
  # the spellings that are not. Brought down here because the file was open.
  if (!.ctFitIsJulia(fit)) return(NULL)

  # No more processes than `cores`: one worker would gain nothing over this
  # session but its own startup, so then the chains run here.
  workers <- .ctBackendSampleWorkers(chains, cores)
  if (workers < 2L) return(NULL)
  # Contiguous blocks, one engine call per worker; pooled back in chain order.
  blocks <- split(seq_len(chains), sort(rep_len(seq_len(workers), chains)))
  # Every core goes to some worker: a worker's own subject split uses its
  # share, which is the same nesting the in-process path does, except that
  # the chains no longer share an allocator.
  threads <- .ctBackendWorkerThreads(cores, workers)

  if (is.null(handles)) {
    handles <- .ctBackendWarmWorkers(fit, workers = workers,
      values = target$estimate, threads = threads)
    if (is.null(handles)) return(NULL)
  }
  # `resolved()` does not block, so this can say what the pause is for before
  # paying it. The pool is started before the optimisation so that the compile
  # overlaps it, and usually nothing is left to wait for -- but a fast
  # optimisation finishes first, and then the parent sits silent for the
  # remainder of a 26-43 s compile with no indication of why.
  pending <- sum(!vapply(handles, .ctBackendChainOver, logical(1)))
  if (isTRUE(progress) && pending > 0L) {
    message(pending, " of ", length(handles), " chain worker(s) still ",
      "compiling for this model shape.")
  }
  # Nothing warmed means no worker in the pool can build this model's
  # objective, so none of them can run a chain either. Returning here samples
  # in this session instead, which is where a broken pool used to arrive
  # anyway -- but by way of every chain failing first, so the fallback cost a
  # chain's startup per chain and reported itself as "2 of 2 chains failed in
  # their worker process" rather than as a pool that was not there.
  if (.ctBackendWarmWait(handles, verbose = verbose) < 1L) {
    warning("No sampling worker process could prepare this model, so the ",
      "chains ran in this session. ctJuliaKill('workers') clears the pool if it ",
      "was left behind by an interrupted run.", call. = FALSE)
    return(NULL)
  }

  # Each chain runs in its own process, so its printed output sits in that
  # process's own stdout buffer and only reaches the parent -- all at once,
  # after the fact -- when `future::value()` collects it.
  # `ctFitUncertainty(fit, uncertainty = 'sample', verbose = TRUE)` under
  # `processes = TRUE` used to print nothing at all for exactly this reason:
  # the reporting existed, in the worker, and had nowhere to go
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
  progress_files <- if (report) vapply(seq_len(workers), function(i)
    tempfile(pattern = sprintf("ctsem_sample_worker%d_", i), fileext = ".progress"),
    character(1)) else NULL
  if (report) on.exit(unlink(progress_files, force = TRUE), add = TRUE)

  # `seed + k - 1`, which makes this path reproduce the in-process one exactly.
  #
  # The engine gives chain `c` the stream `Xoshiro(seed + c)`, and continues it
  # in batch `a` on `Xoshiro(seed + 1000 a + c)`. A worker's block starting at
  # chain `k` is handed `S = seed + k - 1`, so its `j`th chain draws
  # `Xoshiro(seed + k - 1 + j)` -- the stream in-process chain `k + j - 1`
  # would have used -- and, since the parent asks every worker for the same
  # batches the session's own rule would have, its continuations match too.
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
  #
  # The first batch is the run's, sized from all its chains as the session's
  # sampler sizes it (`.ctBackendSampleBudget()`), so every chain starts the
  # same length; the parent's coordinator decides every batch after it.
  settings <- .ctBackendSampleControl(control)
  budget <- .ctBackendSampleBudget(draws, chains, settings)
  coordinator <- .ctBackendChainCoordinator(fit, settings, budget, chains,
    workers, warmup)
  if (!is.null(coordinator)) on.exit(coordinator$close(), add = TRUE)
  workercontrol <- .ctBackendWorkerControl(control)
  results <- lapply(seq_along(blocks), function(w) {
    worker_file <- if (report) progress_files[w] else NULL
    tryCatch(
      # The namespace lookup is explicit because the expression is evaluated in
      # a worker process, where only the installed ctsem exists. `ctsem:::`
      # would do the same job but draws a CRAN NOTE for ::: on our own objects.
      future::future(utils::getFromNamespace(".ctBackendSampleChainBlock",
        "ctsem")(fit, target, warmup,
        budget$first, threads[w], workercontrol, saveEffects,
        as.integer(seed) + blocks[[w]][1L] - 1L, length(blocks[[w]]),
        progress_file = worker_file,
        coordinate = if (is.null(coordinator)) NULL else
          list(dir = coordinator$dir, worker = w)),
        seed = TRUE),
      error = function(e) NULL)
  })
  if (report || !is.null(coordinator)) {
    .ctBackendReportProcesses(results, if (report) progress_files else NULL,
      chains = workers, overwrite = .ctProgressOverwrite(verbose),
      interval = if (is.null(coordinator)) 1 else 0.25,
      tick = if (is.null(coordinator)) NULL else coordinator$tick)
  }
  # One entry per worker, in chain order; a worker that failed outright
  # leaves NULL.
  drawn <- lapply(results, function(handle) if (is.null(handle)) NULL else
    tryCatch(future::value(handle), error = function(e) NULL))
  # A worker that failed returns an empty list carrying an `error` attribute,
  # not NULL, so testing for NULL alone would let it through and the failure
  # would surface later as an empty matrix in the pooling.
  ok <- vapply(drawn, function(d)
    !is.null(d) && !is.null(d$draws) && length(d$draws) > 0,
    logical(1))
  if (sum(ok) < workers) {
    failed <- sum(lengths(blocks)[!ok])
    warning(failed, " of ", chains, " chains failed in their worker ",
      "process, so sampling fell back to this session. The first error was: ",
      .ctBackendFirstError(drawn), call. = FALSE)
    return(NULL)
  }
  .ctBackendPoolChains(fit, target, drawn, chains, warmup, draws, saveEffects,
    ess_target = .ctBackendSampleTargetESS(control),
    rhat_target = .ctBackendSampleTargetRhat(control))
}

# The parent's half of the stopping rule when chains are in worker processes.
#
# `NULL` when the run has no effective-size target: every chain then draws the
# count asked for, and there is nothing to decide. Otherwise a directory the
# workers write each batch into and read each verdict from, and `tick()`, which
# the progress poll calls: once every worker has written the current batch, it
# pools them in chain order and asks `ctsem_sample_verdict` -- the rule
# `_sample_until_target` applies in one session -- how many more draws per
# chain the run needs, zero to stop. Those draws are shared out by each
# worker's speed over its last batch, so a worker with more cores, or fewer
# chains, takes more and all finish the round together; chain lengths then
# differ, and the verdict weights each chain by its length. A worker that has
# ended without writing gets a zero for everyone else, so no worker waits on a
# chain that is never coming, and its failure surfaces at collection.
# `close()` deletes the directory, which also releases any worker still
# waiting.
#' @keywords internal
.ctBackendChainCoordinator <- function(fit, settings, budget, chains, workers,
  warmup = 0L) {
  min_ess <- as.numeric(.ctJuliaOr(settings$min_ess, 0))
  mean_ess <- as.numeric(.ctJuliaOr(settings$mean_ess, 0))
  if (!isTRUE(max(min_ess, mean_ess) > 0)) return(NULL)
  dir <- tempfile("ctsem_chains_")
  dir.create(dir)
  first <- as.integer(budget$first)
  max_draws <- as.integer(.ctJuliaOr(budget$max_draws, first))
  rhat_target <- as.numeric(.ctJuliaOr(settings$rhat_target, 1.01))
  module <- .ctJuliaModule(fit$model_spec$project)
  state <- new.env()
  state$round <- 1L
  state$was_met <- FALSE
  state$done <- FALSE
  # Each worker's draws for the round, one count per worker: zero stops it.
  answer <- function(wanted) {
    wanted <- rep_len(as.integer(wanted), workers)
    for (w in seq_len(workers)) {
      path <- file.path(dir, sprintf("verdict_w%d_r%d.rds", w, state$round))
      saveRDS(wanted[w], paste0(path, ".tmp"))
      file.rename(paste0(path, ".tmp"), path)
    }
    if (all(wanted <= 0L)) state$done <- TRUE
    state$round <- state$round + 1L
  }
  tick <- function(resolved) {
    if (state$done) return(NULL)
    files <- file.path(dir, sprintf("w%d_r%d.rds", seq_len(workers), state$round))
    have <- file.exists(files)
    if (!all(have)) {
      if (any(resolved & !have)) answer(0L)
      return(NULL)
    }
    batches <- lapply(files, readRDS)
    pooled <- do.call(cbind, lapply(batches, function(b) as.matrix(b$draws)))
    nblock <- vapply(batches, function(b) as.integer(b$nchains), integer(1))
    lengths <- unlist(lapply(batches, function(b)
      rep(as.integer(b$total), as.integer(b$nchains))))
    if (sum(nblock) != chains || ncol(pooled) != sum(lengths)) {
      answer(0L)
      return(NULL)
    }
    verdict <- .ctJuliaGet(module$ctsem_sample_verdict(.ctJuliaPut(pooled),
      lengths, ndraws = first, min_ess = min_ess, mean_ess = mean_ess,
      max_draws = max_draws, rhat_target = rhat_target,
      was_met = state$was_met))
    state$was_met <- isTRUE(verdict$met)
    perchain <- if (length(unique(lengths)) > 1L) paste0(sum(lengths),
      " draws (", min(lengths), "-", max(lengths), " per chain)") else
      paste0(lengths[1L], " draws per chain")
    said <- paste0("  ", perchain, ": min ESS (bulk and tail) ",
      round(verdict$worst, 1), ", mean ESS ", round(verdict$average, 1),
      ", worst R-hat ", round(verdict$rhat, 3),
      if (isTRUE(verdict$confirmed)) " -- targets met" else
        if (isTRUE(verdict$met)) " -- targets met, confirming" else "")
    unlink(files)
    answer(.ctBackendShareDraws(as.integer(verdict$wanted) * chains,
      rates = vapply(batches, function(b) as.numeric(b$rate), numeric(1)),
      nchains = nblock))
    said
  }
  list(dir = dir, tick = tick,
    close = function() unlink(dir, recursive = TRUE, force = TRUE))
}

# The worker's half: a function the engine calls after each batch with the
# block's population draws so far and the draws per chain, which writes them
# for the parent -- with how fast this worker's chains ran, iterations per
# chain per second over the batch, warmup counted in the first -- and returns
# the parent's answer: how many more draws per chain, zero to stop. The clock
# runs only while the chains do, not while this waits. A directory that has
# gone means the parent has finished or been interrupted, and stops the chain
# rather than leaving it waiting.
#' @keywords internal
.ctBackendChainCoordinate <- function(dir, worker, warmup = 0L, interval = 0.05) {
  round <- 0L
  started <- Sys.time()
  previous <- 0L
  function(pooled, total) {
    round <<- round + 1L
    total <- as.integer(total)
    pooled <- as.matrix(pooled)
    iterations <- total - previous + if (round == 1L) as.integer(warmup) else 0L
    secs <- max(1e-3, as.numeric(difftime(Sys.time(), started, units = "secs")))
    previous <<- total
    path <- file.path(dir, sprintf("w%d_r%d.rds", worker, round))
    ok <- tryCatch({
      suppressWarnings(saveRDS(list(draws = pooled, total = total,
        nchains = ncol(pooled) %/% max(1L, total), rate = iterations / secs),
        paste0(path, ".tmp")))
      file.rename(paste0(path, ".tmp"), path)
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) return(0L)
    verdict <- file.path(dir, sprintf("verdict_w%d_r%d.rds", worker, round))
    while (!file.exists(verdict)) {
      if (!dir.exists(dir)) return(0L)
      Sys.sleep(interval)
    }
    started <<- Sys.time()
    as.integer(tryCatch(readRDS(verdict), error = function(e) 0L))
  }
}

# The draws a round asks for, `total` over every chain, shared out among the
# workers so that each takes as long as the others: a worker's chains get
# draws in proportion to the rate they ran at (iterations per chain per
# second), so a worker with twice the cores per chain draws about twice as
# many. None gets fewer than one per chain while the run goes on, because
# zero is the signal to stop.
#' @keywords internal
.ctBackendShareDraws <- function(total, rates, nchains) {
  total <- as.integer(total)[1L]
  if (is.na(total) || total <= 0L) return(integer(length(rates)))
  rates <- as.numeric(rates)
  rates[!is.finite(rates) | rates <= 0] <- NA_real_
  if (all(is.na(rates))) rates[] <- 1
  rates[is.na(rates)] <- min(rates, na.rm = TRUE)
  pmax(1L, as.integer(round(total * rates / sum(nchains * rates))))
}

# The control list a worker's chains sample by: the run's, with no stopping
# target of its own.
#
# A worker judges nothing. The parent decides every batch after the first over
# all the chains (`.ctBackendChainCoordinator()`), and the first batch arrives
# as the worker's draw count, so the worker's own budget must leave that count
# as it is -- which a target or a `maxDraws` here would not. Workers that
# stopped on their own share of the target are what this replaced: their
# chains came back at different lengths and pooled cut to the shortest.
#' @keywords internal
.ctBackendWorkerControl <- function(control) {
  control$minESS <- 0
  control$meanESS <- NULL
  control$maxDraws <- NULL
  control
}

# A worker's block of chains, as one engine call.
#
# Deliberately the shared runner's own engine call: a second sampler for the
# process path would be a second thing to keep correct, and it would have to be
# told the same things anyway. What the worker does not run is the assembly --
# constraining a block's draws only to throw them away when the pool is
# assembled is work nobody reads. One call rather than one per chain because
# the parent extends every chain together, so a block's chains have to be
# running at once, not in turn; `seed` is the block's, from which the engine
# gives its `j`th chain the stream that chain has in one session.
#
# `coordinate`, when given, is where the parent decides when the chains stop
# (`.ctBackendChainCoordinate()`).
#
# It used to call `ctFitUncertainty(fit, uncertainty = 'sample')`, which fixed
# the target as well as the code: the joint entry, and a refusal for any fit
# without a Laplace spec. That is right for its own callers and wrong for two
# of the three routes `ctFit(optimize = FALSE)` can take, which is why the
# target now travels explicitly.
#
# `progress_file`, when given, replaces `control$callback` for this call: a
# user's own callback is an R closure over the parent session (a plot device,
# a Shiny reactive) and calling it from here would try to reach across a
# process boundary that does not carry it. What can cross is a path, and the
# parent polls what gets written there -- see `.ctBackendReportProcesses`.
#' @keywords internal
.ctBackendSampleChainBlock <- function(fit, target, warmup, draws, threads,
  control, saveEffects, seed, nchains = 1L, progress_file = NULL,
  coordinate = NULL) {
  tryCatch({
    if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = threads)
    # The pool's sessions were started at the workers' shares of `cores`, and
    # which worker takes which block is the pool's choice, so a block runs at
    # the width of the session it lands in; `threads` started one only if the
    # pool had not.
    width <- suppressWarnings(as.integer(.ctBackendSessionThreads())[1L])
    if (is.na(width) || width < 1L) width <- as.integer(threads)
    callback <- if (is.null(progress_file)) NULL else
      .ctBackendProgressFileWriter(progress_file)
    result <- suppressWarnings(suppressMessages(
      .ctBackendSampleEngine(fit, target, chains = as.integer(nchains),
        warmup = warmup, draws = draws, cores = width,
        saveEffects = saveEffects, seed = seed, control = control,
        verbose = FALSE, progress = FALSE, callback = callback,
        coordinate = if (is.null(coordinate)) NULL else
          .ctBackendChainCoordinate(coordinate$dir, coordinate$worker,
            warmup = warmup))))
    .ctBackendChainResult(result)
  }, error = function(e) structure(list(), error = conditionMessage(e)))
}

# Each worker's share of `cores`: as even as whole threads allow, the
# remainder one each to the first workers -- 6 over 4 is 2, 2, 1, 1 -- so
# every core goes to some worker. Chains on a wider worker run faster and the
# parent gives them more draws (`.ctBackendShareDraws()`).
#' @keywords internal
.ctBackendWorkerThreads <- function(cores, workers) {
  cores <- suppressWarnings(as.integer(cores)[1L])
  workers <- suppressWarnings(as.integer(workers)[1L])
  if (is.na(workers) || workers < 1L) workers <- 1L
  if (is.na(cores) || cores < workers) return(rep(1L, workers))
  base <- cores %/% workers
  base + as.integer(seq_len(workers) <= cores - base * workers)
}

# How many worker processes a run uses: one per chain, at most `cores`.
#' @keywords internal
.ctBackendSampleWorkers <- function(chains, cores) {
  cores <- suppressWarnings(as.integer(cores)[1L])
  if (is.na(cores) || cores < 1L) cores <- 1L
  min(as.integer(chains), cores)
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

#' One line describing every chain
#'
#' The whole report on one line, so that it can be overwritten in place: a
#' carriage return goes to the start of the last visual row, so a block of
#' `chains` lines cannot be updated without cursor movement that not every
#' console honours. One line can, and the same line is what the
#' single-process path prints.
#'
#' What is on it, and why in this order. The phase and the per-chain iteration
#' counts answer "is this going to finish"; the counts share their total and
#' their phase label whenever the chains agree on them, which is nearly
#' always and is what keeps the line inside a terminal width. Then the time
#' remaining, from the *slowest* chain, because that is when the run ends.
#' Then the log posterior per chain, which is what says whether the draws are
#' going anywhere sensible -- and it is per chain rather than summarised
#' because a single chain stuck in a bad region is the failure this reveals,
#' and an average hides exactly that. Divergences last, and only once there
#' are some.
#'
#' @param infos One entry per chain from [.ctBackendReadProgressFile()], with
#'   `NULL` for a chain that has not written yet.
#' @param chains Number of chains, so that the line can say when it is
#'   describing fewer of them than are running.
#' @param elapsed Seconds each chain has spent in its current phase, for the
#'   rate the estimate extrapolates from.
#' @param width Console width to truncate to. A line that wraps cannot be
#'   overwritten in place -- the earlier rows are left behind as debris.
#' @return One string, or `NULL` when no chain has reported yet.
#' @keywords internal
.ctBackendProcessLine <- function(infos, chains, elapsed,
  width = getOption("width", 80L)) {
  have <- which(!vapply(infos, is.null, logical(1)))
  if (!length(have)) return(NULL)
  phase <- vapply(infos[have], function(i) i$phase, character(1))
  iteration <- vapply(infos[have], function(i) i$iteration, integer(1))
  total <- vapply(infos[have], function(i) i$total, integer(1))
  logp <- vapply(infos[have], function(i) i$logp, numeric(1))
  divergent <- vapply(infos[have], function(i) i$divergent, integer(1))

  # Chains in the same phase with the same target are described once: the
  # alternative repeats `warmup` and `/500` per chain, and four of those do
  # not fit on a line that has to hold four log posteriors as well.
  groups <- unique(paste(phase, total))
  counts <- vapply(groups, function(g) {
    use <- paste(phase, total) == g
    sprintf("%s %s/%d", phase[use][1L],
      paste(iteration[use], collapse = ", "), total[use][1L])
  }, character(1), USE.NAMES = FALSE)

  # The slowest chain's estimate, not the mean of them: the pooled draws are
  # not there until every chain has finished.
  rate <- ifelse(iteration > 0 & elapsed[have] > 0, iteration / elapsed[have], 0)
  remaining <- ifelse(rate > 0 & total > iteration, (total - iteration) / rate, 0)

  # Fit the line to the console by dropping what matters least, rather than by
  # cutting the end off. Truncation removed the log posteriors -- four chains
  # of a 500-draw run do not fit in eighty columns with everything on -- and
  # those are the whole reason this line is per chain, so they are the last
  # thing to go and the phrase around the time estimate is the first.
  compose <- function(digits, phrase, showdiv) {
    parts <- counts
    if (any(remaining > 0)) {
      parts <- c(parts, paste0(.ctDuration(max(remaining)),
        if (phrase) " at this rate" else ""))
    }
    parts <- c(parts, paste("logp", paste(
      sprintf(paste0("%.", digits, "g"), logp), collapse = ", ")))
    if (showdiv && any(divergent > 0L)) {
      parts <- c(parts, paste("div", paste(divergent, collapse = ", ")))
    }
    # Every field above is a positional list over the chains that have
    # reported, so a chain missing from them shifts the rest silently. Said
    # only when one is missing, which is the first second of a run -- and a
    # worker that started and never reported, where it is the whole story.
    if (length(have) < chains) {
      parts <- c(parts, sprintf("%d of %d chains reporting", length(have),
        chains))
    }
    paste0("  ", paste(parts, collapse = " | "))
  }
  variants <- list(compose(6, TRUE, TRUE), compose(6, FALSE, TRUE),
    compose(4, FALSE, TRUE), compose(4, FALSE, FALSE))
  width <- as.integer(width)[1L]
  if (is.na(width) || width <= 10L) return(variants[[1L]])
  for (line in variants) if (nchar(line) <= width - 1L) return(line)
  substr(variants[[length(variants)]], 1L, width - 1L)
}

#' Report per-chain progress from the parent while chain processes run
#'
#' Polls each chain's progress file on a short interval and reports every
#' chain on one line -- the shape asked for when chains are processes: each
#' worker's own printed progress sits in output that never reaches the parent
#' until the chain is already done, so the parent reports instead, from what
#' the workers wrote rather than from what they printed.
#'
#' Emitted through [.ctProgressSink()], which is what the engine's own
#' progress line is delivered by -- so the two look alike wherever a console
#' distinguishes a message from plain output, and the rules about carriage
#' returns and padding live in one place rather than in two copies that drift.
#'
#' Overwritten in place where a carriage return means something, exactly as
#' `CTSEMProgress` does it in `progress.jl` and under the same detection --
#' `.ctProgressOverwrite()`, so a log file, a knitr chunk or a captured stream
#' gets its updates on their own lines and rarely instead. Printing a fresh
#' block every poll is what this replaced: a five-minute sample left hundreds
#' of lines of scrollback saying nothing the last of them did not.
#'
#' A chain that has finished keeps its last reported state on the line rather
#' than dropping out of it. The line would otherwise shorten as chains end,
#' which reads as chains disappearing, and the padding that makes in-place
#' overwriting work would leave the tail of the longer line behind.
#'
#' @param results Future handles from [.ctBackendSampleProcesses()], one per
#'   chain, possibly containing `NULL` for a chain that never started.
#' @param progress_files One path per chain, written by
#'   [.ctBackendProgressFileWriter()], or `NULL` to print nothing and only
#'   poll for `tick`.
#' @param chains Number of chains.
#' @param interval Seconds between polls.
#' @param overwrite Update one line in place rather than printing each report
#'   on its own line.
#' @param tick Called with which workers have finished on every poll: the
#'   parent's half of the stopping rule (`.ctBackendChainCoordinator()`). Text
#'   it returns is printed on its own line.
#' @return `NULL`, invisibly. Called for its printing.
#' @keywords internal
.ctBackendReportProcesses <- function(results, progress_files, chains,
  interval = 1, overwrite = .ctProgressOverwrite(1), tick = NULL) {
  now <- Sys.time()
  phase_started <- rep(now, chains)
  phase_seen <- rep(NA_character_, chains)
  # Last-known state per chain, so a finished chain stays on the line and a
  # poll that catches a file mid-write does not blank one.
  infos <- vector("list", chains)
  emit <- .ctProgressSink(overwrite)
  shown <- NULL
  repeat {
    resolved <- vapply(results, .ctBackendChainOver, logical(1))
    if (!is.null(tick)) {
      said <- tick(resolved)
      if (!is.null(said) && !is.null(progress_files)) {
        emit("", "break")
        message(said)
        shown <- NULL
      }
    }
    if (is.null(progress_files)) {
      if (all(resolved)) break
      Sys.sleep(interval)
      next
    }
    for (k in seq_len(chains)) {
      info <- .ctBackendReadProgressFile(progress_files[k])
      if (is.null(info)) next
      if (is.na(phase_seen[k]) || !identical(phase_seen[k], info$phase)) {
        phase_started[k] <- Sys.time()
        phase_seen[k] <- info$phase
      }
      infos[[k]] <- info
    }
    elapsed <- as.numeric(difftime(Sys.time(), phase_started, units = "secs"))
    line <- .ctBackendProcessLine(infos, chains, elapsed)
    # An update that says exactly what is already on the line is not written.
    # Every chain has finished by the last poll, so the closing state would
    # otherwise be written twice -- invisible in a console, which overwrites
    # itself, and a duplicated line everywhere a carriage return is a character.
    if (!is.null(line) && !identical(line, shown)) {
      emit(line, "update")
      utils::flush.console()
      shown <- line
    }
    if (all(resolved)) break
    Sys.sleep(interval)
  }
  # End the line so whatever prints next -- the pooling, a diagnostic warning
  # -- starts on its own rather than inside this one. A no-op when nothing was
  # emitted, or when each update already ended its own line.
  emit("", "break")
  utils::flush.console()
  invisible(NULL)
}

# Whether a worker's future has finished, one way or the other.
#
# An error from `resolved()` counts as finished. It is one way `future` reports
# a worker process that has gone, and the `value()` pass after the poll turns
# that into a failed chain and a fallback to this session. Uncaught, it
# abandoned the whole sample instead; read as "not yet", it would have kept the
# poll waiting on a process that no longer exists.
#' @keywords internal
.ctBackendChainOver <- function(handle) {
  is.null(handle) ||
    isTRUE(tryCatch(future::resolved(handle), error = function(e) TRUE))
}

# What a worker's block of chains sends home.
#
# The engine's own result, minus the two things the pool recomputes -- R-hat and
# effective sample size are properties of the whole run and cannot be averaged
# from per-chain values -- and minus anything that would be a second copy of the
# model. Returning the whole fit would send the data and the model back across
# for every chain, having already sent them out.
#' @keywords internal
.ctBackendChainResult <- function(result) {
  ndraws <- as.integer(result$ndraws)
  nchains <- as.integer(.ctJuliaOr(result$nchains, 1L))[1L]
  list(
    # `kept x (nchains * ndraws)`, chain-major, which is the engine's own
    # layout, so pooling the blocks is a `cbind` and the assembler reshapes the
    # pool exactly as it reshapes a single-process result.
    draws = matrix(as.numeric(result$draws), ncol = nchains * ndraws),
    npar = as.integer(result$npar), ndim = as.integer(result$ndim),
    ndraws = ndraws, nchains = nchains,
    ndivergent = as.integer(result$ndivergent),
    forward_gradients = as.integer(.ctJuliaOr(result$forward_gradients, 0L))[1L],
    explosive_passes = as.integer(.ctJuliaOr(result$explosive_passes, 0L))[1L],
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
    effect_sd = as.numeric(result$effect_sd),
    effect_cov = as.numeric(result$effect_cov),
    effect_cov_rows = as.integer(result$effect_cov_rows),
    effect_cov_cols = as.integer(result$effect_cov_cols),
    sampler = if (is.null(result$sampler)) "nuts" else as.character(result$sampler),
    # Each worker placed its own chain; the first worker's placement stands for
    # the run's in the assembled fit.
    placement = result$placement,
    scale_accept = as.numeric(result$scale_accept),
    ncp_accept = as.numeric(result$ncp_accept),
    # The chain's coreset, where it sampled on one (sample_surrogate.jl).
    coreset = isTRUE(result$coreset), nunits = result$nunits,
    coreset_size = as.integer(result$coreset_size),
    coreset_proposed = as.integer(result$coreset_proposed),
    coreset_accepted = as.integer(result$coreset_accepted),
    coreset_abandoned = as.logical(result$coreset_abandoned))
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
  saveEffects = FALSE, ess_target = NA_real_, rhat_target = 1.01) {
  mats <- lapply(drawn, function(d) as.matrix(d$draws))
  npar <- as.integer(target$npar)
  kept <- nrow(mats[[1]])
  if (!all(vapply(mats, nrow, integer(1)) == kept) || kept < npar) return(NULL)

  # Every chain is kept whole. A worker's chains share a length; workers'
  # need not, since the parent extends each in proportion to its speed
  # (`.ctBackendChainCoordinator()`), and the diagnostics below weight each
  # chain by its length. Chains that stopped separately used to be cut to the
  # shortest here, which threw most of a run away; a block whose draws do not
  # add up is a fault, and the run is sampled again in this session.
  nper <- vapply(drawn, function(d) as.integer(d$ndraws), integer(1))
  nblock <- vapply(drawn, function(d) as.integer(.ctJuliaOr(d$nchains, 1L))[1L],
    integer(1))
  lengths <- rep(nper, nblock)
  if (any(nper < 1L) || sum(nblock) != chains ||
      any(vapply(mats, ncol, integer(1)) != nblock * nper)) return(NULL)
  ndraws <- if (length(unique(lengths)) > 1L) lengths else nper[1L]
  perdraw <- function(field) unlist(lapply(drawn, function(d) as.numeric(d[[field]])))

  # `cbind` puts chain 1's draws first, then chain 2's, which is the chain-major
  # layout `ctsem_sample_diagnostics` indexes, given the chain lengths. Each
  # block is chain-major already and the blocks are in chain order.
  # Getting this wrong would not error -- it would silently mix the chains and
  # report R-hat over the mixture, which is always reassuring.
  pooled <- do.call(cbind, mats)
  ndim <- as.integer(drawn[[1]]$ndim)
  if (!isTRUE(is.finite(ndim))) ndim <- kept
  keepeffects <- isTRUE(saveEffects) && kept > npar
  # A chain that sent effect draws nobody asked for: the assembler reshapes to
  # `npar` rows in that case, so the extra rows have to go rather than be
  # interleaved into plausible-looking nonsense.
  if (!keepeffects && kept > npar) pooled <- pooled[seq_len(npar), , drop = FALSE]

  # The effect summaries recombine rather than concatenate. Each block reports
  # the mean and sd over its own `n_b` draws; the pooled mean weights the
  # blocks' means by `n_b`, and the pooled variance is the within-block sum of
  # squares plus the between-block one, over the pooled degrees of freedom.
  # Averaging the blocks' standard deviations instead would understate the
  # spread by exactly the part between chains -- the part R-hat is about.
  means <- do.call(rbind, lapply(drawn, function(d) as.numeric(d$effect_mean)))
  sds <- do.call(rbind, lapply(drawn, function(d) as.numeric(d$effect_sd)))
  effectmean <- numeric(0)
  effectsd <- numeric(0)
  if (!is.null(means) && ncol(means) > 0L && identical(dim(means), dim(sds))) {
    nb <- nblock * nper
    effectmean <- colSums(means * nb) / sum(nb)
    total <- colSums((nb - 1) * sds^2) +
      colSums(nb * sweep(means, 2, effectmean)^2)
    effectsd <- sqrt(total / max(1, sum(nb) - 1))
  }
  # Each group's covariance pools the same way, entry by packed entry: the
  # blocks' own products plus the products of their means' departures.
  covs <- do.call(rbind, lapply(drawn, function(d) as.numeric(d$effect_cov)))
  rows <- as.integer(drawn[[1]]$effect_cov_rows)
  cols <- as.integer(drawn[[1]]$effect_cov_cols)
  effectcov <- numeric(0)
  if (length(effectmean) && !is.null(covs) && ncol(covs) > 0L &&
      ncol(covs) == length(rows)) {
    nb <- nblock * nper
    departure <- sweep(means, 2, effectmean)
    effectcov <- (colSums((nb - 1) * covs) + colSums(nb * departure[, rows, drop = FALSE] *
      departure[, cols, drop = FALSE])) / max(1, sum(nb) - 1)
  }

  module <- .ctJuliaModule(fit$model_spec$project)
  # The population block alone: R-hat over every random effect as well would
  # cost more than it says, and the assembler reads only the first `npar`.
  # Not `diag`, which would shadow `base::diag` for the rest of the function.
  pooldiag <- .ctJuliaGet(module$ctsem_sample_diagnostics(
    .ctJuliaPut(pooled[seq_len(npar), , drop = FALSE]), lengths))

  result <- list(
    draws = pooled, npar = npar, ndim = ndim, ndraws = ndraws,
    rhat = as.numeric(pooldiag$rhat), ess = as.numeric(pooldiag$ess),
    ess_tail = as.numeric(pooldiag$ess_tail),
    # Summed across chains, because they count events; the step size and E-BFMI
    # are per chain and stay per chain.
    ndivergent = sum(vapply(drawn, function(d) as.integer(d$ndivergent), integer(1))),
    forward_gradients = sum(vapply(drawn,
      function(d) as.integer(.ctJuliaOr(d$forward_gradients, 0L))[1L], integer(1))),
    explosive_passes = sum(vapply(drawn,
      function(d) as.integer(.ctJuliaOr(d$explosive_passes, 0L))[1L], integer(1))),
    warmup_divergent = sum(unlist(lapply(drawn,
      function(d) as.integer(d$warmup_divergent)))),
    nsaturated = sum(vapply(drawn, function(d) as.integer(d$nsaturated), integer(1))),
    max_depth = max(vapply(drawn, function(d) as.integer(d$max_depth), integer(1))),
    stepsize = unlist(lapply(drawn, function(d) as.numeric(d$stepsize))),
    ebfmi = unlist(lapply(drawn, function(d) as.numeric(d$ebfmi))),
    accept = perdraw("accept"), depth = perdraw("depth"),
    energy = perdraw("energy"),
    effect_mean = effectmean, effect_sd = effectsd,
    effect_cov = effectcov, effect_cov_rows = rows, effect_cov_cols = cols,
    # Per chain, like the step size.
    sampler = drawn[[1]]$sampler,
    placement = drawn[[1]]$placement,
    scale_accept = unlist(lapply(drawn, function(d) as.numeric(d$scale_accept))),
    ncp_accept = unlist(lapply(drawn, function(d) as.numeric(d$ncp_accept))),
    coreset = any(vapply(drawn, function(d) isTRUE(d$coreset), logical(1))),
    nunits = drawn[[1]]$nunits,
    coreset_size = unlist(lapply(drawn, function(d) d$coreset_size)),
    coreset_proposed = unlist(lapply(drawn, function(d) d$coreset_proposed)),
    coreset_accepted = unlist(lapply(drawn, function(d) d$coreset_accepted)),
    coreset_abandoned = unlist(lapply(drawn, function(d) d$coreset_abandoned)))

  out <- .ctBackendSampleAssemble(fit, result, npar, keepeffects,
    as.integer(chains), warmup, ndraws, target$hessian,
    target$estimate[seq_len(npar)], marginal = isTRUE(target$marginal),
    ess_target = ess_target, rhat_target = rhat_target)
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
