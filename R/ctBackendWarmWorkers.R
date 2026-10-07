# Starting the worker processes before they are needed.
#
# A worker that is going to run a chain has to reach the state the parent is
# already in: R started, ctsem loaded, Julia running, the objective built, and
# -- the expensive part -- the engine specialised for this model's dimensions.
# Measured on this machine that is roughly 3.6 s of R, 6.3 s of `ctJuliaSetup`,
# 1.8 s of building, and 26-43 s of compiling, and the compiling is unavoidable
# per process: it is Julia specialising its matrix kernels on the latent and
# manifest dimensions, which are carried in the workspace type, and there is no
# supported way to move compiled code between processes.
# `inst/julia/ContinuousTimeSEM/src/workspace_buffers.jl` records what was
# measured against that cost and why it is one worth paying rather than a defect.
#
# What can be avoided is paying it *serially*. A fit optimises before it
# samples, and on the models where process-parallel chains are worth having that
# optimisation runs for minutes. A worker needs nothing from it -- only the
# model, which exists before optimisation starts -- so it can compile during it
# and be waiting when sampling begins.
#
# Nothing here is required for a fit to work. Every failure path returns `NULL`
# or `0`, and the caller samples with threads instead.

# Holds what this file changes globally, so it can be put back.
.ct_warm_state <- new.env(parent = emptyenv())

# `future` is in Suggests rather than Imports, so every entry point checks for
# it. Chains as processes is worth having and not worth making the package fail
# to install without.
#' @keywords internal
.ctBackendCanWarm <- function() {
  requireNamespace("future", quietly = TRUE)
}

#' Start worker sessions and get them ready to sample
#'
#' Launches `workers` background R sessions and sets each one evaluating one
#' gradient for `object`, which is what forces Julia to build the objective and
#' compile the engine for this model's shape. Returns immediately: the sessions
#' work while the caller carries on.
#'
#' The sessions persist. `future::multisession` keeps its pool across calls, so
#' the warmed sessions are the same ones a later chain runs in -- which is the
#' point, and the reason this sets a plan rather than making a one-shot cluster.
#'
#' @param object A prepared `ctJuliaModel` or a `ctJuliaFit`.
#' @param workers How many sessions to start.
#' @param values Parameter vector to evaluate at. Any finite point compiles the
#'   same code; a fit's own estimate is the natural choice when there is one.
#' @param threads Julia threads each worker's session starts with. Fixed when
#'   the session starts, so a pool warmed at another width is replaced rather
#'   than reused.
#' @return A list of future handles, or `NULL` when warming is unavailable.
#' @keywords internal
.ctBackendWarmWorkers <- function(object, workers, values = NULL, threads = 1L) {
  threads <- max(1L, suppressWarnings(as.integer(threads)[1L]), na.rm = TRUE)
  workers <- suppressWarnings(as.integer(workers))
  if (!.ctBackendCanWarm() || !isTRUE(workers >= 1L)) return(NULL)
  if (!inherits(object, c("ctJuliaModel", "ctJuliaFit"))) return(NULL)
  npar <- .ctBackendNpar(object)
  if (!isTRUE(npar >= 1L)) return(NULL)
  if (is.null(values)) values <- rep(0, npar)
  values <- as.numeric(values)
  values <- rep_len(values, npar)
  values[!is.finite(values)] <- 0

  # Every worker must be idle before this asks for one, and if any is not then
  # the pool is left over from a run that was interrupted: R is single threaded,
  # so at the moment a fit starts warming there is nothing of ours that could
  # legitimately still be running in it.
  #
  # That state is not merely untidy. `future::future()` waits for a free worker
  # rather than failing, so warming would queue behind a chain from the
  # abandoned run -- and a worker the parent was mid-read of when the interrupt
  # landed answers with the previous call's reply, the same desynchronisation
  # the Julia bridge suffers. `.ctJuliaInterruptSafe` clears the pool when it
  # sees the interrupt; this is the same repair for every other way it can
  # happen, including an interrupt that reached R somewhere the handler did not.
  #
  # Re-planning the identical plan is a no-op in future -- it keeps the pool,
  # which is the whole point of warming across fits -- so the reset has to be
  # explicit.
  if (.ctBackendWarmPoolStale(workers)) {
    message("Sampling workers from an interrupted run are still busy, so they ",
      "are being replaced. This model's shape recompiles in the new ones.")
    .ctBackendWarmStop(NULL)
  }
  # A Julia session's thread count is fixed when it starts, so a pool whose
  # sessions started at another width cannot serve this call: replaced, or a
  # chain asking for four threads runs on the one its worker was warmed with.
  # That was every processes-mode chain: workers warmed at one thread, so 16
  # cores and 4 chains used 4 (bigre, 16 cores: 352 s against 106 s for the
  # same chains as threads in one session, dev1).
  if (!is.null(.ct_warm_state$threads) && !identical(.ct_warm_state$threads, threads)) {
    .ctBackendWarmStop(NULL)
  }
  .ct_warm_state$threads <- threads

  # Workers start with a one-thread OpenBLAS, unless the caller chose one: an R
  # linked to OpenBLAS otherwise starts a thread per CPU in every worker, and
  # those spin while it loads. A worker's Julia grows its own pool to the
  # worker's share of `cores` (`ctsem_set_max_chunks!`).
  blas <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  if (is.na(blas)) {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    on.exit(Sys.unsetenv("OPENBLAS_NUM_THREADS"), add = TRUE)
  }
  started <- tryCatch({
    future::plan(future::multisession, workers = workers)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(started)) return(NULL)

  # From here on a pool exists. If launching every worker below fails, or
  # anything throws unexpectedly before this function returns, release it
  # rather than leaving a `multisession` plan and an overridden
  # connections-misuse option in place for a pool that never warmed anything --
  # that would be exactly the undocumented state this function otherwise
  # avoids. `ok` is set only once warming actually produced a usable worker; a
  # *successful* warm is deliberately never torn down here, because the pool is
  # meant to persist across fits (see the file header and `.ctBackendWarmStop`).
  ok <- FALSE
  on.exit(if (!ok) .ctBackendWarmStop(NULL), add = TRUE)

  # `future` warns when an expression leaves a connection open, because that is
  # usually a leak. Here it is the entire point: the worker opens a connection
  # to its Julia process and must keep it, so that the session which paid for
  # the compilation is the one that later runs a chain. Restored by
  # `.ctBackendWarmStop`.
  .ct_warm_state$connections <- getOption("future.connections.onMisuse")
  options(future.connections.onMisuse = "ignore")

  # All the workers are warmed identically, so nothing needs pinning: whichever
  # one is free can take whichever chain later.
  handles <- lapply(seq_len(workers), function(k) {
    tryCatch(
      # The namespace lookup is explicit because the expression is evaluated
      # in a worker process, where only the installed ctsem exists. `ctsem:::`
      # would do the same job but draws a CRAN NOTE for ::: on our own objects.
      future::future(utils::getFromNamespace(".ctBackendWarmSession",
        "ctsem")(object, values, threads), seed = TRUE),
      error = function(e) NULL)
  })
  if (!length(handles) || all(vapply(handles, is.null, logical(1)))) return(NULL)
  ok <- TRUE

  # The pool outlives this call by design (see above), so nothing in this file
  # stops it after an ordinary fit. What does: `ctJuliaWorkersStop()`, for a
  # caller that wants the memory back sooner, and this finalizer -- registered
  # once, on the environment this file already uses to remember what it
  # changed -- which releases it when the R session ends, so a long-lived host
  # process (a ctsemGUI/Shiny server, an interactive session left open) is not
  # the only thing standing between a warmed pool and never being cleaned up.
  if (!isTRUE(.ct_warm_state$finalizer_registered)) {
    reg.finalizer(.ct_warm_state,
      function(e) tryCatch(.ctBackendWarmStop(NULL), error = function(e) NULL),
      onexit = TRUE)
    .ct_warm_state$finalizer_registered <- TRUE
  }

  attr(handles, "started") <- Sys.time()
  handles
}

# Is the current pool unusable for a run that is about to start?
#
# True when a `multisession` plan is in force and any of its workers is busy:
# see the caller for why that can only be an orphan. Deliberately narrow -- it
# says nothing about a sequential plan, or about a pool of a different size,
# both of which the caller is about to replace anyway.
#
# Every question here is asked through `tryCatch`, because the answer only
# decides whether to spend a few seconds rebuilding: a `future` version whose
# `nbrOfFreeWorkers()` behaves differently, or a plan this cannot read, should
# leave warming to proceed exactly as it did before this existed.
#' @keywords internal
.ctBackendWarmPoolStale <- function(workers) {
  if (!.ctBackendCanWarm()) return(FALSE)
  current <- tryCatch(future::plan("list")[[1L]], error = function(e) NULL)
  if (is.null(current) || !inherits(current, "multisession")) return(FALSE)
  free <- tryCatch(future::nbrOfFreeWorkers(), error = function(e) NA_integer_)
  total <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
  if (!isTRUE(is.finite(free)) || !isTRUE(is.finite(total))) return(FALSE)
  free < total
}

# Runs inside a worker. Wrapped, because a worker that cannot warm should leave
# the parent to sample with threads rather than take the fit down with it.
#' @keywords internal
.ctBackendWarmSession <- function(object, values, threads = 1L) {
  tryCatch({
    # Only when this process has no Julia yet. Asking again warns that the
    # thread count cannot be changed on a running session, which is true and
    # not worth saying: a warmed worker being asked for more work is the
    # expected case, not a problem.
    # Suppressed, as `.ctBackendSampleChainBlock` suppresses its own: `future`
    # relays a worker's messages to the parent when the parent collects the
    # worker's value, which here is in the middle of the parent's own run and
    # long after they were true. One warmed worker per chain put "Starting
    # Julia ..." and "Compiling the julia engine for this model shape" on the
    # parent's console once per worker, directly after the parent's Hessian --
    # which reads as the parent restarting Julia, and the last of them landed
    # inside the progress line. `.ctBackendWarmWait()` says how many warmed,
    # which is the parent's business; how each one got there is not.
    suppressMessages({
      if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = threads)
      # One evaluation does both jobs: `ctJuliaEvaluate` builds the objective,
      # which marshals the data, and evaluating it forces the specialisation.
      # The objective cache is per process and keyed on content, and this
      # process has an empty one, so nothing here is shared with the parent.
      invisible(ctsem::ctJuliaEvaluate(object, values, gradient = TRUE,
        cores = 1L))
    })
    # A list rather than a bare TRUE, so the parent can check that this worker
    # is running the same build it is. `.ctBackendWarmWait()` accepts both
    # shapes: a worker old enough to predate this returns TRUE, and that
    # absence is itself a difference worth reporting.
    #
    # Not `structure(TRUE, fingerprint = ...)`, which would look like success
    # and not be: `isTRUE()` is FALSE for a TRUE carrying attributes, so the
    # warmed count would silently drop to zero and every run would fall back to
    # in-process sampling.
    list(ok = TRUE, fingerprint = .ctBackendWarmFingerprint())
  }, error = function(e) structure(FALSE, message = conditionMessage(e)))
}

# What identifies the build a process is running, for this path.
#
# The package version and the engine's content hash.
# `.ctJuliaEngineVersion()` hashes `system.file("julia", ...)`, which resolves
# to the working tree under `devtools::load_all()` and to the installed copy
# inside a worker -- so these strings differ exactly when the two processes
# would compute different numbers. The version alone would not: a tree and its
# installed sibling usually carry the same version, which is the whole
# difficulty.
#
# Tolerant of an engine it cannot hash. This decides whether to warn, and a
# fingerprint that could not be taken is reported as unknown rather than
# allowed to take a sampling run down.
#' @keywords internal
.ctBackendWarmFingerprint <- function() {
  version <- tryCatch(as.character(utils::packageVersion("ctsem")),
    error = function(e) "?")
  engine <- tryCatch(.ctJuliaEngineVersion(), error = function(e) "?")
  paste0(version, "/", engine)
}

#' Wait for warmed workers
#'
#' @param handles From [.ctBackendWarmWorkers()].
#' @param verbose Report how many are ready.
#' @return Number of workers that warmed successfully.
#' @keywords internal
.ctBackendWarmWait <- function(handles, verbose = FALSE) {
  if (is.null(handles) || !length(handles)) return(0L)
  started <- attr(handles, "started")
  results <- lapply(handles, function(h) {
    if (is.null(h)) return(FALSE)
    tryCatch(future::value(h), error = function(e) FALSE)
  })
  warmed <- vapply(results, function(r)
    isTRUE(r) || (is.list(r) && isTRUE(r$ok)), logical(1))
  ok <- sum(warmed)
  .ctBackendWarmCheckBuild(results[warmed])
  if (isTRUE(verbose)) {
    elapsed <- if (is.null(started)) NA_real_ else
      as.numeric(difftime(Sys.time(), started, units = "secs"))
    message(ok, " of ", length(handles), " worker(s) ready",
      if (is.finite(elapsed)) paste0(" (", round(elapsed), "s since they started, ",
        "most of it alongside the fit)") else "", ".")
  }
  as.integer(ok)
}

# Warn when a warmed worker is not running this build.
#
# The failure this exists for is silent and expensive: the worker samples a
# different model from the one the parent optimised, the pooled draws come back
# plausible, and `sample$processes` says TRUE. Measured once here at 7.5e-05
# against a documented 6.2e-10, which reads as a pooling fault and is not one.
#
# The same shape as `.ctClusterCheckBuild()` on the stan path, and it is here
# for the same reason: `devtools::load_all()` puts the tree on no library path,
# so a worker's `library(ctsem)` finds the *installed* package. The fallback
# the sampling path relies on -- no installed ctsem, no warm -- does not cover
# an installed ctsem that merely differs.
#
# A worker with no fingerprint at all is reported too. That means a build
# predating this function, which is itself a difference.
#' @keywords internal
.ctBackendWarmCheckBuild <- function(results) {
  if (!length(results)) return(invisible(NULL))
  master <- tryCatch(.ctBackendWarmFingerprint(), error = function(e) NULL)
  if (is.null(master)) return(invisible(NULL))
  workers <- vapply(results, function(r)
    if (is.list(r) && is.character(r$fingerprint) &&
      length(r$fingerprint) == 1L) r$fingerprint else NA_character_,
    character(1))
  bad <- unique(workers[is.na(workers) | workers != master])
  if (!length(bad)) return(invisible(NULL))
  bad[is.na(bad)] <- "unknown (a build predating this check)"
  warning("Sampling worker process(es) loaded a different ctsem build (",
    paste(bad, collapse = ", "), " vs ", master,
    "). Their draws are from a different model than this session optimised. ",
    "Install this tree to a library and put it on .libPaths(), or sample ",
    "with processes = FALSE.", call. = FALSE)
  invisible(NULL)
}

#' Shut down a warmed pool
#'
#' @param handles From [.ctBackendWarmWorkers()].
#' @keywords internal
.ctBackendWarmStop <- function(handles) {
  if (!is.null(.ct_warm_state$connections)) {
    options(future.connections.onMisuse = .ct_warm_state$connections)
  } else {
    options(future.connections.onMisuse = NULL)
  }
  .ct_warm_state$connections <- NULL
  if (!.ctBackendCanWarm()) return(invisible(NULL))
  # Each worker holds a Julia process, so leaving the pool up leaves those up
  # too. `sequential` releases the sessions.
  tryCatch(future::plan(future::sequential), error = function(e) NULL)
  .ct_warm_state$threads <- NULL
  invisible(NULL)
}

# Independent jobs in the warmed worker pool: cross-validation folds
# (`ctLOO(parallelFolds = TRUE)`), profile walks (`ctFitProfile(processes =
# TRUE)`). Separate fits parallelise better as processes than as threads inside
# one fit -- the subject loop's threads share one objective and its
# synchronisation, a process shares nothing -- so a job here runs whole in one
# worker with a one-thread Julia, warmed for the model exactly as a sampling
# chain's worker is (`.ctBackendWarmWorkers()`), and the pool persists for the
# next call the same way.
#
# `fn` names a ctsem namespace function of one argument, the job, which must
# carry everything the function needs: the worker has the installed package and
# the job, nothing else. A job that can fail should catch its own failure and
# say so in its value. Returns the values in job order, or `NULL` when the pool
# could not be used or a worker did not answer, for the caller to run the jobs
# in this session instead -- the contract the sampler's processes path keeps.
#' @keywords internal
.ctBackendWorkerMap <- function(object, jobs, fn, workers, values = NULL,
  verbose = FALSE) {
  workers <- min(suppressWarnings(as.integer(workers)[1L]), length(jobs))
  if (!isTRUE(workers >= 2L) || !.ctBackendCanWarm()) return(NULL)
  handles <- .ctBackendWarmWorkers(object, workers, values = values)
  if (.ctBackendWarmWait(handles, verbose = verbose) < 1L) return(NULL)
  futures <- lapply(jobs, function(job) tryCatch(
    future::future(utils::getFromNamespace(fn, "ctsem")(job), seed = TRUE),
    error = function(e) NULL))
  if (any(vapply(futures, is.null, logical(1L)))) return(NULL)
  out <- lapply(futures, function(f) tryCatch(future::value(f),
    error = function(e) structure(list(), class = "ctWorkerMapFailure",
      message = conditionMessage(e))))
  if (any(vapply(out, inherits, logical(1L), "ctWorkerMapFailure"))) return(NULL)
  out
}

#' Release the warmed sampling worker pool
#'
#' \code{\link{ctFitUncertainty}} with \code{uncertainty = 'sample'} and
#' \code{ctFit(backend = 'julia', optimize = FALSE)} warm a pool of background
#' R processes ahead of a multi-chain sample, each
#' one compiled for the model's shape before its chain starts, so the compile
#' cost overlaps the optimisation that runs first rather than being paid
#' serially once sampling begins. That pool is deliberately left running
#' afterwards: \pkg{future}'s \code{multisession} plan persists across calls,
#' so a later multi-chain sample reuses the same warmed workers instead of
#' recompiling. Each worker holds its own Julia session, so it is also memory a
#' caller may want back sooner than the end of the R session, when it is
#' released automatically.
#'
#' Calling this in the middle of a sample that is still using the pool stops
#' the workers that sample is running in.
#'
#' @return \code{NULL}, invisibly.
#' @seealso \code{\link{ctFitUncertainty}}, \code{\link{ctJuliaSetup}}
#' @export
ctJuliaWorkersStop <- function() {
  invisible(.ctBackendWarmStop(NULL))
}
