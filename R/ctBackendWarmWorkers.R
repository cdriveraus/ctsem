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
#' @return A list of future handles, or `NULL` when warming is unavailable.
#' @keywords internal
.ctBackendWarmWorkers <- function(object, workers, values = NULL) {
  workers <- suppressWarnings(as.integer(workers))
  if (!.ctBackendCanWarm() || !isTRUE(workers >= 1L)) return(NULL)
  if (!inherits(object, c("ctJuliaModel", "ctJuliaFit"))) return(NULL)
  npar <- .ctBackendNpar(object)
  if (!isTRUE(npar >= 1L)) return(NULL)
  if (is.null(values)) values <- rep(0, npar)
  values <- as.numeric(values)
  values <- rep_len(values, npar)
  values[!is.finite(values)] <- 0

  started <- tryCatch({
    future::plan(future::multisession, workers = workers)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(started)) return(NULL)

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
      future::future(ctsem:::.ctBackendWarmSession(object, values),
        seed = TRUE),
      error = function(e) NULL)
  })
  if (!length(handles) || all(vapply(handles, is.null, logical(1)))) return(NULL)
  attr(handles, "started") <- Sys.time()
  handles
}

# Runs inside a worker. Wrapped, because a worker that cannot warm should leave
# the parent to sample with threads rather than take the fit down with it.
#' @keywords internal
.ctBackendWarmSession <- function(object, values) {
  tryCatch({
    # Only when this process has no Julia yet. Asking again warns that the
    # thread count cannot be changed on a running session, which is true and
    # not worth saying: a warmed worker being asked for more work is the
    # expected case, not a problem.
    if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = 1)
    # One evaluation does both jobs: `ctJuliaEvaluate` builds the objective,
    # which marshals the data, and evaluating it forces the specialisation. The
    # objective cache is per process and keyed on content, and this process has
    # an empty one, so nothing here is shared with the parent.
    invisible(ctsem::ctJuliaEvaluate(object, values, gradient = TRUE))
    TRUE
  }, error = function(e) structure(FALSE, message = conditionMessage(e)))
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
  ok <- sum(vapply(results, isTRUE, logical(1)))
  if (isTRUE(verbose)) {
    elapsed <- if (is.null(started)) NA_real_ else
      as.numeric(difftime(Sys.time(), started, units = "secs"))
    message(ok, " of ", length(handles), " worker(s) ready",
      if (is.finite(elapsed)) paste0(" (", round(elapsed), "s since they started, ",
        "most of it alongside the fit)") else "", ".")
  }
  as.integer(ok)
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
  invisible(NULL)
}

# Was a sixth copy of the parameter count. `.ctBackendNpar` is the one
# definition; this stays as a name because tests and scratch scripts call it.
#' @keywords internal
.ctBackendWarmNpar <- function(object) .ctBackendNpar(object)
