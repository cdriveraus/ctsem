# Random restarts for a julia fit that did not converge.
#
# Only for a fit the certification says is not a maximum, or could not certify:
# a fit that converged, or reached a maximum some coordinate of which the data
# do not determine, is left alone. And not when the user chose the start or
# capped the iterations, because then a non-converged fit is the fit that was
# asked for -- an earlier after-the-fact retry moved exactly such a fit, from
# raw 24 to raw 12.8, and reported it without comment.
#
# The case this is for: a likelihood with more than one basin, where the start
# decides which one the optimiser reaches. A rank-deficient two-random-effect
# laplace fixture (test-julia-multivariate-mixed.R) stops on a ridge 2 nats
# below a maximum in the other basin; nothing local reaches it, and a handful
# of dispersed starts might. Off unless asked for (optimcontrol$restarts = n),
# sd 1 by default, and spread over processes when there are cores to spare.

.ctBackendRestartsWanted <- function(result, certification, optimcontrol,
  inits, intoverstates) {
  # Opt-in: each restart is a whole optimisation, which on a slow laplace fit
  # is minutes, and on the fixture below five of them found nothing better.
  n <- suppressWarnings(as.integer(.ctJuliaOr(optimcontrol$restarts, 0L))[1L])
  if (!isTRUE(n >= 1L)) return(0L)
  if (!isTRUE(intoverstates)) return(0L)
  if (!is.null(inits) || !is.null(optimcontrol$maxiter)) return(0L)
  if (!is.null(certification)) {
    if (isTRUE(certification$certified) ||
        identical(certification$status, "unidentified")) return(0L)
  } else if (isTRUE(result$converged)) return(0L)
  n
}

# One restart, in whichever process runs it: the engine optimiser from `from`
# with the fit's own controls, reporting nothing.
.ctBackendRestartOne <- function(spec, from, optimcontrol, gradient) {
  tryCatch({
    if (is.null(.ct_julia_cache$module)) ctsem::ctJuliaSetup(threads = 1)
    r <- suppressWarnings(suppressMessages(.ctJuliaOptimise(spec, from,
      optimcontrol = optimcontrol, gradient = gradient, cores = 1L,
      verbose = 0L)))
    r
  }, error = function(e) structure(list(), error = conditionMessage(e)))
}

# Run the restarts; returns the best result (or NULL), the table, and whether
# the user stopped them. The table always has one row per start asked for,
# NA where a start did not run.
.ctBackendRestarts <- function(spec, start, npar, n, optimcontrol, gradient,
  cores = 1L, current = -Inf, report = FALSE) {
  sd <- as.numeric(.ctJuliaOr(optimcontrol$restartsd, 1))[1L]
  starts <- lapply(seq_len(n), function(i) as.numeric(start) + stats::rnorm(npar, 0, sd))
  # Nothing a worker cannot carry, and nothing that would report from inside it.
  control <- optimcontrol
  control$callback <- NULL
  control$progress <- FALSE
  control$restarts <- 0L
  table <- data.frame(start = seq_len(n), logposterior = NA_real_,
    converged = NA, iterations = NA_integer_, error = NA_character_,
    stringsAsFactors = FALSE)
  results <- vector("list", n)
  record <- function(i, r) {
    if (is.null(r) || !length(r) || is.null(r$maximum_loglik)) {
      # Kept, not dropped: a restart that fails silently is how a pool that
      # could not run any of them went unnoticed.
      table$error[i] <<- if (!is.null(attr(r, "error"))) attr(r, "error") else
        "no result"
      return()
    }
    results[[i]] <<- r
    table$logposterior[i] <<- as.numeric(r$maximum_loglik)[1L]
    table$converged[i] <<- isTRUE(r$converged)
    table$iterations[i] <<- as.integer(.ctJuliaOr(r$iterations, NA_integer_))
    if (report) message(sprintf("  restart %d of %d: logpost %.4f%s", i, n,
      table$logposterior[i], if (isTRUE(r$converged)) "" else ", not converged"))
  }
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  workers <- min(as.integer(n), as.integer(cores))
  cancelled <- FALSE
  pool <- NULL
  # Where the restarts are running when an interrupt lands: in workers (this
  # session's Julia is idle) or here (it may be mid-call).
  here <- TRUE
  tryCatch({
    if (workers > 1L && .ctBackendCanWarm()) {
      pool <- .ctBackendWarmWorkers(handle, workers = workers, values = start)
    }
    if (!is.null(pool) && .ctBackendWarmWait(pool) >= 1L) {
      here <- FALSE
      jobs <- lapply(seq_len(n), function(i) tryCatch(future::future(
        utils::getFromNamespace(".ctBackendRestartOne", "ctsem")(spec,
          starts[[i]], control, gradient), seed = TRUE),
        error = function(e) NULL))
      for (i in seq_len(n)) {
        if (!is.null(jobs[[i]])) record(i, tryCatch(future::value(jobs[[i]]),
          error = function(e) structure(list(), error = conditionMessage(e))))
      }
      # None ran in the workers -- a worker that loaded a different ctsem
      # build than this session, for one -- so run them here rather than
      # report five failures as a result.
      if (all(is.na(table$logposterior))) {
        message("Restarts could not run in worker processes (",
          table$error[1L], "); running them in this session.")
        table$error <- NA_character_
        here <- TRUE
        for (i in seq_len(n)) record(i, .ctBackendRestartOne(spec, starts[[i]],
          control, gradient))
      }
    } else {
      for (i in seq_len(n)) record(i, .ctBackendRestartOne(spec, starts[[i]],
        control, gradient))
    }
  }, interrupt = function(cnd) {
    cancelled <<- TRUE
    if (!is.null(pool)) try(.ctBackendWarmStop(NULL), silent = TRUE)
    if (!here) {
      # This session's Julia was idle throughout, so only the workers go.
      message("Restarts stopped; keeping the fit as found.")
    } else {
      # Interrupted inside a Julia call in this session, which leaves the
      # bridge waiting for a reply that belongs to the call it abandoned.
      try(.ctJuliaClearSession(), silent = TRUE)
      message("Restarts stopped; keeping the fit as found. The Julia session ",
        "was restarted, so the rest of this fit recompiles for its model shape.")
    }
  })
  best <- which.max(replace(table$logposterior, is.na(table$logposterior), -Inf))
  keep <- length(best) && is.finite(table$logposterior[best]) &&
    table$logposterior[best] > current + 1e-3
  table$chosen <- FALSE
  if (keep) table$chosen[best] <- TRUE
  list(best = if (keep) results[[best]] else NULL, table = table,
    cancelled = cancelled)
}
