# Doubly adaptive importance sampling: the normal `'checked'` holds, moved to
# the posterior's moments before anything is sampled where its check named a
# direction. The method, its stopping rule and what each rule cost to learn
# are in the engine (`ctsem_dais`, inst/julia/ContinuousTimeSEM/src/dais.jl),
# which runs it, and nothing crosses the bridge until the answer does. Where
# `'checked'` uses a coreset (`.ctCoresetLpg()`), DAIS's draws go through one
# of its own, on its units once one holds; why its own and not the importance
# sampler's is in the engine's comment. Otherwise they are evaluated at once
# on the pool or on copies of the objective.
#
# It ran here, in R, until 2026-10-24, beside the rest of `'checked'`: one
# bridge call an iteration for its draws, its arithmetic in R, and each
# iteration's draws fed to the importance sampler's coreset, whose refits were
# then R's QR least squares. On the 266-parameter scale fit 87 per cent of
# `'checked'`'s time went to those refits and `'checked'` had not finished in
# nine and a half hours. The coreset and its refits are in the engine now too
# (coreset_sampler.jl).
#
# So it needs the engine: a julia fit's density carries it (`'engine'`,
# `.ctBackendLpgFunc()`), and a density without one -- a stan fit's, or a
# closed form a test builds in R -- has no DAIS stage, and `'checked'` goes on
# from the check's normal as it did before.

# The DAIS stage on `lpg`'s engine objective from `N(centre, cov)`, the
# directions in `hold` held at the normal's width; `NULL` where the density
# has no engine. Where `lpg` is a coreset density (`.ctCoresetLpg()`) the
# draws go through a coreset of DAIS's own. The record has the fields
# `'checked'` reads: `used`, the normal (`mean`, `cov`), iterations, averaged,
# settled, the last damping and effective sample, the unusable draws, the
# draws an iteration, gradients, `cost` in full evaluations (a gradient as
# two), what the coreset did, and the trace.
#' @keywords internal
.ctDaisEngine <- function(lpg, centre, cov, control = list(), verbose = 0, hold = NULL) {
  engine <- attr(lpg, 'engine')
  if (!is.function(engine)) return(NULL)
  e <- engine()
  overwrite <- .ctProgressOverwrite(verbose)
  args <- list(e$objective, .ctJuliaNumericVector(as.numeric(centre)),
    .ctJuliaPut(as.matrix(cov)),
    draws = as.integer(.ctJuliaOr(control$daisDraws, 0L))[1L],
    maxiter = max(1L, as.integer(.ctJuliaOr(control$daisMaxIter, 150L))[1L]),
    # From R's stream, so that set.seed() makes 'checked' reproducible.
    seed = sample.int(.Machine$integer.max, 1L),
    rtol = .ctFlatDirectionRtol(), verbose = isTRUE(verbose > 0),
    # A coreset density's 'coreset' is its record: there is a coreset to use.
    coreset = is.function(attr(lpg, 'coreset')),
    progress_overwrite = overwrite,
    progress_sink = if (isTRUE(verbose > 0)) .ctProgressSink(overwrite) else NULL)
  if (length(hold) && NCOL(hold)) args$hold <- .ctJuliaPut(as.matrix(hold))
  r <- .ctBackendJuliaValue(do.call(e$module$ctsem_dais, args))
  trace <- as.matrix(r$trace)
  if (length(trace)) colnames(trace) <- as.character(unlist(r$trace_names))
  out <- list(used = isTRUE(r$used), diverged = isTRUE(r$diverged),
    iterations = as.integer(r$iterations), unusable = as.integer(r$unusable),
    draws = as.integer(r$draws), gradients = as.integer(r$gradients), values = 0L,
    cost = as.numeric(r$cost), trace = trace,
    coreset = list(used = isTRUE(r$coreset), size = as.integer(r$coreset_size),
      chosen = as.integer(r$coreset_chosen), failed = as.integer(r$coreset_failed),
      record = if (!is.null(r$coreset_record)) .ctCoresetRecord(r$coreset_record)))
  if (isTRUE(r$diverged)) out$growth <- as.numeric(r$growth)
  if (isTRUE(r$used)) {
    out$mean <- as.numeric(r$mean)
    out$cov <- as.matrix(r$cov)
    out$averaged <- as.integer(r$averaged)
    out$settled <- isTRUE(r$settled)
    out$gamma <- as.numeric(r$gamma)
    out$ess <- as.numeric(r$ess)
  }
  out
}
