# Doubly adaptive importance sampling: the normal `'checked'` holds, moved to
# the posterior's moments before anything is sampled where its check named a
# direction. The method, its stopping rule and what each rule cost to learn
# are in the engine (`ctsem_dais`, inst/julia/ContinuousTimeSEM/src/dais.jl),
# which runs it: draws evaluated at once on the pool or on copies of the
# objective, and nothing crossing the bridge until the answer does.
#
# It ran here, in R, until 2026-10-24, beside the rest of `'checked'`: one
# bridge call an iteration for its draws, its arithmetic in R, and each
# iteration's draws fed to the importance sampler's coreset, whose refits are
# R's QR least squares. On the 266-parameter scale fit 87 per cent of
# `'checked'`'s time went to those refits and `'checked'` had not finished in
# nine and a half hours. The coreset now learns from the importance sampler's
# own batches, as it did before DAIS existed.
#
# So it needs the engine: a julia fit's density carries it (`'engine'`,
# `.ctBackendLpgFunc()`), and a density without one -- a stan fit's, or a
# closed form a test builds in R -- has no DAIS stage, and `'checked'` goes on
# from the check's normal as it did before.

# The DAIS stage on `lpg`'s engine objective from `N(centre, cov)`, the
# directions in `hold` held at the normal's width; `NULL` where the density
# has no engine. The record has the fields `'checked'` reads: `used`, the
# normal (`mean`, `cov`), iterations, averaged, settled, the last damping and
# effective sample, the unusable draws, the draws an iteration, gradients,
# `cost` in full evaluations (a gradient as two), and the trace.
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
    progress_overwrite = overwrite,
    progress_sink = if (isTRUE(verbose > 0)) .ctProgressSink(overwrite) else NULL)
  if (length(hold) && NCOL(hold)) args$hold <- .ctJuliaPut(as.matrix(hold))
  r <- .ctBackendJuliaValue(do.call(e$module$ctsem_dais, args))
  trace <- as.matrix(r$trace)
  if (length(trace)) colnames(trace) <- as.character(unlist(r$trace_names))
  out <- list(used = isTRUE(r$used), diverged = isTRUE(r$diverged),
    iterations = as.integer(r$iterations), unusable = as.integer(r$unusable),
    draws = as.integer(r$draws), gradients = as.integer(r$gradients), values = 0L,
    cost = 2 * as.numeric(r$gradients), trace = trace)
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
