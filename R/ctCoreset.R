# A coreset for importance sampling on the Laplace route.
#
# A Laplace fit's log density is a sum of per-unit terms, each an inner mode
# solve, and an importance sampler pays for every unit at every draw. Over the
# region the sampler visits, the units' terms vary in about as many dimensions
# as the sampler moves -- to first order a unit's term is its score times the
# displacement -- however many units there are. So M of the N units with
# non-negative weights and a free intercept, fitted by least squares to the
# full sum, follow it, and the sampler can pay M/N of a full evaluation per
# draw. On the strong-shell model of review/UNCERTAINTY-coreset-test-2026-10-21.md,
# 30-40 units carried the sum of 100, 250 or 1,000 along 'checked''s directions,
# and the sampler reached the full sum's posterior in 2.7, 4.9 and 10 times
# fewer evaluations.
#
# The coreset is the engine's (`CTSEMISCoreset`, coreset_sampler.jl): its
# training rows -- every unit's term at the draws evaluated in full, and at a
# few of each full gradient batch's draws every unit's score -- stay where they
# are computed, it chooses, checks and corrects there, and its selection is
# the one the optimiser's stages and SMC use (coreset_select.jl). Every decision
# is by the expected cost of the sampler's work in full evaluations, never in
# seconds, so the same fit makes the same choices anywhere. What is here wraps
# it as the sampler's density, a batch a call.

# The sampler's density with a coreset behind it: `lpg` (`.ctBackendLpgFunc()`)
# with 'batch', 'gradbatch', 'fullgradbatch' (the full density's gradient, for
# a check of the answer, `.ctSteinCheck()`, that a coreset's error must not
# reach), 'cost' (full evaluations spent, a gradient as two), 'finalize' (the
# end check, which corrects the weights to the full density where the coreset
# strayed; see `imis_is()`) and 'coreset' (the record
# `$details$importance_sampling$coreset` reports). `lpg` itself comes back
# where a coreset cannot apply: off the Laplace route, `control$coreset` FALSE,
# or a density that is not a sum over the fit's units (its 'unitsum' mark).
# The coreset is built on the density's own engine object (its 'engine'), which
# says what a unit's term is: the fit's Laplace terms (`.ctBackendLpgFunc()`),
# the one density that carries the mark.
.ctCoresetLpg <- function(fit, lpg, control = list()) {
  if (isFALSE(control$coreset) || !.ctBackendIsLaplace(fit) ||
    !isTRUE(attr(lpg, 'unitsum')) || !is.function(attr(lpg, 'batch')) ||
    !is.function(attr(lpg, 'gradbatch')) || !is.function(attr(lpg, 'engine'))) {
    return(lpg)
  }
  engine <- local({
    built <- NULL
    function() {
      if (is.null(built)) {
        density <- attr(lpg, 'engine')()
        built <<- list(module = density$module,
          cs = density$module$ctsem_is_coreset(density$objective))
      }
      built
    }
  })
  # The engine takes draws as columns.
  call <- function(name, ...) {
    e <- engine()
    .ctBackendJuliaValue(e$module[[name]](e$cs, ...))
  }
  columns <- function(X) t(as.matrix(X))
  gradients <- function(name, X) {
    X <- as.matrix(X)
    res <- call(name, columns(X))
    list(value = as.numeric(res$value),
      gradient = t(matrix(as.numeric(res$gradient), ncol = nrow(X))))
  }
  batch <- function(draws) as.numeric(call('ctsem_is_coreset_batch!', columns(draws)))
  gradbatch <- function(draws) gradients('ctsem_is_coreset_gradbatch!', draws)
  gradfull <- function(draws) gradients('ctsem_is_coreset_fullgrad!', draws)
  finalize <- function(samples, logp, logw) {
    out <- call('ctsem_is_coreset_finalize!', columns(samples),
      .ctJuliaNumericVector(as.numeric(logp)), .ctJuliaNumericVector(as.numeric(logw)))
    if (is.null(out)) NULL else as.numeric(out)
  }
  record <- function() .ctCoresetRecord(.ctBackendJuliaValue(
    engine()$module$ctsem_is_coreset_record(engine()$cs)))
  fn <- function(parm) batch(matrix(parm, 1L))
  attr(fn, 'batch') <- batch
  # The full objective's engine, for a stage that runs there on the full
  # density (`.ctDaisEngine()`); the coreset serves only the batches here.
  attr(fn, 'engine') <- attr(lpg, 'engine')
  attr(fn, 'gradbatch') <- gradbatch
  attr(fn, 'fullgradbatch') <- gradfull
  attr(fn, 'cost') <- function() as.numeric(call('ctsem_is_coreset_cost'))
  attr(fn, 'finalize') <- finalize
  attr(fn, 'coreset') <- record
  fn
}

# What a coreset did (`ctsem_is_coreset_record`, as the engine returns it), as
# the importance sampler's record of it and DAIS's report it.
.ctCoresetRecord <- function(r) {
  used <- isTRUE(r$used)
  choices <- as.matrix(r$choices)
  list(
    used = used,
    units = if (r$units > 0) as.integer(r$units),
    size = if (used) as.integer(r$size),
    sizes = as.integer(unlist(r$sizes)),
    phase = r$phase,
    switches = as.integer(r$switches),
    rounds_on_coreset = as.integer(r$rounds_on_coreset),
    failed_checks = as.integer(r$failed_checks),
    expected_round_cost = if (used) as.numeric(r$expected_round_cost),
    free_checks = as.integer(r$free_checks),
    chosen = as.integer(r$chosen),
    attempts = as.integer(r$attempts),
    chosen_passed = as.integer(r$chosen_passed),
    seconds = unlist(r$seconds),
    draws = unlist(r$draws),
    full_for = unlist(r$full_for),
    choices = if (length(choices)) `colnames<-`(choices, c('rows', 'M', 'why', 'fail',
      'cost', 'full', 'bestM', 'relheld')),
    training_draws = as.integer(r$training_draws),
    score_draws = as.integer(r$score_draws),
    reused_values = as.integer(r$reused_values),
    check_rms = as.numeric(r$check_rms),
    check_max = as.numeric(r$check_max),
    final = if (isTRUE(r$has_final)) list(checked = as.integer(r$final_checked),
      worst = as.numeric(r$final_worst), corrected = as.integer(r$final_corrected)),
    cost = as.numeric(r$cost))
}
