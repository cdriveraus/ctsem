# The Laplace corrections: the continuation (`.ctLaplaceContinue()`), which
# climbs the quadrature hybrid from the estimate and draws the fit's
# uncertainty on that density, and the step correction before it. A request
# to sample against the quadrature -- `ctFitUncertainty(fit, 'is' or
# 'checked', control = list(likelihood = 'quadrature'))`,
# `.ctBackendUncertaintyQuadrature()` -- is the continuation with the method
# asked for. The separate importance sampler against the adaptive
# Gauss-Hermite posterior at every unit (`.ctLaplaceCorrect()`, to
# 2026-10-23) paid `nodes^k` member evaluations a unit a draw where the hybrid
# pays them only for the units the screen flags, moved the estimate by one
# Newton step where the continuation climbs, and reported no error for the
# density it sampled.

# The correction every Laplace fit gets by default ---------------------------
#
# `optimcontrol$laplace_correct` (default TRUE) runs `ctsem_laplace_autocorrect`
# at the end of an optimised `intoverpop='laplace'` fit: after the optimiser,
# after the curvature certification and its Newton continuation have converged
# the Laplace objective, and after the uncertainty stage has built the fit's
# Hessian and draws. It computes no Hessian of its own -- the step's metric is
# the one the standard errors came from, at the point it was evaluated -- and
# pays in tiers, each skipped when the one before finds nothing:
#
#   screen  one quadrature evaluation: the per-unit gap summed in absolute
#           value. Where every unit's likelihood is quadratic in its effects
#           Laplace is exact, the gap is rounding, and the fit is left
#           untouched to the bit.
#   step    `2 * npar` quadrature evaluations for the gradient, then a line
#           search on the quadrature objective. Accepted only on an increase.
#   repeat  up to `maxsteps` (3), same metric, until the predicted gain is
#           small or a step moved no parameter by `step_tol` (0.1) of its
#           standard error. On a 40-subject random-DRIFT fixture one step
#           took 80% of the quadrature objective's gain and three took 99%;
#           see review/LAPLACE-default-correction-2026-09-24.md.
#
# What moves when it applies, and what deliberately does not:
#
#   fit$estimate$raw             the corrected point
#   fit$estimate$loglik          the quadrature log likelihood there (and
#                                `logposterior`, `subject_loglik` with it);
#                                `loglik_laplace` keeps the Laplace value at
#                                the Laplace optimum, `loglik_method` says which
#   fit$estimate$rawposterior    shifted by the step: recentred, not reshaped
#   fit$estimate$cov, $se        unchanged -- the Laplace curvature, at
#                                `fit$uncertainty$evaluated_at`
#   fit$laplace$correction       what was done, and where
#
# `'quadrature'`, the default, runs a different stage in the same place: see
# `.ctLaplaceContinue()` below.
#
# Provenance of the constants (Appendix B of
# review/OPTIM-consolidation-plan-2026-09-25.md):
#
#   nodes      5   the quadrature's own default since ctLaplaceCheck: enough to
#                  locate the maximum in the cases tested, 9 to settle the
#                  value (quadrature.jl's docstring measured the gap on one
#                  40-subject dataset).
#   tolerance  0.01 nats, the screen's pass mark, summed absolutely over units.
#                  A tolerance in objective units rather than a measured value:
#                  the linear fixture screens at 3e-14, so any model with a
#                  nonlinear effect clears it by orders of magnitude.
#   maxsteps   3   measured, on the 40-subject nonlinear fixture only: three
#                  steps take 99% of the quadrature objective's gain and stop
#                  0.1 se short of the refine optimum
#                  (review/LAPLACE-default-correction-2026-09-24.md 3a).
#   gain_tol   1e-3 nats predicted, and
#   step_tol   0.1 se: the stopping rules, chosen so the next gradient is not
#                  paid for; set on the same fixture.
#   step       1e-3 raw, the finite-difference step of the gap gradient. Not
#                  tuned; differences of the gap, which is smooth and small.
#   material   0.1 se: when print() says the estimate moved. Presentation.
.ctLaplaceCorrectDefaults <- list(nodes = 5L, tolerance = 0.01, maxsteps = 3L,
  gain_tol = 1e-3, step_tol = 0.1, step = 1e-3, material = 0.1)

# Which correction this fit gets -- FALSE, "step" or "quadrature" -- refusing
# by name a request that cannot apply. FALSE is accepted anywhere, since it
# describes what every other route does. TRUE, and the default, are
# "quadrature", the continuation: on the known-broken cases of
# review/OPTIM-consolidation-plan-2026-09-25.md section 10 it ended higher in
# exact log likelihood than the step correction on 9 of 11 and level on the
# other 2, and never below the Laplace optimum, where the step correction
# ended 1.4 nats below it on the gated-gaps A14 config.
.ctLaplaceCorrectResolve <- function(optimcontrol, intoverpop, optimize,
  intoverstates) {
  value <- optimcontrol$laplace_correct
  explicit <- !is.null(value)
  valid <- (is.logical(value) && length(value) == 1L && !is.na(value)) ||
    (is.character(value) && length(value) == 1L &&
      value %in% c("step", "quadrature"))
  if (explicit && !valid) {
    stop("optimcontrol$laplace_correct must be TRUE or FALSE, or 'step' or ",
      "'quadrature'.", call. = FALSE)
  }
  if (explicit && isFALSE(value)) return(FALSE)
  method <- if (is.character(value)) value else "quadrature"
  why <- if (!identical(as.character(intoverpop)[1L], "laplace")) {
    paste0("applies to intoverpop='laplace' only; with intoverpop='",
      intoverpop, "' there is no Laplace term to correct")
  } else if (!isTRUE(optimize)) {
    paste0("corrects an optimised estimate, and a sampled fit ",
      "(optimize=FALSE) has none to move")
  } else if (!isTRUE(intoverstates)) {
    paste0("corrects the marginal Laplace estimate, which an ",
      "intoverstates=FALSE fit does not make")
  } else if (isTRUE(optimcontrol$estonly)) {
    "steps against the fit's Hessian, which optimcontrol$estonly skips"
  } else NULL
  if (is.null(why)) return(method)
  if (explicit) {
    stop("optimcontrol$laplace_correct ", why, ". Drop it.", call. = FALSE)
  }
  FALSE
}

# The raw estimate the Laplace objective was maximised at: the fit's estimate,
# unless the default correction moved it. Anything that is about the Laplace
# objective itself -- cross-validating it, importance sampling its posterior --
# reads this rather than `fit$estimate$raw`.
.ctLaplaceOptimum <- function(fit) {
  corr <- fit$laplace$correction
  if (isTRUE(corr$applied) &&
      length(corr$laplace_estimate) == length(fit$estimate$raw)) {
    return(as.numeric(corr$laplace_estimate))
  }
  as.numeric(fit$estimate$raw)
}

# Whether the default correction moved this fit's estimate.
.ctLaplaceIsCorrected <- function(fit) isTRUE(fit$laplace$correction$applied)

# `fit` with the uncertainty `measure(fit)` gives on the Laplace objective,
# for a fit whose continuation did not draw it on the hybrid (`ctFit()`).
# Where the continuation moved the estimate without reaching its fixed point
# (`.ctLaplaceContinue()`, `draws = "recentred"`) the stored curvature is the
# Laplace optimum's, so the method runs there and its draws move with the
# estimate, as the continuation moves the normal's; otherwise it runs at the
# estimate.
.ctLaplaceUncertaintyAtOptimum <- function(fit, measure) {
  at <- .ctLaplaceOptimum(fit)
  est <- as.numeric(fit$estimate$raw)
  if (identical(at, est)) return(measure(fit))
  laplace <- fit
  laplace$estimate$raw <- at
  laplace <- measure(laplace)
  note <- fit$uncertainty$details$laplace_correction
  fit$estimate$cov <- laplace$estimate$cov
  fit$estimate$se <- laplace$estimate$se
  fit$estimate$rawposterior <- sweep(laplace$estimate$rawposterior, 2L, est - at, "+")
  fit$uncertainty <- laplace$uncertainty
  fit$uncertainty$details$laplace_correction <- note
  fit$identifiability <- laplace$identifiability
  fit
}

.ctLaplaceAutoCorrect <- function(fit, cores = 1L, verbose = 0L,
  control = .ctLaplaceCorrectDefaults) {
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  hessian <- fit$uncertainty$hessian
  # The curvature at the estimate: evaluated there, or within the fit's
  # `hessianreuse` standard errors, where the optimiser's finish keeps the
  # Hessian it took at the hand-over (see `.ctBackendStoredHessian()`).
  usable <- !is.null(.ctBackendStoredHessian(fit, est))
  # A 1x1 NaN rather than an empty matrix: an empty one deadlocks the bridge,
  # and the engine reads any size mismatch as "no Hessian".
  if (!usable) hessian <- matrix(NaN, 1L, 1L)
  failed <- function(phrase) {
    warning("Laplace correction skipped: ", phrase, ". The uncorrected fit is ",
      "returned; see fit$laplace$correction.", call. = FALSE)
    fit$laplace$correction <- list(status = "failed", applied = FALSE,
      message = phrase, nodes = as.integer(control$nodes))
    fit
  }
  # The standard errors scale the stopping rule; never an empty vector across
  # the bridge, so an absent one is a vector of NaN, which the engine ignores.
  se <- as.numeric(fit$estimate$se)
  if (length(se) != npar) se <- rep(NaN, npar)
  se[!is.finite(se)] <- NaN
  chunks <- suppressWarnings(as.integer(fit$optim$chunks)[1L])
  if (is.na(chunks) || chunks < 1L) chunks <- max(1L, as.integer(cores)[1L])
  if (verbose > 0) message("Laplace correction: quadrature screen (",
    control$nodes, " nodes)")
  res <- try(.ctBackendWithMaxChunks(chunks, .ctJuliaGet(
    .ctJuliaModule(fit$model_spec$project)$ctsem_laplace_autocorrect(
      .ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(est),
      .ctJuliaPut(as.matrix(hessian)),
      nodes = as.integer(control$nodes), step = as.numeric(control$step),
      tolerance = as.numeric(control$tolerance),
      maxsteps = as.integer(control$maxsteps),
      gain_tol = as.numeric(control$gain_tol),
      scale = .ctJuliaNumericVector(se),
      step_tol = as.numeric(control$step_tol)))), silent = TRUE)
  if (inherits(res, "try-error")) {
    return(failed(paste0("the quadrature could not be evaluated (",
      .ctLaplaceTryReason(res), ")")))
  }
  status <- as.character(res$status)
  if (identical(status, "quadrature_failed")) {
    return(failed("the quadrature was not finite at the estimate"))
  }
  if (identical(status, "nonfinite_step")) {
    return(failed("the correction step was not finite"))
  }
  steps <- as.integer(res$steps)
  newest <- as.numeric(res$estimate)
  delta <- newest - est
  delta_se <- if (length(se) == npar) ifelse(se > 0, delta / se, NA_real_) else
    rep(NA_real_, npar)
  names(delta) <- names(delta_se) <- names(fit$estimate$se)
  record <- list(method = "step", status = status, applied = FALSE,
    nodes = as.integer(res$nodes), tolerance = as.numeric(res$tolerance),
    # sum over units of |quadrature - Laplace| at the Laplace optimum, and the
    # signed total the fit's log likelihood was off by there.
    screen = as.numeric(res$screen), gap = as.numeric(res$gap_start),
    laplace_estimate = est,
    loglik_laplace = as.numeric(fit$estimate$loglik),
    steps = steps,
    predicted_gain = if (steps > 0L) as.numeric(res$predicted)[seq_len(steps)] else
      numeric(0),
    alpha = if (steps > 0L) as.numeric(res$alpha)[seq_len(steps)] else numeric(0),
    first_delta = as.numeric(res$first_delta),
    delta = delta, delta_se = delta_se,
    dropped_directions = as.integer(res$dropped_directions),
    # Where the step's metric was evaluated: the Laplace optimum, whose
    # Hessian the standard errors also come from.
    hessian_at = if (usable) "laplace_estimate" else "none",
    seconds = c(screen = as.numeric(res$screen_seconds),
      total = as.numeric(res$seconds)))
  moved <- identical(status, "corrected") && all(is.finite(newest))
  # The log likelihood is the quadrature one wherever the screen found a gap,
  # whether or not a step was accepted. On AnomAuth's spurious optimum no step
  # raised the quadrature objective, and the Laplace value left on the fit was
  # 17 nats above the true optimum's while its quadrature value was 10 below:
  # the reported number is what makes such a fit comparable, so it cannot wait
  # on the estimate moving.
  units <- as.numeric(res$quadrature_units)
  quadrature_ok <- status %in% c("corrected", "no_gain", "no_hessian") &&
    is.finite(as.numeric(res$quadrature)) && all(is.finite(units))
  if (quadrature_ok) {
    record$loglik_quadrature <- sum(units)
    record$logposterior_quadrature <- as.numeric(res$quadrature)
    record$logposterior_laplace <- as.numeric(res$laplace)
    # Quadrature minus Laplace at the point the fit now reports.
    record$gap_reported <- as.numeric(res$quadrature) - as.numeric(res$laplace)
    fit$estimate$loglik_laplace <- record$loglik_laplace
    fit$estimate$loglik <- record$loglik_quadrature
    fit$estimate$logposterior <- record$logposterior_quadrature
    subjects <- as.numeric(res$quadrature_subjects)
    if (length(subjects) == length(fit$estimate$subject_loglik) &&
        all(is.finite(subjects))) {
      fit$estimate$subject_loglik <- subjects
    }
    fit$estimate$loglik_method <- "quadrature"
  }
  if (!moved) {
    fit$laplace$correction <- record
    return(fit)
  }
  # Applied. The draws are moved by the step and not redrawn, so they keep the
  # Laplace curvature's shape -- recentred, not reshaped, which is what
  # `ctFitUncertainty(fit, 'is', control = list(likelihood = 'quadrature'))`
  # is for.
  record$applied <- TRUE
  record$material <- isTRUE(max(abs(delta_se), na.rm = TRUE) >= control$material)
  record$draws <- "recentred"
  fit$estimate$raw <- newest
  post <- fit$estimate$rawposterior
  if (!is.null(post) && ncol(post) == npar) {
    fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
  }
  if (!is.null(fit$uncertainty)) {
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("recentred on the corrected estimate; the covariance is ",
        "the Laplace curvature at fit$uncertainty$evaluated_at"),
      nodes = record$nodes)
  }
  fit$laplace$correction <- record
  fit
}

# The quadrature continuation (`optimcontrol$laplace_correct = 'quadrature'`,
# the default) --
#
# Where `'step'` takes up to three Newton steps on the quadrature objective by
# finite differences, `'quadrature'` climbs it: from the Laplace optimum, with the
# engine's own optimiser, on an objective whose nodes are held fixed so that
# its gradient is exact (the Fisher identity; see laplace_continuation.jl), in
# rounds between which the nodes are re-placed. bigIRT's `laplaceRefine` is the
# model (../bigIRT/optimisation-review-2026-09.md, R3.4 and R3.5); what differs
# here, and why, is in review/OPTIM-consolidation-plan-2026-09-25.md, section
# 10.
#
#   screen   the rule at the Laplace optimum against Laplace, unit by unit, as
#            `'step'` screens. A fit that passes is left untouched to the bit.
#   flags    only the units the screen flags carry nodes: largest gap first,
#            until what is left carries at most the screen's tolerance. Every
#            other unit keeps its Laplace term and its exact gradient.
#   rounds   each a trust region in the Laplace fit's own standard errors:
#            the round's L-BFGS runs in the fit's identified directions,
#            whitened by its curvature, so it starts from the Newton step the
#            fit's Hessian implies and is confined to `radius` standard errors
#            of the round's centre -- none until a round has been rejected --
#            and ends where its rules have drifted (a rule's effective number
#            of nodes halved since its placement, `_continuation_drifted`),
#            which is as far as a run on frozen rules can be trusted; then the
#            nodes are re-placed where it ended.
#   target   the FIXED POINT of that iteration: the estimate at which the
#            gradient of the objective with its nodes placed there is zero.
#            With the nodes held, that gradient is the quadrature's own
#            estimate of the exact score -- a node-weighted average of the
#            per-node scores -- so the fixed point solves the quadrature
#            version of the exact likelihood equation. It is not the maximum
#            of the quadrature VALUE with the nodes re-placed at every point
#            (what the finite-difference refinement found), whose gradient
#            also carries how the rule's error moves with its centre. The two
#            differ by that term, and on the 40-subject fixture the fixed point
#            is the nearer to the exact optimum: 0.07 of a standard error
#            against 0.14 (the exact marginal by a dense grid; see
#            test-julia-laplace-continue.R). bigIRT's refinement targets the
#            same fixed point. It is not always the better point, which is
#            why the rounds only approach it while the re-placed value holds
#            (`accept`): where the two part, the estimate is the last point
#            that did not lose value, status `stalled`, and the fit's
#            certification says it is not stationary.
#   accept   a round is kept when the fixed-point residual falls -- half the
#            squared whitened gradient at the re-placed point, the Newton gain
#            it predicts -- AND the quadrature value with the nodes re-placed
#            there falls by no more than `value_tol`. Either failing, the
#            round's nodes go back and the radius shrinks to a quarter of the
#            step. Both conditions are measured, one at a time:
#              - the value alone stalls short of the fixed point, because near
#                it the re-placed value can fall along the score (by 2.6e-4
#                against a promised +5e-5 on the 40-subject fixture);
#              - the residual alone reaches a worse point. On the gated-gaps
#                A14 config (random drift, T0MEANS and CINT: three effects a
#                subject, so the soft rule) it took the residual from 1.3 to
#                0.007 while the re-placed value fell by 0.9 nats, and ended
#                at an exact log posterior of -398.81, where both conditions
#                stop at -397.97 (Laplace -399.23, the step correction -400.60).
#                Before `_continuation_stiff_rule` gave the soft rule's stiff
#                complement a full product rule it was worse: the complement
#                sat at one node, whose fixed-node gradient leaves out its log
#                determinant, and the estimate ended 4.9 exact nats below the
#                Laplace optimum.
#   no worse the estimate is only reported if the re-placed value there is not
#            below the value at the Laplace optimum by more than `value_tol`;
#            otherwise the fit keeps the Laplace optimum (status `no_gain`),
#            with the quadrature log likelihood there reported, as the step
#            correction does when no step raises its objective.
#   skip     before any round: one gradient gives the gain the fixed-node
#            model predicts from the Laplace optimum (the residual there);
#            below `skip_gain` nothing moves (status `skipped`).
#   stop     when that residual is below the fit's certification tolerance,
#            or after `rounds` kept rounds or `attempts` rounds in all, or
#            (status `resolution`) when what a round can promise is under the
#            frozen draws' resolution, twice their combined standard error,
#            and the last round rose by less: a gain the keep rule cannot
#            read. Zero where no unit is on draws.
#            `stop_gain` > 0 stops instead on predicted gain in nats: the
#            residual below it once the last kept round realised less (a
#            round raising the re-placed objective by more being kept even if
#            the residual rose), a region shrunk until it promises less, or a
#            kept round gaining less. It is off (0) by default: on the bench's
#            paired grid (review/bench, quadcost 2026-09-27) it cut the
#            correction's time to a quarter and gained up to 2.5 exact nats on
#            gB8 and 0.7 on gC8, but lost 0.46 to 0.52 on gD3, 0.026 on gC2,
#            0.031 on gN3 and 0.06 to 0.16 on the AnomAuth default starts
#            against the rounds run to the tolerance.
#   guard    a continuation that moves the quadrature objective by more than
#            max(50, N/2) nats is reverted to the Laplace optimum with a
#            warning, keeping the rejected point, as bigIRT does.
#   Hessian  of the hybrid at the final point, its nodes placed there: exact
#            for the flagged units (the Louis identity over their fixed
#            nodes, each node's Hessian by forward-over-reverse sweeps),
#            central differences for the units left to Laplace and the prior
#            (see `ctsem_laplace_continuation_hessian`), and only where the
#            rounds reached the fixed point; the
#            covariance and draws then come from it, not from the Laplace
#            curvature at another point. Where they stopped short, the
#            estimate is not stationary on that objective and its Hessian need
#            not be concave, so none is taken: the fit's own uncertainty stays
#            and its draws are recentred, as the step correction's are.
#
# What moves when it applies: `fit$estimate$raw` is the continuation's
# estimate; `loglik`, `logposterior` and `subject_loglik` the quadrature values
# there (every unit by the rule, nodes placed there), with `loglik_laplace`
# kept; `cov`, `se` and the draws, when the rounds converged, from the
# continuation's Hessian, which is `fit$uncertainty$hessian` with
# `evaluated_at` the new estimate, and whose certification is the fit's
# (`hessian_at` in the record says which, as the step correction's does). `fit$laplace$correction` says what
# ran.
#
# Provenance of the constants. `nodes`, `tolerance` and `material` are the
# step correction's (see `.ctLaplaceCorrectDefaults`). `product_maxdim` (2),
# `soft_tau` (3.5) and `soft_maxdirs` (2) are the engine's, explained in
# laplace_continuation.jl. `rounds` (10) is twice bigIRT's cap of 5: its units
# are flat, and a nested unit here holds its leaves' conditional modes fixed
# with its nodes, so the iteration contracts more slowly there (the nested
# fixture's residual fell by a factor of 0.3 to 0.7 a round, where the
# one-level fixture's fell by 20 to 300). `attempts` (15) bounds the rounds a
# trust region may reject. `radius` (`Inf`: no bound until a round is
# rejected, from which the region is a quarter of that step, doubling to at
# most `radius_max` = 16 se when a kept round ends on the boundary having
# gained at least three quarters of what its model promised and cut to a
# quarter of the step when a round is rejected, the textbook schedule). A
# first round held to one se stopped short wherever the hybrid's optimum lay
# several se off: on the bench's gB8 ten rounds were still gaining when the
# budget ran out, 1.85 exact nats short, and gN1's rounds shrank onto
# re-placement differences of 0.01 nats 1.3 short, while a climb on the frozen
# rules reached both in one run (the anytime grids of 2026-10-23). The
# unbounded round is that run, kept only where the re-placed value holds,
# and ended where the rules have drifted: unbounded without that, gA14's
# first round ran 40 se out on rules that no longer described the integrand
# and was kept 0.3 exact nats below where bounded rounds ended.
# Growing the region on the
# residual instead -- only when a round cut it to a quarter -- held gated-gaps
# D3 at 0.25 se for six rounds whose gains matched their model's to 2%, and
# the rounds ran out 0.7 exact nats short of the best-known point. `maxiter`
# (100) caps one round; bigIRT's whole continuation took 30 to 60
# evaluations. Five is enough when `stop_gain` is on (gated-gaps A1, N1, C8:
# gradients from 54, 52 and 90 to 34, 52 and 79, end points within 0.04
# exact nats), not with the rounds run to the certification tolerance.
# `guard` and `guard_per_subject` are bigIRT's max(50, N/2).
# `skip_gain` (5e-3 nats) was set on the optimiser bench's default Laplace
# cells (dev/lapcontinue/calibrate-cost.R, dev1, 2026-09-27, one start each,
# the AnomAuth cells from their spurious maxima). Below it were cf_mixed,
# cf_ordinal, mvmix, ordinal, cf_binary and gD1, whose whole correction gained
# 6e-6 to 4.9e-3 exact nats -- at most 1.7 times the prediction -- and moved
# the estimate at most 0.075 se; the smallest prediction above it was jflat's
# 0.07 (gain 0.078). It is also the gain of a whitened Newton step of 0.1 se,
# the move `material` calls worth a line in print(). On the paired grid the
# twelve cells it skips kept their exact log likelihood to 0.005 nats at 0.05
# to 0.24 of the correction's time.
# `value_tol` (1e-3 nats) is the step correction's `gain_tol`: a change in the
# objective below it is not one to act on either way.
# `rtol` (1e-8) is the identifiability report's: a direction the Laplace
# curvature does not identify is held where the fit left it. `maxdim` (5) is
# the widest product rule the correction places, explained with the engine's
# default; a unit wider than that can be scored by frozen draws instead, per
# block of its tree when it nests, `draws` (256) the first batch, doubling to
# a standard error of `draws_se` (0.05 nats) or `draws_max` (16384) draws,
# the blocks beneath another each to its share of that, also explained there.
# `draws_scope` says which: "forced", the default, only the units the climb
# named at a jump (`fit$laplace$jump`), every other wide unit
# keeping its Laplace term as before; "all" every one. Not "all"
# by default: a flagged unit's gradient is a reverse sweep per draw, and on the
# round-9 fit (40 people of eight effects) scoring every unit so made the
# correction far dearer than the fit
# (review/LAPLACE-shoulders-frozen-draws-2026-10-21.md, section 9).
# `draws_df` (4) is the degrees of freedom of the proposal's components, a
# Student t whose tails cover what a Gaussian's at the mode missed (0.6 to 1.0
# nats over gA1's 40 units with the replicates' spread saying 0.1), `Inf` for
# Gaussian; only a unit's top block takes it.
# `check` (TRUE) runs the shift check first (`.ctLaplaceShiftCheck()`): every
# unit, whatever its width, on `check_draws` (128) frozen draws at each of
# `check_seeds` (4) seeds. Where the step correcting Laplace would take is
# resolved under `material` the rounds do not run; otherwise the units
# carrying it are forced into the rounds, on product rules up to `maxdim`
# effects and on draws beyond. 128 draws resolved the step to 0.05 se
# on the bench's narrow fixtures and 32 did not (0.17-0.45 se of noise where
# nothing moved); four seeds give the noise, which at ten effects a unit with
# a weakly identified covariance was larger than the step itself
# (scripts/lapcheck-2026-10/highdim_check.R).
.ctLaplaceContinueDefaults <- list(nodes = 5L, tolerance = 0.01,
  product_maxdim = 2L, soft_tau = 3.5, soft_maxdirs = 2L, rounds = 10L,
  attempts = 15L, radius = Inf, radius_max = 16, maxiter = 100L,
  material = 0.1, guard = 50, guard_per_subject = 0.5, rtol = 1e-8,
  value_tol = 1e-3, maxdim = 5L, stop_gain = 0, skip_gain = 5e-3,
  draws = 256L, draws_se = 0.05, draws_max = 16384L, draws_df = 4,
  draws_scope = "forced", check = TRUE, check_draws = 128L, check_seeds = 4L)

# How far the hybrid's objective is from the exact marginal, from the engine's
# info at a placement (`ctsem_laplace_continuation_info`): `nats`, the frozen
# rules' standard errors in quadrature plus the product rules' node errors
# plus what the screen leaves on Laplace; `se`, how far the estimate can be
# from the exact-marginal optimum in standard errors, by the certification's
# own arithmetic (a displacement costing `e` nats on the quadratic is at most
# `sqrt(2 e)` long); the three parts; each unit's own error (`NaN` where its
# term is Laplace's); and `unscored`, units too wide to be scored at all,
# which the number does not cover.
.ctLaplaceIntegrationError <- function(info, unscored = 0L) {
  nats <- as.numeric(info$error)[1L]
  # Rounding is not an error: on a linear Gaussian fit the screen's gaps are
  # 1e-14 nats of arithmetic, and reported as such they read as a finding.
  if (is.finite(nats) && nats < 1e-8) nats <- 0
  list(nats = nats,
    se = if (is.finite(nats) && nats >= 0) sqrt(2 * nats) else NA_real_,
    draws = as.numeric(info$error_draws)[1L],
    nodes = as.numeric(info$error_nodes)[1L],
    residual = as.numeric(info$error_residual)[1L],
    units = as.numeric(info$centre_error),
    unscored = as.integer(.ctJuliaOr(unscored, 0L))[1L])
}

# The engine arguments of a continuation built under `control`: the product
# rules' constants, the frozen importance rule's for units too wide for one
# (`draws`, `draws_se`, `draws_max`; laplace_continuation.jl's module
# docstring), the units the caller names to flag whatever their gap (`force`),
# points it knows for them (`points`, a list of `units` and a `values` matrix
# with a column each), and whether every unit is scored (`screen`). Optional
# entries are left out rather than sent empty, which would hang the bridge.
.ctLaplaceContinuationArgs <- function(control, force = NULL, points = NULL,
  screen = TRUE) {
  args <- list(nodes = as.integer(control$nodes),
    tolerance = as.numeric(control$tolerance),
    product_maxdim = as.integer(control$product_maxdim),
    soft_tau = as.numeric(control$soft_tau),
    soft_maxdirs = as.integer(control$soft_maxdirs),
    maxdim = as.integer(control$maxdim),
    draws = as.integer(.ctJuliaOr(control$draws, 0L)),
    draws_se = as.numeric(.ctJuliaOr(control$draws_se, 0.05)),
    draws_max = as.integer(.ctJuliaOr(control$draws_max, 16384L)),
    draws_df = as.numeric(.ctJuliaOr(control$draws_df, 4)),
    draws_scope = as.character(.ctJuliaOr(control[["draws_scope"]], "forced")),
    screen = isTRUE(screen))
  force <- suppressWarnings(as.integer(force))
  force <- force[!is.na(force) & force > 0L]
  if (length(force)) args$force <- .ctJuliaVector(force)
  if (!is.null(points)) {
    units <- as.integer(points$units)
    keep <- !is.na(units) & units > 0L
    if (any(keep)) {
      values <- as.matrix(points$values)[, keep, drop = FALSE]
      storage.mode(values) <- "double"
      args$point_units <- .ctJuliaVector(units[keep])
      args$point_values <- values
    }
  }
  args
}

# Warns, once, when more than a quarter of the units are integrated by draws.
# A unit on draws costs about a dozen Laplace units a gradient (2048 draws of
# one member against the inner solves of 40 units on the round-9 fit, local),
# so past that share each gradient is several times the Laplace route's, and
# sampling the posterior (`optimize = FALSE`) is the honest route. Returns
# whether it has warned.
.ctLaplaceDrawsShareWarn <- function(units, nunits, warned) {
  if (isTRUE(warned) || nunits < 4L || length(units) <= 0.25 * nunits) {
    return(isTRUE(warned))
  }
  warning(sprintf(paste0("%d of %d units' Laplace terms jumped and are ",
    "integrated by draws, which costs more a gradient than the Laplace route ",
    "saves; sampling the fit (optimize = FALSE) may be the better route."),
    length(units), nunits), call. = FALSE)
  TRUE
}


# The directions a round moves in: the Laplace curvature's identified ones,
# each scaled to one of its standard errors, so the round's L-BFGS starts from
# that curvature and its trust region is measured in standard errors.
#
# Left out, and held where the fit left them: the directions the likelihood was
# measured flat along (`flat`, the columns of
# `fit$uncertainty$details$flatdirections$vectors`), and those whose curvature
# is below `rtol` of the largest -- the identifiability report's two rules for
# a direction the data do not identify. A direction of negative curvature is
# kept, whitened by its magnitude: the quadrature objective may rise along it,
# which is how a continuation leaves a saddle of the Laplace objective.
.ctLaplaceContinueBasis <- function(hessian, flat = NULL,
  rtol = .ctLaplaceContinueDefaults$rtol) {
  information <- -(hessian + t(hessian)) / 2
  n <- nrow(information)
  projector <- diag(n)
  measured <- 0L
  if (is.matrix(flat) && nrow(flat) == n && ncol(flat) > 0L &&
      all(is.finite(flat))) {
    decomposition <- qr(flat)
    q <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
    projector <- projector - q %*% t(q)
    measured <- decomposition$rank
  }
  kept <- projector %*% information %*% projector
  e <- eigen((kept + t(kept)) / 2, symmetric = TRUE)
  scale <- max(abs(e$values))
  if (!is.finite(scale) || scale <= 0) return(NULL)
  keep <- abs(e$values) > rtol * scale
  basis <- e$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(abs(e$values[keep])), nrow = sum(keep), ncol = sum(keep))
  list(basis = basis, kept = sum(keep), frozen = n - sum(keep),
    measured_flat = measured, negative = sum(e$values[keep] < 0))
}

# The Laplace term checked on every unit, whatever its width
# (review/PROPOSAL-laplace-check-every-unit-2026-10-27.md): each unit on
# `draws` frozen importance draws at the estimate, at each of `seeds`
# independently (`ctsem_laplace_check`), and the Newton step the corrected
# gradient takes in the rounds' basis -- how far correcting Laplace would move
# each estimate, in its standard errors. The seeds are replicates, so the step
# is their mean and its noise their spread over sqrt(seeds): at ten effects a
# unit with a weakly identified covariance one seed's step was mostly noise.
# `status` is "passed" where every parameter's step is resolved under `bar`
# (its size plus twice its noise), "flagged" where one is resolved over it,
# and "unresolved" where the noise decides neither. `carrying` are the units
# whose own steps account for the step, largest first, until what is left is
# under `bar`. `gap` is the corrected minus the Laplace log likelihood, with
# its standard error. NULL where the fit has no curvature to step on or fewer
# than two seeds scored.
.ctLaplaceShiftCheck <- function(fit, draws = 128L, seeds = 1:4, bar = 0.1,
  control = .ctLaplaceContinueDefaults) {
  est <- as.numeric(fit$estimate$raw)
  hessian <- .ctBackendStoredHessian(fit, est)
  if (is.null(hessian) || !all(is.finite(hessian))) return(NULL)
  basis <- .ctLaplaceContinueBasis(hessian,
    fit$uncertainty$details$flatdirections$vectors, rtol = control$rtol)
  if (is.null(basis)) return(NULL)
  # In the curvature's widths along the directions it identifies, the basis
  # the step is taken in. A coordinate the fit reports with a near-zero
  # standard error because it lies where there is no curvature has no width
  # here either, and no step: divided by the reported one, a 20-effect fit
  # with most directions held showed steps of thousands of standard errors.
  B <- basis$basis
  se <- sqrt(rowSums(B^2))
  held <- !is.finite(se) | se <= 1e-8 * max(se[is.finite(se)], 0)
  se[held] <- 1
  module <- .ctJuliaModule(fit$model_spec$project)
  started <- proc.time()[["elapsed"]]
  r <- .ctJuliaGet(module$ctsem_laplace_check(.ctJuliaObjective(fit, hybrid = TRUE),
    .ctJuliaNumericVector(est), draws = as.integer(draws),
    seeds = .ctJuliaVector(as.integer(seeds)), draws_df = as.numeric(control$draws_df)))
  G <- matrix(as.numeric(r$gradient), nrow = length(est))
  scored <- apply(G, 2L, function(g) all(is.finite(g)))
  if (sum(scored) < 2L) return(NULL)
  project <- B %*% t(B)
  project[held, ] <- 0
  Z <- (project %*% G[, scored, drop = FALSE]) / se
  shift <- rowMeans(Z)
  noise <- apply(Z, 1L, stats::sd) / sqrt(ncol(Z))
  status <- if (all(abs(shift) + 2 * noise < bar)) "passed" else
    if (any(abs(shift) - 2 * noise >= bar)) "flagged" else "unresolved"
  # Each unit's own step, units in rows; the units carrying the step taken in
  # order of their largest component until what is left is under `bar`.
  U <- matrix(as.numeric(r$unit_gradient), ncol = length(est))
  steps <- t((project %*% t(U)) / se)
  order <- order(-apply(abs(steps), 1L, max))
  left <- shift
  carrying <- integer()
  for (u in order) {
    if (max(abs(left)) < bar) break
    carrying <- c(carrying, u)
    left <- left - steps[u, ]
  }
  laplace_units <- as.numeric(r$laplace_unit)
  checked <- matrix(as.numeric(r$check_unit), nrow = length(laplace_units))
  gaps <- colSums(checked[, scored, drop = FALSE]) - sum(laplace_units)
  list(status = status, shift = shift, noise = noise, bar = bar,
    max_shift = max(abs(shift)), worst = which.max(abs(shift)),
    carrying = carrying, gap = mean(gaps), gap_se = stats::sd(gaps) / sqrt(length(gaps)),
    units = length(laplace_units), draws = as.integer(draws), seeds = sum(scored),
    width = as.integer(r$unit_width),
    seconds = proc.time()[["elapsed"]] - started)
}

# What a fit keeps of the shift check (`fit$laplace$correction$check`): the
# verdict, the largest step and its noise by parameter name, the Laplace gap,
# and how many units carried the step. NULL where no check ran.
.ctLaplaceShiftRecord <- function(check, fit) {
  if (is.null(check)) return(NULL)
  names <- .ctBackendRawParameterNames(fit, length(check$shift))
  list(status = check$status, max_shift = check$max_shift,
    worst = names[check$worst], noise = check$noise[check$worst],
    shift = stats::setNames(check$shift, names),
    shift_noise = stats::setNames(check$noise, names), bar = check$bar,
    gap = check$gap, gap_se = check$gap_se, carrying = length(check$carrying),
    draws = check$draws, seeds = check$seeds, seconds = check$seconds)
}

# The rounds, shared by the fit (`.ctLaplaceContinue()`) and by
# `ctLaplaceCheck(refine = TRUE)`. `cont` is the engine's continuation object,
# placed at `est`; on return it is placed at the point returned. See the notes
# above `.ctLaplaceContinueDefaults` for what a round is, what it aims at and
# when one is kept. `residual` is half the squared whitened gradient with the
# nodes placed at the point: the Newton gain it predicts, zero at the answer.
.ctLaplaceContinueRun <- function(module, cont, est, basis, tol,
  control = .ctLaplaceContinueDefaults, verbose = 0L, skip_gain = 0,
  sink = NULL) {
  get <- .ctJuliaGet
  # `.ctLaplaceContinue()` passes its own sink, so the rounds continue the
  # same in-place line the screen started. `ctLaplaceCheck(refine = TRUE)`
  # passes none, so one is built here -- both report at the same default
  # verbosity rather than only when `verbose > 0`.
  if (is.null(sink)) sink <- .ctBackendProgressSink(verbose)
  x <- as.numeric(est)
  # Every stopping decision is a predicted gain in nats against one bar:
  # `stop_gain`, or the fit's certification tolerance where `stop_gain` is 0,
  # which is how the rounds ran before 2026-09-27 (they then stopped only on
  # that tolerance, 1e-6 by default, or by running out of rounds or attempts).
  stop_gain <- as.numeric(.ctJuliaOr(control$stop_gain, 0))
  bar <- if (stop_gain > 0) stop_gain else as.numeric(tol)
  start <- get(module$ctsem_laplace_continuation_info(cont))
  value <- as.numeric(start$quadrature)
  radius <- as.numeric(control$radius)
  # Where flagged units are integrated by frozen draws the value carries their
  # error, and a change within twice its combined standard error is no
  # change the rounds can read (the draws are the same at every placement,
  # so this overstates the error of a difference): the keep rule's allowance,
  # and the least gain worth a round.
  resolution_at <- function(info) 2 * sqrt(sum(as.numeric(info$centre_se)^2,
    na.rm = TRUE))
  resolution <- resolution_at(start)
  kept <- 0L
  attempts <- 0L
  status <- "rounds"
  rows <- list()
  B <- .ctJuliaPut(as.matrix(basis))
  # The engine reads a radius of 0 as no bound (`_continuation_inside`).
  optimise <- function(from, stationary = FALSE, tol = bar) get(
    module$ctsem_laplace_continuation_optimize(cont, .ctJuliaNumericVector(from),
      B, if (is.finite(radius)) radius else 0,
      maxiter = as.integer(control$maxiter), tol = tol,
      stationary_only = stationary))
  row <- function(round, residual_after, keep, flagged, gain = NA_real_) data.frame(
    round = attempts, radius = radius,
    residual = as.numeric(round$start_gain), residual_after = residual_after,
    fixed_gain = as.numeric(round$value) - as.numeric(round$start_value),
    gain = gain,
    moved = as.numeric(round$moved), iterations = as.integer(round$iterations),
    f_calls = as.integer(round$f_calls), g_calls = as.integer(round$g_calls),
    drifted = isTRUE(round$drifted), kept = keep, flagged = as.integer(flagged))
  # The gain the fixed-node model predicts from the start: half the squared
  # whitened gradient, the Newton gain in nats. One gradient, and it decides
  # whether any round runs at all (`skip_gain`).
  residual <- as.numeric(optimise(x, stationary = TRUE)$start_gain)
  predicted <- residual
  # What the last kept round raised the re-placed objective by; none yet.
  realised <- 0
  if (!is.finite(residual)) {
    status <- "failed"
  } else if (residual < skip_gain) {
    status <- "skipped"
  } else repeat {
    # The residual is a Newton gain in the Laplace fit's metric, and where the
    # quadrature objective is flatter than that it understates what is left:
    # gated-gaps C2 stopped at a residual of 0.003 with 0.1 exact nats still to
    # gain. So a residual under the bar is converged only when the last round
    # also realised less than it; while rounds keep realising more, they go
    # on, each taking its step on the curvature its own L-BFGS measures.
    rising <- stop_gain > 0 && realised >= stop_gain
    if (residual < bar && !rising) { status <- "converged"; break }
    if (kept >= control$rounds) { status <- "rounds"; break }
    if (attempts >= control$attempts) { status <- "attempts"; break }
    # What the next round can promise inside its region: the Newton gain when
    # the whitened Newton step fits, else the gain of the whitened gradient's
    # step to the boundary on the same model. A region that rejections have
    # shrunk until it promises less than the bar is as done as a residual
    # below it; before this, those rounds ran on until the attempts ran out.
    reach <- sqrt(2 * residual)
    promise <- if (reach <= radius) residual else radius * reach - radius^2 / 2
    if (stop_gain > 0 && promise < stop_gain && !rising) {
      status <- "stalled"; break
    }
    if (promise < resolution && realised < resolution) {
      status <- "resolution"; break
    }
    attempts <- attempts + 1L
    # Below the bar, the round's own tolerance goes below the residual, so
    # that it takes its step rather than reporting itself stationary.
    round <- optimise(x, tol = if (residual < bar) residual / 100 else bar)
    residual <- as.numeric(round$start_gain)
    if (isTRUE(round$stationary)) {
      status <- "converged"
      rows[[length(rows) + 1L]] <- row(round, NA_real_, FALSE, start$nflagged)
      break
    }
    if (!is.finite(residual)) { status <- "failed"; break }
    xn <- as.numeric(round$minimizer)
    fixed <- as.numeric(round$value) - as.numeric(round$start_value)
    if (!is.finite(fixed) || fixed <= 0) {
      # The fixed-node model had a gradient to offer and no step along it: a
      # line search that found no decrease. There is nothing to re-place for.
      status <- "no_progress"
      rows[[length(rows) + 1L]] <- row(round, NA_real_, FALSE, start$nflagged)
      break
    }
    placed <- get(module[["ctsem_laplace_continuation_recentre!"]](cont,
      .ctJuliaNumericVector(xn)))
    after <- as.numeric(optimise(xn, stationary = TRUE)$start_gain)
    gain <- as.numeric(placed$quadrature) - value
    # Kept when the residual fell without the re-placed value falling, or --
    # whatever the residual did -- when that value rose by more than
    # `value_tol` (or `stop_gain` where that is on). On gated-gaps D3 two
    # rounds that raised it by 0.087 and 0.020 were rejected for a residual
    # that rose, the region shrank, and the rounds stopped 0.6 exact nats
    # short; with this, the paired grid of 2026-09-27 gained 0.9 to 2.5 exact
    # nats on gB8, 0.6 to 0.7 on gC8 and 2.3 on AnomAuth S1 from its spurious
    # start. A kept rise is a rise in the objective the fit reports, so this
    # cannot walk downhill, and each such round gains at least the bar, so it
    # cannot cycle.
    noise <- resolution_at(placed)
    rise <- max(stop_gain, as.numeric(control$value_tol)) + noise
    keep <- is.finite(after) && is.finite(gain) &&
      ((after < residual && gain >= -(as.numeric(control$value_tol) + noise)) ||
        gain > rise)
    rows[[length(rows) + 1L]] <- row(round, after, keep, placed$nflagged, gain)
    # At the same default verbosity as the rest of the correction: a round can
    # run for minutes (Charles's ordinal fixture, 2026-09-28), and this used to
    # print only at verbose > 0, so a default fit's screen showed nothing for
    # the whole of it. Kept as one detailed line rather than split into a
    # separate terse default and a separate detailed verbose>0 one -- it is
    # already a phrase, not a paragraph.
    if (!is.null(sink)) {
      sink(sprintf(paste0("Laplace continuation round %d: radius %.3g, ",
        "residual %.3g -> %.3g, value %+.4g, %s"), attempts, radius, residual,
        after, gain, if (keep) "kept" else "rejected"), "update")
    }
    realised <- if (keep) gain else 0
    if (keep) {
      kept <- kept + 1L
      x <- xn
      value <- as.numeric(placed$quadrature)
      resolution <- noise
      # Grown when the round ended on the boundary and the re-placed value
      # rose by three quarters or more of what the fixed-node model promised:
      # the model is right, and the region is what holds the rounds back.
      if (isTRUE(round$boundary) && gain >= 0.75 * fixed) {
        radius <- min(2 * radius, as.numeric(control$radius_max))
      }
      residual <- after
      # Kept, but the objective with its nodes re-placed rose by less than
      # the bar: the rounds are closing on the fixed point without gaining
      # anything a fit could report. On gated-gaps A1 eight such rounds in a
      # row each took 0.015 se and lost 2e-4 to 9e-4 nats, until the rounds
      # ran out.
      if (stop_gain > 0 && gain < stop_gain && residual >= bar) {
        status <- "stalled"; break
      }
    } else {
      get(module[["ctsem_laplace_continuation_revert!"]](cont))
      radius <- 0.25 * max(as.numeric(round$moved), 1e-12)
      if (radius < 1e-3) { status <- "stalled"; break }
    }
  }
  info <- get(module$ctsem_laplace_continuation_info(cont))
  list(x = x, status = status, rounds = kept, attempts = attempts,
    radius = radius, residual = residual, predicted = predicted, bar = bar,
    resolution = resolution,
    trace = if (length(rows)) do.call(rbind, rows) else NULL,
    start = start, info = info, value = value)
}

# Newton steps on the continuation's own Hessian, where it says more is left
# than the rounds' stopping rule did.
#
# The rounds stop on a residual whitened by the Laplace fit's curvature
# (`.ctLaplaceContinueRun()`), and the certification measures the gap on the
# continuation's own Hessian, formed once they stop. The two estimate the same
# quantity and can disagree by an order of magnitude: on the affect pilot's
# slow-component fit (528 units, 352 parameters) the rounds stopped at a
# residual of 6.4e-07 against a bar of 1e-06, the exact gap was 5.6e-06, and
# the fit reported `converged = FALSE` for a continuation that had met its own
# rule. The Hessian also predicts the step that closes the gap, so that step is
# taken, as the optimiser's finish takes one on its Hessian: at most
# `maxsteps`, each kept only where the objective with its nodes re-placed did
# not fall and the gap fell, and only while the estimate stays within
# `.ctBackendHessianReuse()` standard errors of where `hc` was evaluated, so
# the matrix still describes the curvature there. A step is a placement and a
# gradient, against the `2 npar` gradients the Hessian took.
#
# `cont` is placed at `x`, `final` is the objective and its gradient there and
# `value` its quadrature value. Returns where it ended, placed there, with
# `final`, `info` (NULL when no step was kept) and the gap before and after.
#' @keywords internal
.ctLaplaceContinueFinish <- function(module, cont, x, hc, final, value, bar,
  value_tol, maxsteps = 2L, reuse = .ctBackendHessianReuse()) {
  get <- .ctJuliaGet
  at <- x
  gap <- .ctBackendOptimGap(hc, as.numeric(final$gradient))
  before <- gap$gap
  info <- NULL
  steps <- 0L
  while (steps < maxsteps && isTRUE(gap$ok) && identical(gap$nnegative, 0L) &&
      isTRUE(gap$gap > bar)) {
    xn <- x + gap$step
    if (!isTRUE(.ctBackendHessianDistance(hc, at, xn) <= reuse)) {
      break
    }
    placed <- get(module[["ctsem_laplace_continuation_recentre!"]](cont,
      .ctJuliaNumericVector(xn)))
    after <- get(module$ctsem_laplace_continuation_evaluate(cont,
      .ctJuliaNumericVector(xn), gradient = TRUE))
    next_gap <- .ctBackendOptimGap(hc, as.numeric(after$gradient))
    rose <- as.numeric(placed$quadrature) - value
    if (!isTRUE(is.finite(rose) && rose >= -value_tol &&
        next_gap$gap < gap$gap)) {
      get(module[["ctsem_laplace_continuation_revert!"]](cont))
      break
    }
    steps <- steps + 1L
    x <- xn
    final <- after
    info <- placed
    value <- as.numeric(placed$quadrature)
    gap <- next_gap
  }
  list(x = x, final = final, info = info, value = value, steps = steps,
    gap_before = before, gap_after = gap$gap)
}

# A log-posterior function on the continuation's hybrid, in the `lpgFunc`
# contract `.ctBackendUncertainty()` takes, with the engine, batch and
# gradient-batch attributes `.ctBackendLpgFunc()`'s density carries, so that
# 'checked' runs DAIS and its Stein check here too and draws cross the bridge
# a batch at a time. Every point is evaluated with the flagged units' rules
# placed there (`ctsem_continuation_density`): read where they were placed at
# the estimate, the rules said less the further a draw moved, and gD1's draws
# put a population sd 2.9 posterior sds high.
.ctLaplaceContinueLpg <- function(module, cont) {
  engine <- local({
    built <- NULL
    function() {
      if (is.null(built)) built <<- list(module = module,
        objective = module$ctsem_continuation_density(cont))
      built
    }
  })
  batches <- .ctBackendEngineBatches(engine)
  function(fit, gradient = TRUE) {
    wantgrad <- isTRUE(gradient)
    fn <- function(parm) {
      at <- matrix(as.numeric(parm), nrow = 1L)
      if (!wantgrad) return(batches$batch(at))
      result <- batches$gradbatch(at)
      value <- result$value[1L]
      grad <- as.numeric(result$gradient)
      # `.ctBackendLpgFunc()`'s guard: a point the hybrid cannot use is a large
      # finite penalty with no gradient.
      if (!is.finite(value) || value <= -1e100 || any(!is.finite(grad))) {
        value <- -1e100
        grad <- rep(0, length(parm))
      }
      attributes(value) <- list(gradient = grad)
      value
    }
    # A sum over the units too -- the flagged units' rules, the rest's Laplace
    # terms -- which a coreset of them may stand in for (`.ctCoresetLpg()`),
    # on the density's own terms: every draw a placement of every flagged
    # unit's rule was most of 'checked''s time here (anomS1: 804 s against
    # 224 with the rules read where they were placed).
    attr(fn, "unitsum") <- TRUE
    attr(fn, "engine") <- engine
    attr(fn, "batch") <- batches$batch
    attr(fn, "gradbatch") <- batches$gradbatch
    fn
  }
}

# The log posterior the quadrature correction would report at `values`, from
# its own screen (`ctsem_laplace_continuation`): quadrature for every unit no
# wider than `maxdim` random effects, the Laplace term for the rest, plus the
# prior. `NA` when the quadrature cannot be evaluated. Used to compare two
# maxima of one Laplace objective on the objective the fit reports
# (`.ctJuliaOptimiseFit()`'s guard on `saem = TRUE`).
#' @keywords internal
.ctLaplaceCorrectedValue <- function(model_spec, values,
  control = .ctLaplaceContinueDefaults) {
  spec <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  cont <- try(do.call(module$ctsem_laplace_continuation,
    c(list(.ctJuliaObjective(spec, hybrid = TRUE), .ctJuliaNumericVector(values)),
      .ctLaplaceContinuationArgs(control))), silent = TRUE)
  if (inherits(cont, "try-error")) return(c(corrected = NA_real_,
    laplace = NA_real_, wide = NA_real_, units = NA_real_))
  info <- .ctJuliaGet(module$ctsem_laplace_continuation_info(cont))
  c(corrected = as.numeric(info$quadrature)[1L],
    laplace = as.numeric(info$laplace)[1L],
    wide = as.numeric(info$nwide)[1L], units = as.numeric(info$nunits)[1L])
}

# A try-error's message in a phrase: its first line, without the Julia
# stacktrace that follows an engine error.
# The first line that says what went wrong. An engine error arrives as
# "Evaluation in Julia failed." and then the engine's own message, which is
# the line a reader needs; the generic one is what a record used to keep.
.ctLaplaceTryReason <- function(x) {
  lines <- trimws(strsplit(paste(conditionMessage(attr(x, "condition")),
    collapse = "\n"), "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return("")
  own <- lines[!grepl("^Evaluation in Julia failed|^Original (Julia )?error message",
    lines)]
  reason <- if (length(own)) own[1L] else lines[1L]
  if (nchar(reason) > 300) reason <- paste0(substr(reason, 1, 300), "...")
  reason
}

.ctLaplaceContinue <- function(fit, cores = 1L, verbose = 0L,
  control = .ctLaplaceContinueDefaults, method = NULL, uncertaintycontrol = NULL) {
  started <- proc.time()[["elapsed"]]
  seconds <- function() proc.time()[["elapsed"]] - started
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  nsubjects <- length(fit$model_spec$subject_starts)
  module <- .ctJuliaModule(fit$model_spec$project)
  get <- .ctJuliaGet
  # The whole correction's progress line: the screen, each round and the final
  # Hessian all report through this one sink, at the same default verbosity as
  # the optimiser's own line, rather than only when `verbose > 0` -- which is
  # how a 22-minute correction (screen plus one round 4.8 min, then the
  # Hessian 18.8 min) printed nothing at all by default. One sink for the
  # whole function, not one per phase, so a shorter later phrase does not
  # leave the tail of a longer earlier one behind it when a line is
  # overwritten in place.
  sink <- .ctBackendProgressSink(verbose)
  # Unconditional and registered before any early return below: whatever the
  # last thing printed through `sink` was, this closes it so a warning or the
  # ordinary R prompt does not land on the same line. Safe to call with
  # nothing open -- `.ctProgressSink()`'s "break" is then a no-op.
  if (!is.null(sink)) on.exit(sink("", "break"), add = TRUE)
  failed <- function(phrase) {
    if (!is.null(sink)) sink("", "break")
    warning("Laplace continuation skipped: ", phrase, ". The uncorrected fit is ",
      "returned; see fit$laplace$correction.", call. = FALSE)
    fit$laplace$correction <- list(method = "quadrature", status = "failed",
      applied = FALSE, message = phrase, nodes = as.integer(control$nodes))
    fit
  }
  chunks <- suppressWarnings(as.integer(fit$optim$chunks)[1L])
  if (is.na(chunks) || chunks < 1L) chunks <- max(1L, as.integer(cores)[1L])
  previous <- .ctBackendSetMaxChunks(chunks)
  on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  # `sink` is non-NULL exactly when `.ctBackendReporting(verbose)` is, which
  # `verbose > 0` alone already satisfies -- so there is no separate case left
  # for a bare `verbose > 0` to reach for a plainer `message()` here.
  if (!is.null(sink)) {
    sink(sprintf("Laplace continuation screen (%d nodes) | %8s",
      control$nodes, .ctDuration(seconds())), "update")
  }
  # Units whose Laplace terms jumped in the fit (`fit$laplace$jump`) are
  # flagged whatever the screen finds, on draws, with the points the line
  # search saw as components of their proposals; the rounds then climb from
  # the jump point, which is the fit's estimate.
  jump <- fit$laplace$jump
  # Where every unit's integrand is Gaussian in its effects
  # (`laplace$quadratic`), each Laplace term is that unit's exact integral, at
  # any width: nothing to score and nothing to move. The screen's product
  # rules used to run here for nothing, and a unit wider than `maxdim` was
  # reported as left uncorrected.
  if (is.null(jump) && isTRUE(fit$model_spec$laplace$quadratic)) {
    fit$laplace$correction <- list(method = "quadrature", status = "exact",
      applied = FALSE, nodes = as.integer(control$nodes),
      tolerance = as.numeric(control$tolerance), screen = 0, gap = 0,
      laplace_estimate = est, loglik_laplace = as.numeric(fit$estimate$loglik),
      flagged = 0L, wide = 0L, maxdim = as.integer(control$maxdim),
      error = list(nats = 0, se = 0, draws = 0, nodes = 0, residual = 0,
        units = NULL, unscored = 0L),
      seconds = c(total = seconds()))
    return(fit)
  }
  # The shift check first, on every unit whatever its width: where correcting
  # Laplace would move no estimate by `material` standard errors, resolved
  # past the check's own noise, the fit is left alone; otherwise the units
  # carrying the move are flagged on draws as jump units are, and the rounds
  # climb as before. Not at a jump, whose units are named already and whose
  # point is no optimum to take a step from.
  check <- if (is.null(jump) && !isFALSE(control$check)) {
    if (!is.null(sink)) sink(sprintf("Laplace check, every unit on %d draws | %8s",
      as.integer(.ctJuliaOr(control$check_draws, 128L)), .ctDuration(seconds())), "update")
    tryCatch(.ctLaplaceShiftCheck(fit,
      draws = as.integer(.ctJuliaOr(control$check_draws, 128L)),
      seeds = seq_len(as.integer(.ctJuliaOr(control$check_seeds, 4L))),
      bar = as.numeric(control$material), control = control),
      error = function(e) NULL)
  }
  checkrecord <- .ctLaplaceShiftRecord(check, fit)
  if (identical(check$status, "passed") ||
      (identical(check$status, "unresolved") && !length(check$carrying))) {
    # Left alone to the bit, as a passed screen leaves it: "exact" where the
    # corrected terms equal Laplace's to the screen's tolerance, "checked"
    # where they differ without moving the estimate. Every unit was looked
    # at, so none is unscored, and the estimate's distance from the
    # exact-marginal optimum is the check's resolved step, not the bound a
    # value gap of that size would give.
    gap <- abs(check$gap)
    if (!is.finite(gap) || gap < 1e-8) gap <- 0
    record <- list(method = "quadrature",
      status = if (gap <= control$tolerance) "exact" else "checked",
      applied = FALSE, nodes = as.integer(control$nodes),
      tolerance = as.numeric(control$tolerance), screen = gap, gap = check$gap,
      laplace_estimate = est, loglik_laplace = as.numeric(fit$estimate$loglik),
      flagged = 0L, units = as.integer(check$units), wide = 0L,
      maxdim = as.integer(control$maxdim), check = checkrecord)
    bound <- if (gap > 0) max(abs(check$shift) + 2 * check$noise) else 0
    record$error <- list(nats = gap, se = bound, draws = check$gap_se, nodes = 0,
      residual = 0, units = NULL, unscored = 0L)
    record$seconds <- c(check = seconds(), total = seconds())
    fit$laplace$correction <- record
    return(fit)
  }
  # Past a quarter of the units on draws a gradient costs more than the Laplace
  # route saves and sampling is the better route (`.ctLaplaceDrawsShareWarn()`).
  # Only units too wide for a product rule go on draws, and the check spreads
  # the step over most units of a wide model -- 97 to 99 of 100 at twenty
  # effects a unit -- where the rounds on draws would run for hours, so it
  # says so rather than starting them. Narrow units' product rules are cheap.
  drawn <- check$carrying[check$width[check$carrying] > as.integer(control$maxdim)]
  if (!is.null(check) && length(drawn) > 0.25 * check$units) {
    warning(sprintf(paste0("Laplace's error moves the estimates by up to %.1f ",
      "standard errors here (%s), carried by %d of %d units too wide for a ",
      "product rule: correcting it by draws costs more than sampling the ",
      "joint posterior, which ctFitUncertainty(fit, 'sample', control = ",
      "list(intoverpop = 'none')) does. See fit$laplace$correction$check."),
      check$max_shift, checkrecord$worst, length(drawn), check$units),
      call. = FALSE)
    fit$laplace$correction <- list(method = "quadrature", status = "flagged",
      applied = FALSE, laplace_estimate = est,
      loglik_laplace = as.numeric(fit$estimate$loglik),
      units = as.integer(check$units), check = checkrecord,
      seconds = c(check = seconds(), total = seconds()))
    return(fit)
  }
  cont <- try(do.call(module$ctsem_laplace_continuation,
    c(list(.ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(est)),
      .ctLaplaceContinuationArgs(control, force = union(jump$units, check$carrying),
        points = jump$points))),
    silent = TRUE)
  # The engine's reason, which is what tells an allocation it refused from a
  # point the model cannot evaluate.
  if (inherits(cont, "try-error")) return(failed(paste0(
    "the quadrature could not be evaluated (", .ctLaplaceTryReason(cont), ")")))
  start <- get(module$ctsem_laplace_continuation_info(cont))
  .ctLaplaceDrawsShareWarn(intersect(which(as.logical(start$importance)),
    as.integer(start$flagged)[seq_len(start$nflagged)]), as.integer(start$nunits)[1L],
    FALSE)
  # `loglik_laplace` is the Laplace objective's log likelihood at the point the
  # correction starts from, from the placement's own Laplace terms.
  record <- list(method = "quadrature", status = "exact", applied = FALSE,
    nodes = as.integer(start$nodes), tolerance = as.numeric(start$tolerance),
    screen = as.numeric(start$screen), gap = as.numeric(start$gap),
    laplace_estimate = est, loglik_laplace = sum(as.numeric(start$laplace_units)),
    flagged = as.integer(start$nflagged), units = as.integer(start$nunits),
    rule_failures = as.integer(start$rule_failures),
    wide = as.integer(start$nwide), maxdim = as.integer(start$maxdim))
  record$check <- checkrecord
  # What the objective is known to at the screen; the paths that place again
  # replace it with the final placement's.
  record$error <- .ctLaplaceIntegrationError(start, unscored = record$wide)
  # Wall seconds by stage as each is passed -- the screen, the rounds, the
  # continuation's Hessian, its certification, the uncertainty redrawn from it
  # -- and `total` when the record is written.
  stages <- c(screen = seconds())
  timed <- function() c(stages, total = seconds())
  since <- function(from) proc.time()[["elapsed"]] - from
  # Units wider than `maxdim` are not scored and keep the Laplace term (see
  # laplace_continuation.jl for the measured cost). Said in one line, since
  # the correction then covers less than the fit, or none of it -- unless the
  # shift check ran, which looked at them all.
  if (record$wide > 0L && is.null(check)) {
    if (!is.null(sink)) sink("", "break")
    message(sprintf(paste0("Quadrature correction: %d of %d %s kept the ",
      "Laplace term, having more than %d random effects."), record$wide,
      record$units, if (record$units == nsubjects) "subjects" else "groups",
      record$maxdim))
  }
  if (record$wide == record$units) {
    record$status <- "too_wide"
    record$seconds <- timed()
    fit$laplace$correction <- record
    return(fit)
  }
  if (!is.finite(record$screen)) {
    return(failed("the quadrature was not finite at the estimate"))
  }
  # Passed, as the step correction's screen passes: Laplace is exact here to
  # the tolerance, and the fit is left alone to the bit. Not at a jump: the
  # screen says the hybrid agrees with Laplace at the point, which is no
  # optimum of either (round 9's jump point passed the screen with a gradient
  # of 2.5e5), and the rounds are what climb from it. Nor where the shift
  # check found a move to make: the screen reads values, the check the move.
  if (is.null(jump) && is.null(check) && record$screen <= record$tolerance) {
    record$seconds <- timed()
    fit$laplace$correction <- record
    return(fit)
  }
  # The quadrature log likelihood wherever the screen found a gap, as the step
  # correction reports it, whether or not anything moves.
  report_quadrature <- function(fit, info) {
    units <- as.numeric(info$quadrature_units)
    if (!all(is.finite(units))) return(fit)
    fit$estimate$loglik_laplace <- record$loglik_laplace
    fit$estimate$loglik <- sum(units)
    fit$estimate$logposterior <- as.numeric(info$quadrature)
    subjects <- as.numeric(info$quadrature_subjects)
    if (length(subjects) == length(fit$estimate$subject_loglik) &&
        all(is.finite(subjects))) fit$estimate$subject_loglik <- subjects
    fit$estimate$loglik_method <- "quadrature"
    fit
  }
  # The fit's curvature at the estimate is the rounds' first metric and
  # basis. Asked for no rounds (`rounds = 0`: the rules placed at the
  # estimate and the uncertainty drawn there, for a fit already corrected or
  # one a caller wants left where it is) the continuation needs neither, and
  # its own Hessian at the end is the one the uncertainty is drawn on.
  stay <- as.integer(.ctJuliaOr(control$rounds, 0L)) == 0L
  hessian <- .ctBackendStoredHessian(fit, est)
  flat <- fit$uncertainty$details$flatdirections$vectors
  # A fit that ended on a jump was not certified and has no Hessian: the
  # rounds' basis is the continuation's own curvature at the jump point,
  # which is continuous there where the Laplace objective is not. Its Krylov
  # model (`ctsem_krylov_hessian`: the scores' curvature corrected by
  # Hessian-vector products), which serves to whiten and bound the rounds'
  # steps; the exact one costs `2 npar` gradients, and the continuation forms
  # its exact Hessian where the rounds end anyway. The exact one where no
  # model can be made.
  if (!is.null(jump) && (is.null(hessian) || !all(is.finite(hessian)))) {
    if (!is.null(sink)) {
      sink(sprintf("Laplace continuation curvature at the jump | %8s",
        .ctDuration(seconds())), "update")
    }
    square <- function(h) {
      h <- if (inherits(h, "try-error")) numeric(0) else as.numeric(h)
      if (length(h) == npar^2 && all(is.finite(h))) matrix(h, npar, npar) else NULL
    }
    hessian <- square(try(.ctBackendJuliaValue(module$ctsem_krylov_hessian(cont,
      .ctJuliaNumericVector(est))), silent = TRUE))
    if (is.null(hessian)) {
      hessian <- square(try(.ctBackendJuliaValue(
        module$ctsem_laplace_continuation_hessian(cont, .ctJuliaNumericVector(est))),
        silent = TRUE))
    }
    flat <- NULL
  }
  if (is.null(hessian) || !all(is.finite(hessian))) {
    if (!stay) {
      record$status <- "no_hessian"
      record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
      fit <- report_quadrature(fit, start)
      fit$laplace$correction <- record
      return(fit)
    }
    hessian <- NULL
  }
  basis <- if (is.null(hessian)) list(basis = diag(npar), kept = npar) else
    .ctLaplaceContinueBasis(hessian, flat, rtol = control$rtol)
  if (is.null(basis) || basis$kept < 1L) {
    if (!stay) {
      record$status <- "no_directions"
      record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
      fit <- report_quadrature(fit, start)
      fit$laplace$correction <- record
      return(fit)
    }
    basis <- list(basis = diag(npar), kept = npar)
  }
  tol <- .ctBackendGapTolerance(fit)
  # What the continuation cost, in the engine's own counts -- whole-objective
  # values and gradients, placements of the rule, and the member likelihood
  # values and reverse sweeps they took -- and its seconds in each kind of
  # call (`engine_seconds`, the screen's placement among the placements),
  # beside the wall seconds by stage.
  account <- function(record) {
    after <- get(module$ctsem_laplace_continuation_info(cont))
    record$error <- .ctLaplaceIntegrationError(after, unscored = record$wide)
    record$evaluations <- c(values = as.integer(after$value_calls),
      gradients = as.integer(after$gradient_calls),
      placements = as.integer(after$recentres),
      member_values = as.numeric(after$member_values),
      member_sweeps = as.numeric(after$member_sweeps),
      refused = as.integer(after$refused))
    record$engine_seconds <- c(values = as.numeric(after$seconds_values),
      gradients = as.numeric(after$seconds_gradients),
      placements = as.numeric(after$seconds_placements),
      hessian = as.numeric(after$seconds_hessian))
    record$seconds <- timed()
    record
  }
  rounds_started <- proc.time()[["elapsed"]]
  run <- try(.ctLaplaceContinueRun(module, cont, est, basis$basis, tol,
    control = control, verbose = verbose, sink = sink,
    skip_gain = as.numeric(.ctJuliaOr(control$skip_gain, 0))), silent = TRUE)
  stages["rounds"] <- since(rounds_started)
  if (inherits(run, "try-error")) {
    return(failed(paste0("a round could not be evaluated (",
      .ctLaplaceTryReason(run), ")")))
  }
  record$predicted_gain <- as.numeric(run$predicted)
  record$skip_gain <- as.numeric(.ctJuliaOr(control$skip_gain, 0))
  record$stop_gain <- as.numeric(run$bar)
  record$resolution <- as.numeric(run$resolution)
  # Skipped: the fixed-node model promised less than `skip_gain` from the
  # Laplace optimum, so no round ran. The estimate stays, and the log
  # likelihood reported is the quadrature one there, as where a round ran and
  # gained nothing. A run asked to stay goes on to draw there.
  if (!stay && identical(run$status, "skipped")) {
    record$status <- "skipped"
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    record$logposterior_quadrature <- as.numeric(start$quadrature)
    record$gap_reported <- as.numeric(start$quadrature) - as.numeric(start$laplace)
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- account(record)
    return(fit)
  }
  info <- run$info
  change <- as.numeric(info$quadrature) - as.numeric(start$quadrature)
  limit <- max(control$guard, control$guard_per_subject * nsubjects)
  record$rounds <- run$rounds
  record$attempts <- run$attempts
  record$radius <- run$radius
  record$stationarity <- run$residual
  record$trace <- run$trace
  record$continuation <- run$status
  record$frozen_directions <- as.integer(basis$frozen)
  record$flagged <- as.integer(info$nflagged)
  record$flagged_units <- as.integer(info$flagged)[seq_len(info$nflagged)]
  record$bar <- as.numeric(info$bar)
  record$residual <- as.numeric(info$residual)
  record$soft_blocks <- as.integer(info$soft_blocks)
  record$guard <- list(limit = limit, change = change,
    fired = !is.finite(change) || abs(change) > limit)
  # The run moved nothing: no round was kept, the rounds ended lower on the
  # objective than they began, or the guard reverts it. The quadrature log
  # likelihood at the Laplace optimum is reported, as the step correction
  # reports it when no step raised the objective. Not a run asked to stay,
  # which moves nothing by design and goes on to its Hessian and draws there;
  # returning here, `likelihood = 'quadrature'` on a corrected fit never drew
  # on the hybrid at all.
  worse <- is.finite(change) && change < -as.numeric(control$value_tol)
  record$ended_lower <- worse
  moved <- run$rounds > 0L && any(run$x != est) && !worse
  if (isTRUE(record$guard$fired) || (!moved && !stay)) {
    if (isTRUE(record$guard$fired)) {
      if (!is.null(sink)) sink("", "break")
      warning(sprintf(paste0("The Laplace continuation moved the quadrature ",
        "objective by %.3g nats, which no approximation error of this size ",
        "explains; the fit keeps the Laplace optimum, and the rejected point is ",
        "fit$laplace$correction$rejected_estimate."), change), call. = FALSE)
      record$rejected_estimate <- run$x
      record$status <- "reverted"
    } else {
      record$status <- "no_gain"
    }
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    record$logposterior_quadrature <- as.numeric(start$quadrature)
    record$gap_reported <- as.numeric(start$quadrature) - as.numeric(start$laplace)
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- account(record)
    return(fit)
  }
  x <- run$x
  # The continuation's own curvature only where it reached its fixed point.
  # Stopped short, the estimate is not stationary on the objective the Hessian
  # is of, and on gated-gaps A14 that Hessian had two directions of positive
  # curvature: ten parameters would have reported no interval, and the fit
  # `converged: FALSE`, where the Laplace curvature at the fit's optimum gives
  # them all one. So a continuation that stopped short keeps the fit's
  # uncertainty and recentres its draws, as the step correction does, and
  # takes no Hessian at all: on wide blocks that is a large share of a
  # correction's cost, and `stationarity` in the record already says how far
  # from stationary the estimate is. A fit that ended on a jump has no other
  # curvature, so it takes this one however the rounds ended: its rounds
  # stall on the draws' re-placement differences (round 9: 0.017 nats of
  # residual against draws' errors of 0.04), and what the fit reports is
  # then the hybrid's certification at that point, gap and all.
  # A run of no rounds reached the point it was asked to stay at.
  reached <- run$status %in% c("converged", "resolution") ||
    (stay && isTRUE(as.integer(run$attempts) == 0L))
  takes_hessian <- reached || !is.null(jump)
  hc <- NULL
  if (takes_hessian) {
    # By forward differences this Hessian was 18.8 minutes of a 22-minute
    # correction on the ordinal fixture that found the silence; it is exact
    # now (laplace_continuation.jl) and can still take minutes. A start and a
    # done line around the call, since the engine reports no count from
    # inside it; the sink reports at every verbosity that reports at all.
    if (!is.null(sink)) {
      sink(sprintf("Laplace continuation hessian | %8s",
        .ctDuration(seconds())), "update")
    }
    hessian_started <- proc.time()[["elapsed"]]
    hc <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_laplace_continuation_hessian(cont, .ctJuliaNumericVector(x)))),
      npar, npar), silent = TRUE)
    stages["hessian"] <- since(hessian_started)
    if (!is.null(sink)) {
      sink(sprintf("Laplace continuation hessian done in %s",
        .ctDuration(stages[["hessian"]])), "done")
    }
  }
  final <- get(module$ctsem_laplace_continuation_evaluate(cont,
    .ctJuliaNumericVector(x), gradient = TRUE))
  # Where the continuation's own curvature says more is left than the rounds'
  # stopping rule saw, the step it predicts is taken before anything is
  # reported (`.ctLaplaceContinueFinish()`). `hessian_at` keeps the point `hc`
  # describes, which the estimate may then be a few thousandths of a standard
  # error from.
  hessian_at <- x
  if (!stay && takes_hessian && is.matrix(hc) && all(is.finite(hc))) {
    finish <- .ctLaplaceContinueFinish(module, cont, x, hc, final,
      value = run$value, bar = run$bar, value_tol = control$value_tol,
      reuse = .ctBackendHessianReuse(.ctBackendFitOptimcontrol(fit)))
    record$newton <- finish[c("steps", "gap_before", "gap_after")]
    if (finish$steps > 0L) {
      x <- finish$x
      final <- finish$final
      info <- finish$info
      record$flagged <- as.integer(info$nflagged)
      record$flagged_units <- as.integer(info$flagged)[seq_len(info$nflagged)]
      record$bar <- as.numeric(info$bar)
      record$residual <- as.numeric(info$residual)
    }
  }
  # A run asked to stay placed its rules at the estimate and moved nothing.
  record$status <- if (stay) "placed" else "continued"
  record$applied <- !stay
  se <- as.numeric(fit$estimate$se)
  delta <- x - est
  delta_se <- if (length(se) == npar) ifelse(se > 0, delta / se, NA_real_) else
    rep(NA_real_, npar)
  names(delta) <- names(delta_se) <- names(fit$estimate$se)
  record$delta <- delta
  record$delta_se <- delta_se
  # No standard errors yet on a fit that ended on a jump (its uncertainty
  # stage waits for this Hessian), so nothing to measure the move against.
  record$material <- any(is.finite(delta_se)) &&
    isTRUE(max(abs(delta_se), na.rm = TRUE) >= control$material)
  record$loglik_quadrature <- sum(as.numeric(info$quadrature_units))
  record$logposterior_quadrature <- as.numeric(info$quadrature)
  record$gap_reported <- as.numeric(info$quadrature) - as.numeric(info$laplace)
  # The objective the estimate maximises: the hybrid, whose unflagged units are
  # Laplace. It differs from the quadrature value above by at most `residual`.
  record$objective <- as.numeric(final$value)
  record$error <- .ctLaplaceIntegrationError(info, unscored = record$wide)
  fit$estimate$raw <- x
  fit <- report_quadrature(fit, info)
  usable <- takes_hessian && is.matrix(hc) && all(is.finite(hc))
  # The uncertainty drawn on the continuation's density: the method asked
  # for (`likelihood = 'quadrature'`, or a fit's), else the fit's own, with
  # the control asked for, else the fit's.
  method <- .ctJuliaOr(method, .ctJuliaOr(fit$uncertainty$settings$method, "hessian"))
  # Present even when empty: `$` partial-matches, and without it
  # `correction$hessian` would answer with `hessian_at`.
  record["hessian"] <- list(NULL)
  if (usable) {
    # The certification of the estimate on the objective it maximises: its
    # gradient against its own Hessian, with the flat-direction probe walking
    # the hybrid.
    certification_started <- proc.time()[["elapsed"]]
    gradient <- as.numeric(final$gradient)
    gap <- .ctBackendOptimGap(hc, gradient)
    probe <- NULL
    # Only along a direction the curvature does not trust. With every one
    # trusted the residual is the rounding left by projecting onto a complete
    # eigenbasis, and the probe walked that: three values of the hybrid, a
    # fifteenth of ord4's correction, and a gain of 0 on every converged
    # continuation of the bench (controbust-c2a67356). A negative direction
    # is walked by the saddle ladder first, and is a saddle unless that found
    # it flat, when the probe covers it as it does a flat one.
    saddle <- NULL
    if (isTRUE(gap$ok) && isTRUE(gap$nnegative > 0L)) {
      saddle <- .ctBackendSaddleProbe(module, cont, x, hc, ascent = gradient)
    }
    if (isTRUE(gap$ok) && isTRUE(gap$nflat > 0L || gap$nnegative > 0L) &&
        isTRUE(gap$residual_norm > 0)) {
      out <- try(get(module$ctsem_flat_probe(cont, .ctJuliaNumericVector(x),
        .ctJuliaNumericVector(gap$residual))), silent = TRUE)
      if (!inherits(out, "try-error")) probe <- .ctBackendProbeFields(out, npar)
    }
    parnames <- try(.ctBackendRawParameterNames(fit, npar), silent = TRUE)
    if (inherits(parnames, "try-error")) parnames <- NULL
    # Certified to the bar the rounds stopped on, not the Laplace fit's: the
    # quadrature rule is not resolved to 1e-6 nats, and a continuation that
    # met its own stopping rule is not a failure to converge. Nor is one
    # within its draws' resolution of the fixed point, which no round can
    # read past; the gap the record keeps says how far it is.
    certification <- .ctBackendCertificationRecord(gap, probe,
      tolerance = max(run$bar, run$resolution),
      saturated = isTRUE(fit$optim$saturated), overshot = FALSE,
      parnames = parnames, saddle = saddle)
    record$certification <- certification
    record$hessian <- hc
    # A fit that ended on a jump has no other certification: this is its
    # verdict, and the gradient it reports is the hybrid's here, not the
    # Laplace objective's at the jump.
    if (!is.null(jump)) {
      fit$optim$converged <- identical(certification$status, "certified")
      fit$optim$gradient_norm <- max(abs(gradient))
      fit$optim$stalled <- FALSE
    }
    stages["certification"] <- since(certification_started)
  }
  own <- usable && method %in% c("hessian", "checked", "is")
  record$hessian_at <- if (own) "estimate" else "laplace_estimate"
  if (own) {
    uncertainty_started <- proc.time()[["elapsed"]]
    settings <- fit$uncertainty$settings
    finishsamples <- .ctJuliaOr(settings$finishsamples,
      if (!is.null(fit$estimate$rawposterior)) nrow(fit$estimate$rawposterior) else 1000L)
    fit$uncertainty <- list(hessian = hc, evaluated_at = hessian_at,
      certification = certification)
    fit <- .ctBackendUncertainty(fit, uncertainty = method,
      draws = .ctResolveDraws(method, "auto"),
      finishsamples = finishsamples, cores = chunks,
      control = .ctJuliaOr(uncertaintycontrol, .ctJuliaOr(settings$control, list())),
      verbose = verbose, lpg = .ctLaplaceContinueLpg(module, cont))
    record$draws <- "redrawn"
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("draws about the continuation's estimate, by '",
        method, "' on its quadrature objective there"),
      nodes = record$nodes,
      # What the density the draws came from is known to, in nats.
      error = record$error)
    stages["uncertainty"] <- since(uncertainty_started)
  } else {
    # Recentred, as the step correction does, when the continuation stopped
    # short of its fixed point, has no Hessian to offer, or the fit's
    # uncertainty did not come from one.
    post <- fit$estimate$rawposterior
    if (!is.null(post) && ncol(post) == npar) {
      fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
    }
    record$draws <- "recentred"
    if (!is.null(fit$uncertainty)) {
      fit$uncertainty$details$laplace_correction <- list(
        draws = paste0("recentred on the continuation's estimate; the ",
          "covariance is the fit's own at fit$uncertainty$evaluated_at"),
        nodes = record$nodes)
    }
  }
  fit$laplace$correction <- account(record)
  fit
}
