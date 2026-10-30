# The random-effects check and the refined fit (random effects by
# quadrature), both `.ctLaplaceContinue()`: the check scores every unit's
# Laplace error at the estimate (the default, `optimcontrol$checks`); the
# refined fit (`ctFit(intoverpop = 'quadrature')`, or `ctFitUncertainty(fit, intoverpop =
# 'quadrature')` after the fit) climbs the quadrature hybrid from the estimate
# with the units the check names integrated, and draws the fit's uncertainty
# on that density. A request to sample against the quadrature --
# `ctFitUncertainty(fit, 'is' or 'checked', intoverpop = 'quadrature')`,
# `.ctBackendUncertaintyQuadrature()` -- is the continuation with the method
# asked for. The separate importance sampler against the adaptive
# Gauss-Hermite posterior at every unit (`.ctLaplaceCorrect()`, to
# 2026-10-23) paid `nodes^k` member evaluations a unit a draw where the hybrid
# pays them only for the units the screen flags, moved the estimate by one
# Newton step where the continuation climbs, and reported no error for the
# density it sampled.

# What a fit does about Laplace's error -- FALSE, "check" (the random-effects
# check) or "quadrature" (the refined fit) -- from `optimcontrol$checks`
# and `ctFit(intoverpop = 'quadrature')` (`optimcontrol$.refit`, set by
# ctFit), refusing a refinement that cannot apply. The default is the check
# alone, which warns where Laplace's error moves an estimate and names the
# refined fit, rather than running it (Charles, 2026-10-27): it adds seconds
# to hours to a fit. Wherever both were
# measured it cost no more than sampling the joint posterior -- 1/25 to 1/200
# of it on one level, as much on nested units (gN3: 244 s against 242 s),
# under a quarter on wide units by draws (k5: 2.8 h against over 13;
# scripts/lapcheck-2026-10/rounds_cost.R, sample_cost.R) -- but whether to
# pay for either is the user's call. A refined fit runs the check first, whatever
# `checks` says.
.ctLaplaceCorrectResolve <- function(optimcontrol, intoverpop, optimize,
  intoverstates) {
  refit <- isTRUE(optimcontrol$.refit)
  check <- "randomeffects" %in% .ctUncertaintyChecks(.ctJuliaOr(optimcontrol$checks, TRUE))
  if (!refit && !check) return(FALSE)
  method <- if (refit) "quadrature" else "check"
  why <- if (!identical(as.character(intoverpop)[1L], "laplace")) {
    paste0("applies to intoverpop='laplace' only; with intoverpop='",
      intoverpop, "' there is no Laplace term")
  } else if (!isTRUE(optimize)) {
    paste0("concerns an optimised estimate, and a sampled fit ",
      "(optimize=FALSE) has none")
  } else if (!isTRUE(intoverstates)) {
    paste0("concerns the marginal Laplace estimate, which an ",
      "intoverstates=FALSE fit does not make")
  } else if (isTRUE(optimcontrol$estonly)) {
    "steps against the fit's Hessian, which optimcontrol$estonly skips"
  } else NULL
  if (is.null(why)) return(method)
  if (refit) stop("intoverpop = 'quadrature' ", why, ".", call. = FALSE)
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
# for a fit whose refined fit did not draw it on its own objective (`ctFit()`).
# Where the refined fit moved the estimate with no Hessian of its own that the
# method reads (`.ctLaplacePlacedRefine()`, `draws = "recentred"`) the stored
# curvature is the Laplace optimum's, so the method runs there and its draws
# move with the estimate, as the refined fit moves the normal's; otherwise it
# runs at the estimate.
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

# The refined fit (`intoverpop = 'quadrature'`) climbs the quadrature
# objective from the Laplace optimum, its flagged units' rules placed at every
# evaluation (`laplace_placed.jl`, `.ctLaplacePlacedRefine()`). bigIRT's
# `laplaceRefine` was the first model (../bigIRT/optimisation-review-2026-09.md,
# R3.4 and R3.5).
#
#   screen   the rule at the Laplace optimum against Laplace, unit by unit.
#            A fit that passes is left untouched to the bit.
#   flags    only the units the screen flags carry rules: largest gap first,
#            until what is left carries at most the screen's tolerance. Every
#            other unit keeps its Laplace term and its exact gradient.
#   climb    `ctsem_optimize`, L-BFGS and the Newton finish, on one objective
#            whose gradient is its value's, the placement's own derivative
#            included: its maximum is that of the quadrature value with the
#            rules placed at every point. Until 2026-10-30 rounds climbed the
#            rules frozen where they were placed and placed them again between
#            rounds, towards the fixed point of the frozen rules' score, which
#            is not that maximum; on the race's seven fixtures five ran out of
#            rounds with no Hessian, where the placed climb certified all seven
#            (review/PLACED-rules-race-2026-10-30.md). At 5 nodes neither
#            point was resolved better than the other on a 9-node reference.
#   draws    a unit on importance draws keeps them for the whole run, placed
#            again only where their effective size has collapsed, the run's
#            L-BFGS memory then dropped (`ctsem_laplace_placed_optimize`).
#   no worse the estimate is only reported if the objective there is not
#            below its value at the Laplace optimum by more than `value_tol`;
#            otherwise the fit keeps the Laplace optimum (status `no_gain`),
#            with the quadrature log likelihood there reported.
#   guard    a run that moves the objective by more than max(50, N/2) nats is
#            reverted to the Laplace optimum with a warning, keeping the
#            rejected point, as bigIRT does.
#   Hessian  of the same objective at the estimate, by differences of its
#            exact gradient over each column's readers; the certification, the
#            covariance and the draws come from it, the draws on that objective.
#   stay     (`stay = TRUE`) places the rules at the estimate and moves
#            nothing: a fit already refined, or `correctEstimate = FALSE`.
#
# What moves when it applies: `fit$estimate$raw` is the climb's estimate;
# `loglik`, `logposterior` and `subject_loglik` the quadrature values there
# (every unit by its rule placed there, `ctsem_laplace_placed_quadrature`),
# with `loglik_laplace` kept; `cov`, `se` and the draws from the objective's
# Hessian, which is `fit$uncertainty$hessian`, and whose certification is the
# fit's (`hessian_at` in the record says which). `fit$laplace$correction` says
# what ran.
#
# Provenance of the constants. `nodes` (5): enough to locate the maximum in
# the cases tested, 9 to settle the value (quadrature.jl's docstring measured
# the gap on one 40-subject dataset). `tolerance` (0.01 nats, summed
# absolutely over units): the screen's pass mark, in objective units; the
# linear fixture screens at 3e-14, so any model with a nonlinear effect clears
# it by orders of magnitude. `material` (0.1 se): when print() says the
# estimate moved. Which rule each unit takes, product Gauss-Hermite or
# lattice-t, and at how many points, is the engine's choice by counted cost
# for the accuracy, at any width (`_continuation_refine!`, point_sets.jl); the
# climb places it at every evaluation (`laplace_placed.jl`). `maxiter` (1000)
# caps the climb's L-BFGS, as a fit's own does.
# `guard` and `guard_per_subject` are bigIRT's max(50, N/2).
# `value_tol` (1e-3 nats): a change in the objective below it is not one to
# act on either way.
# `rtol` (1e-8) is the identifiability report's: a direction the Laplace
# curvature does not identify is held where the fit left it, in the shift
# check's step. `draws_se` (0.05 nats) is the error every flagged unit's
# rule is refined to. A unit no rule reaches it on is scored by frozen draws
# instead, per block of its tree when it nests, `draws` (256) the first batch,
# doubling to that standard error or `draws_max` (16384) draws, the blocks
# beneath another each to its share of that, explained in
# laplace_continuation.jl. `draws_df` (4) is the degrees of freedom of the proposal's components, a
# Student t whose tails cover what a Gaussian's at the mode missed (0.6 to 1.0
# nats over gA1's 40 units with the replicates' spread saying 0.1), `Inf` for
# Gaussian; only a unit's top block takes it.
# `check` (TRUE) runs the shift check first (`.ctLaplaceShiftCheck()`): every
# unit, whatever its width, on `check_draws` (128) frozen draws at each of
# `check_seeds` (4) seeds. Where the step correcting Laplace would take is
# resolved under `material` the climb does not run; otherwise the units
# carrying it are forced into the climb, each on the rule its cost chooses --
# or, under `.ctLaplaceContinue(correct =
# FALSE)`, the default fit's, named in a warning with the calls that correct
# them, and kept in the record (`check$carriers`) for the post-hoc
# correction to force. 128 draws resolved the step to 0.05 se
# on the bench's narrow fixtures and 32 did not (0.17-0.45 se of noise where
# nothing moved); four seeds give the noise, which at ten effects a unit with
# a weakly identified covariance was larger than the step itself
# (scripts/lapcheck-2026-10/highdim_check.R).
.ctLaplaceContinueDefaults <- list(nodes = 5L, tolerance = 0.01,
  maxiter = 1000L, material = 0.1, guard = 50, guard_per_subject = 0.5,
  rtol = 1e-8, value_tol = 1e-3,
  draws = 256L, draws_se = 0.05, draws_max = 16384L, draws_df = 4,
  check = TRUE, check_draws = 128L, check_seeds = 4L, stay = FALSE)

# How far the hybrid's objective is from the exact marginal, from the engine's
# info at a placement (`ctsem_laplace_continuation_info`): `nats`, the frozen
# rules' standard errors in quadrature plus the product rules' node errors
# plus what the screen leaves on Laplace; `se`, how far the estimate can be
# from the exact-marginal optimum in standard errors, by the certification's
# own arithmetic (a displacement costing `e` nats on the quadratic is at most
# `sqrt(2 e)` long); the three parts; and each unit's own error (`NaN` where
# its term is Laplace's).
.ctLaplaceIntegrationError <- function(info) {
  nats <- as.numeric(info$error)[1L]
  # Rounding is not an error: on a linear Gaussian fit the screen's gaps are
  # 1e-14 nats of arithmetic, and reported as such they read as a finding.
  if (is.finite(nats) && nats < 1e-8) nats <- 0
  list(nats = nats,
    se = if (is.finite(nats) && nats >= 0) sqrt(2 * nats) else NA_real_,
    draws = as.numeric(info$error_draws)[1L],
    nodes = as.numeric(info$error_nodes)[1L],
    residual = as.numeric(info$error_residual)[1L],
    units = as.numeric(info$centre_error))
}

# The engine arguments of a continuation built under `control`: the screen's
# node count and tolerance, the frozen importance rule's for units no rule
# converges on (`draws`, `draws_se`, `draws_max`; laplace_continuation.jl's
# module docstring), the units the caller names to flag whatever their gap (`force`),
# points it knows for them (`points`, a list of `units` and a `values` matrix
# with a column each), and whether every unit is scored (`screen`). Optional
# entries are left out rather than sent empty, which would hang the bridge.
.ctLaplaceContinuationArgs <- function(control, force = NULL, points = NULL,
  screen = TRUE) {
  args <- list(nodes = as.integer(control$nodes),
    tolerance = as.numeric(control$tolerance),
    draws = as.integer(.ctJuliaOr(control$draws, 0L)),
    draws_se = as.numeric(.ctJuliaOr(control$draws_se, 0.05)),
    draws_max = as.integer(.ctJuliaOr(control$draws_max, 16384L)),
    draws_df = as.numeric(.ctJuliaOr(control$draws_df, 4)),
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

# The directions the shift check steps in (`.ctLaplaceShiftCheck()`): the
# Laplace curvature's identified ones, each scaled to one of its standard
# errors.
#
# Left out, and held where the fit left them: the directions the likelihood was
# measured flat along (`flat`, the columns of
# `fit$uncertainty$details$flatdirections$vectors`), and those whose curvature
# is below `rtol` of the largest -- the identifiability report's two rules for
# a direction the data do not identify. A direction of negative curvature is
# kept, whitened by its magnitude: the quadrature objective may rise along it.
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
# gradient takes in the identified directions (`.ctLaplaceContinueBasis()`) --
# how far correcting Laplace would move
# each estimate, in its standard errors. The seeds are replicates, so the step
# is their mean and its noise their spread over sqrt(seeds): at ten effects a
# unit with a weakly identified covariance one seed's step was mostly noise.
# `status` is "passed" where every parameter's step is resolved under `bar`
# (its size plus twice its noise), "flagged" where one is resolved over it,
# and "unresolved" where the noise decides neither. `carrying` are the units
# whose own steps account for the step, largest first, until what is left is
# under `bar`. `gap` is the corrected minus the Laplace log likelihood: on
# each unit where Laplace's term with its O(1/n) correction
# (`ctsem_laplace_expansion`) has the tighter error (`.ctLaplaceExpanded()`),
# that, whose bounds `gap_bound` sums; on the rest, the draws, whose seeds give
# `gap_se`. NULL where the fit has no curvature to step on or fewer than two
# seeds scored.
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
    seeds = .ctJuliaVector(as.integer(seeds)), draws_df = as.numeric(control$draws_df),
    expansion_gate = .ctLaplaceExpansionGate))
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
  expanded <- .ctLaplaceExpanded(r, length(laplace_units), scored)
  expansion <- as.numeric(r$expansion_term) + as.numeric(r$expansion_correction) -
    laplace_units
  gaps <- colSums(checked[!expanded, scored, drop = FALSE]) -
    sum(laplace_units[!expanded])
  unit_gap <- rowMeans(checked[, scored, drop = FALSE]) - laplace_units
  unit_gap[expanded] <- expansion[expanded]
  list(status = status, shift = shift, noise = noise, bar = bar,
    max_shift = max(abs(shift)), worst = which.max(abs(shift)),
    carrying = carrying, gap = sum(expansion[expanded]) + mean(gaps),
    gap_se = stats::sd(gaps) / sqrt(length(gaps)),
    gap_bound = sum(as.numeric(r$expansion_bound)[expanded]),
    expanded = sum(expanded), expanded_units = which(expanded), unit_gap = unit_gap,
    units = length(laplace_units), draws = as.integer(draws), seeds = sum(scored),
    seconds = proc.time()[["elapsed"]] - started)
}

# The smallest eigenvalue of a unit's curvature in its effects, relative to the
# prior's 1, from which Laplace's O(1/n) correction (`ctsem_laplace_expansion`)
# may integrate the unit for the reported log likelihood. At 1 the likelihood
# is concave in every direction at the mode. On 17 fixtures, 2,486 units
# (review/LAPLACE-correction-measure-2026-10-30.md, section 8), every unit at or
# above it was within its bound, `|c4| + c3`, of a converged product rule.
# Below it the evidence is AnomAuth's alone, the only fixture whose units went
# below 0.906: its units with three to five waves exceeded the bound from 0.95
# down, and below about 0.8 the correction made them worse, by up to 8 nats.
.ctLaplaceExpansionGate <- 1

# Which units the expansion integrates: those past the gate whose bound is
# under twice the draws' standard error there -- of the unit's two error
# statements the tighter -- and whose expansion evaluated; none where the
# engine returned none. Past the gate alone was worse: the correction errs
# with one sign, by 1.01-1.5 times Laplace's error, so over many units its
# errors add while the draws' average out, and on 17 fixtures the summed gap
# was 8.5 nats off against the draws' 2.3. With the choice it was 1.6, and no
# fixture's was further off than the draws'.
.ctLaplaceExpanded <- function(r, nunits, scored = NULL) {
  term <- as.numeric(r$expansion_term)
  corr <- as.numeric(r$expansion_correction)
  bound <- as.numeric(r$expansion_bound)
  lowest <- as.numeric(r$expansion_mineig)
  if (length(term) != nunits || length(corr) != nunits || length(bound) != nunits ||
      length(lowest) != nunits) return(logical(nunits))
  se <- matrix(as.numeric(r$check_se), nrow = nunits)
  if (is.null(scored)) scored <- apply(se, 2L, function(v) all(is.finite(v)))
  if (!any(scored)) return(logical(nunits))
  # The draws' standard error for their mean over the seeds that scored.
  drawn_se <- rowMeans(se[, scored, drop = FALSE]) / sqrt(sum(scored))
  is.finite(term) & is.finite(corr) & is.finite(bound) & !is.na(lowest) &
    lowest >= .ctLaplaceExpansionGate & bound < 2 * drawn_se
}

# What a fit keeps of the shift check (`fit$laplace$correction$check`): the
# verdict, the largest step and its noise by parameter name, the Laplace gap,
# how many units carried the step and which (`carriers`, the engine's unit
# indices, which a correction made later forces). NULL where no check ran.
.ctLaplaceShiftRecord <- function(check, fit) {
  if (is.null(check)) return(NULL)
  names <- .ctBackendRawParameterNames(fit, length(check$shift))
  list(status = check$status, max_shift = check$max_shift,
    worst = names[check$worst], noise = check$noise[check$worst],
    shift = stats::setNames(check$shift, names),
    shift_noise = stats::setNames(check$noise, names), bar = check$bar,
    gap = check$gap, gap_se = check$gap_se, gap_bound = check$gap_bound,
    expanded = check$expanded, expanded_units = as.integer(check$expanded_units),
    carrying = length(check$carrying),
    carriers = as.integer(check$carrying),
    draws = check$draws, seeds = check$seeds, seconds = check$seconds)
}

# The fit's log likelihood where the random-effects check integrated every unit and
# left the estimate alone: Laplace's plus the check's gap, the summed
# difference between each unit's integral and its Laplace term -- by Laplace's
# O(1/n) correction where that has the tighter error, within `loglik_bound`,
# and on its draws elsewhere (their seeds' mean, `loglik_se` their standard
# error) -- and each subject's term likewise where units are subjects. Laplace's value, the one the
# optimiser maximised, is kept as `loglik_laplace`, as the refined fit keeps it
# (`loglik_method` says which a fit reports): a model's log likelihood is for
# comparing models, and Laplace's error differs between them -- 30 nats on
# anomS1, 0.2 on cf_binary. Comparisons of where optimisations ended read the
# maximised value (`.ctFitOptimisedLoglik()`).
.ctLaplaceCheckedLoglik <- function(fit, check) {
  gap <- as.numeric(check$gap)
  if (length(gap) != 1L || !is.finite(gap) || gap == 0) return(fit)
  fit$estimate$loglik_laplace <- as.numeric(fit$estimate$loglik)
  fit$estimate$loglik <- fit$estimate$loglik_laplace + gap
  if (!is.null(fit$estimate$logposterior)) {
    fit$estimate$logposterior <- as.numeric(fit$estimate$logposterior) + gap
  }
  fit$estimate$loglik_se <- as.numeric(check$gap_se)
  fit$estimate$loglik_bound <- as.numeric(.ctJuliaOr(check$gap_bound, 0))
  units <- as.numeric(check$unit_gap)
  if (length(units) && length(units) == length(fit$estimate$subject_loglik) &&
      all(is.finite(units))) {
    fit$estimate$subject_loglik <- fit$estimate$subject_loglik + units
  }
  fit$estimate$loglik_method <- "checked"
  fit
}

# How far a checked log likelihood can be off, in words: the bound on the units
# the expansion integrated and the draws' standard error on the rest, each
# where it applies. A fit checked before the expansion has only the second.
.ctLoglikCheckedError <- function(estimate) {
  bound <- as.numeric(.ctJuliaOr(estimate$loglik_bound, 0))
  se <- as.numeric(.ctJuliaOr(estimate$loglik_se, NA_real_))
  parts <- c(if (length(bound) && is.finite(bound) && bound > 0)
      paste0("within ", format(bound, digits = 2), " by the expansion"),
    if (length(se) && is.finite(se) && se > 0) paste0("se ", format(se, digits = 2),
      " by draws"))
  if (!length(parts)) return("exact on every unit")
  paste(parts, collapse = ", ")
}

# The log likelihood the optimiser maximised: Laplace's on a fit that reports
# the check's integral (`.ctLaplaceCheckedLoglik()`), else the one reported.
.ctFitOptimisedLoglik <- function(fit) {
  if (identical(fit$estimate$loglik_method, "checked") &&
      !is.null(fit$estimate$loglik_laplace)) {
    return(as.numeric(fit$estimate$loglik_laplace))
  }
  as.numeric(fit$estimate$loglik)
}

# What a fit the check flagged and nothing corrected says, in its warning and
# in print(): the move, Laplace's error in the log likelihood (`gap`, where it
# is resolved past twice `gap_se`), that the estimates and intervals are
# Laplace's while the log likelihood the fit reports is the integrated one
# (`.ctLaplaceCheckedLoglik()`), and the calls that remove the error -- the
# refined fit after the fit, and sampling the joint posterior, which has no
# Laplace error.
#
# The step is the check's linearisation, right in size to about one standard
# error and no guide past it, in either direction: on the bench it read the
# refined fit's moves of 3.5 and 10 se as 1.1 and 1.0 (gD3, anomS1), and one of 1.6
# as 2.9-3.8 (k5; review/PROPOSAL-laplace-check-every-unit-2026-10-27.md). So
# past one standard error the message gives a range, not a size. Quoted as
# "up to" its value, anomS1's warning said 1 standard error where the
# posterior mean was 1.9 posterior sds away and Laplace's interval held 1 per
# cent of the mass (review/REVIEW-fit-uncertainty-2026-10-29.md, section 6).
# The gap is not the trigger -- its level does not predict a move -- but a
# large one says Laplace is off throughout, which is why the fit reports the
# integrated log likelihood (`.ctLaplaceCheckedLoglik()`).
.ctLaplaceShiftMessage <- function(shift, worst, carrying, units,
  gap = NA_real_, gap_se = NA_real_) {
  sample <- "ctFitUncertainty(fit, 'sample', intoverpop = 'none')"
  move <- if (shift < 1) {
    paste0("by about ", format(shift, digits = 2), " standard errors here (", worst, ")")
  } else {
    paste0("by a standard error or more here (", worst, ": ", format(shift, digits = 2),
      " by the check's linear step, which is no guide to the size past one)")
  }
  resolved <- is.finite(gap) && abs(gap) >= 0.05 &&
    (!is.finite(gap_se) || abs(gap) > 2 * gap_se)
  paste0("Laplace's error moves the estimates ", move, ", carried by ", carrying,
    " of ", units, " units",
    if (resolved) paste0(", and puts Laplace's log likelihood ", format(abs(gap), digits = 2),
      " nats ", if (gap < 0) "above" else "below", " the integrated one, which the fit ",
      "reports"),
    ". The estimates and their intervals are Laplace's, ",
    "uncorrected. ctFitUncertainty(fit, intoverpop = 'quadrature') refines the ",
    "fit with those units integrated by quadrature, and ",
    sample,
    " samples the joint posterior.",
    " See fit$laplace$correction$check.")
}

# The log posterior the quadrature correction would report at `values`, from
# its own screen (`ctsem_laplace_continuation`): quadrature for every unit,
# plus the prior. `NA` when the quadrature cannot be evaluated. Used to compare two
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
    laplace = NA_real_, units = NA_real_))
  info <- .ctJuliaGet(module$ctsem_laplace_continuation_info(cont))
  c(corrected = as.numeric(info$quadrature)[1L],
    laplace = as.numeric(info$laplace)[1L], units = as.numeric(info$nunits)[1L])
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

# `correct = FALSE` is the default fit's (the random-effects check): the
# shift check alone, and where it finds Laplace's error moving an estimate, a
# warning naming the corrections and the estimate left where it is. `force`
# adds units to those the refined fit integrates whatever the screen says: a
# correction asked for after such a fit forces the units its check named.
.ctLaplaceContinue <- function(fit, cores = 1L, verbose = 0L,
  control = .ctLaplaceContinueDefaults, method = NULL, uncertaintycontrol = NULL,
  correct = TRUE, force = NULL) {
  started <- proc.time()[["elapsed"]]
  seconds <- function() proc.time()[["elapsed"]] - started
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  nsubjects <- length(fit$model_spec$subject_starts)
  module <- .ctJuliaModule(fit$model_spec$project)
  get <- .ctJuliaGet
  # The whole correction's progress line: the screen, the climb and its
  # Hessian all report through this one sink, at the same default verbosity as
  # the optimiser's own line, rather than only when `verbose > 0` -- which is
  # how a 22-minute correction (screen, a climb of 4.8 min, then the Hessian
  # 18.8 min) printed nothing at all by default. One sink for the
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
  # flagged whatever the screen finds, with the points the line search saw as
  # components of their proposals should they take draws; the climb then
  # starts from the jump point, which is the fit's estimate.
  jump <- fit$laplace$jump
  # Asked only to check, at a jump: no optimum to take the check's step from,
  # and the jump says what the check would, that a unit's term is not one
  # Gaussian basin's. Said, with the call that integrates the units named.
  if (!correct && !is.null(jump)) {
    named <- as.integer(jump$units)
    msg <- sprintf(paste0("The fit stopped where the Laplace term of %s %s ",
      "jumps, at no maximum, so its estimate and standard errors are ",
      "approximate. ctFitUncertainty(fit, intoverpop = 'quadrature') refines ",
      "the fit with %s integrated by quadrature; ctFitUncertainty(fit, 'sample', ",
      "intoverpop = 'none') ",
      "samples the joint posterior."), if (length(named) == 1L) "unit" else "units",
      paste(named, collapse = ", "), if (length(named) == 1L) "it" else "them")
    warning(msg, call. = FALSE)
    fit$laplace$correction <- list(method = "quadrature", status = "jump",
      applied = FALSE, units = named, message = msg, laplace_estimate = est,
      loglik_laplace = as.numeric(fit$estimate$loglik),
      seconds = c(total = seconds()))
    return(fit)
  }
  # Where every unit's integrand is Gaussian in its effects
  # (`laplace$quadratic`), each Laplace term is that unit's exact integral, at
  # any width: nothing to score and nothing to move. The screen's product
  # rules used to run here for nothing.
  if (is.null(jump) && isTRUE(fit$model_spec$laplace$quadratic)) {
    fit$laplace$correction <- list(method = "quadrature", status = "exact",
      applied = FALSE, nodes = as.integer(control$nodes),
      tolerance = as.numeric(control$tolerance), screen = 0, gap = 0,
      laplace_estimate = est, loglik_laplace = as.numeric(fit$estimate$loglik),
      flagged = 0L,
      error = list(nats = 0, se = 0, draws = 0, nodes = 0, residual = 0,
        units = NULL),
      seconds = c(total = seconds()))
    return(fit)
  }
  # The shift check first, on every unit whatever its width: where correcting
  # Laplace would move no estimate by `material` standard errors, resolved
  # past the check's own noise, the fit is left alone; otherwise the units
  # carrying the move are flagged as jump units are, and the refined fit
  # climbs. Not at a jump, whose units are named already and whose
  # point is no optimum to take a step from.
  checkfailed <- NULL
  check <- if (is.null(jump) && (!correct || !isFALSE(control$check))) {
    if (!is.null(sink)) sink(sprintf("Random-effects check, every unit on %d draws | %8s",
      as.integer(.ctJuliaOr(control$check_draws, 128L)), .ctDuration(seconds())), "update")
    tryCatch(.ctLaplaceShiftCheck(fit,
      draws = as.integer(.ctJuliaOr(control$check_draws, 128L)),
      seeds = seq_len(as.integer(.ctJuliaOr(control$check_seeds, 4L))),
      bar = as.numeric(control$material), control = control),
      error = function(e) {
        checkfailed <<- .ctLaplaceTryReason(structure("", condition = e))
        NULL
      })
  }
  checkrecord <- .ctLaplaceShiftRecord(check, fit)
  # Asked only to check, and the check could not be made: said, unless the
  # fit has no curvature to check against, which its uncertainty stage has
  # said already.
  if (!correct && is.null(check)) {
    if (!is.null(sink)) sink("", "break")
    if (!is.null(checkfailed) || !is.null(.ctBackendStoredHessian(fit, est))) {
      warning("Laplace's error could not be checked at this estimate",
        if (!is.null(checkfailed)) paste0(" (", checkfailed, ")"),
        "; ctLaplaceCheck(fit) measures it.", call. = FALSE)
    }
    fit$laplace$correction <- list(method = "quadrature", status = "unchecked",
      applied = FALSE, message = checkfailed, laplace_estimate = est,
      loglik_laplace = as.numeric(fit$estimate$loglik),
      seconds = c(check = seconds(), total = seconds()))
    return(fit)
  }
  # Left alone where the check resolved no move, and, asked only to check,
  # where it resolved none either way: a step that may or may not reach
  # `material` is no reason to warn.
  if (identical(check$status, "passed") ||
      (identical(check$status, "unresolved") && (!correct || !length(check$carrying)))) {
    # Left alone to the bit, as a passed screen leaves it: "exact" where the
    # corrected terms equal Laplace's to the screen's tolerance, "checked"
    # where they differ without moving the estimate. Every unit was looked
    # at, and the estimate's distance from the
    # exact-marginal optimum is the check's resolved step, not the bound a
    # value gap of that size would give.
    gap <- abs(check$gap)
    if (!is.finite(gap) || gap < 1e-8) gap <- 0
    record <- list(method = "quadrature",
      status = if (gap <= control$tolerance) "exact" else "checked",
      applied = FALSE, nodes = as.integer(control$nodes),
      tolerance = as.numeric(control$tolerance), screen = gap, gap = check$gap,
      laplace_estimate = est, loglik_laplace = as.numeric(fit$estimate$loglik),
      flagged = 0L, units = as.integer(check$units), check = checkrecord)
    bound <- if (gap > 0) max(abs(check$shift) + 2 * check$noise) else 0
    record$error <- list(nats = gap, se = bound, draws = check$gap_se, nodes = 0,
      residual = 0, units = NULL)
    record$seconds <- c(check = seconds(), total = seconds())
    fit$laplace$correction <- record
    return(.ctLaplaceCheckedLoglik(fit, check))
  }
  # A move found and asked only to check: said once, with the calls that
  # correct it, and kept in the record for print().
  if (!correct) {
    if (!is.null(sink)) sink("", "break")
    said <- .ctLaplaceShiftMessage(check$max_shift, checkrecord$worst,
      length(check$carrying), check$units, gap = check$gap, gap_se = check$gap_se)
    warning(said, call. = FALSE)
    fit$laplace$correction <- list(method = "quadrature", status = "flagged",
      applied = FALSE, message = said, laplace_estimate = est,
      loglik_laplace = as.numeric(fit$estimate$loglik),
      units = as.integer(check$units), check = checkrecord,
      seconds = c(check = seconds(), total = seconds()))
    return(.ctLaplaceCheckedLoglik(fit, check))
  }
  cont <- try(do.call(module$ctsem_laplace_continuation,
    c(list(.ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(est)),
      .ctLaplaceContinuationArgs(control,
        force = union(union(jump$units, check$carrying), as.integer(force)),
        points = jump$points))),
    silent = TRUE)
  # The engine's reason, which is what tells an allocation it refused from a
  # point the model cannot evaluate.
  if (inherits(cont, "try-error")) return(failed(paste0(
    "the quadrature could not be evaluated (", .ctLaplaceTryReason(cont), ")")))
  start <- get(module$ctsem_laplace_continuation_info(cont))
  # `loglik_laplace` is the Laplace objective's log likelihood at the point the
  # correction starts from, from the placement's own Laplace terms.
  record <- list(method = "quadrature", status = "exact", applied = FALSE,
    nodes = as.integer(start$nodes), tolerance = as.numeric(start$tolerance),
    screen = as.numeric(start$screen), gap = as.numeric(start$gap),
    laplace_estimate = est, loglik_laplace = sum(as.numeric(start$laplace_units)),
    flagged = as.integer(start$nflagged), units = as.integer(start$nunits),
    rule_failures = as.integer(start$rule_failures),
    lattice = as.integer(start$nlattice))
  record$check <- checkrecord
  # What the objective is known to at the screen; the paths that place again
  # replace it with the final placement's.
  record$error <- .ctLaplaceIntegrationError(start)
  # Wall seconds by stage as each is passed -- the screen, the climb, its
  # Hessian, the certification, the uncertainty redrawn from it -- and `total`
  # when the record is written (`.ctLaplacePlacedRefine()` adds its stages).
  stages <- c(screen = seconds())
  timed <- function() c(stages, total = seconds())
  since <- function(from) proc.time()[["elapsed"]] - from
  if (!is.finite(record$screen)) {
    return(failed("the quadrature was not finite at the estimate"))
  }
  # Passed: Laplace is exact here to
  # the tolerance, and the fit is left alone to the bit. Not at a jump: the
  # screen says the hybrid agrees with Laplace at the point, which is no
  # optimum of either (round 9's jump point passed the screen with a gradient
  # of 2.5e5), and the climb is what leaves it. Nor where the shift
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
  # The rules placed at every evaluation and one objective climbed from the
  # estimate (`.ctLaplacePlacedRefine()`); asked to stay (`stay = TRUE`: a fit
  # already refined, or one a caller wants left where it is), the rules placed
  # at the estimate, its Hessian and the uncertainty drawn there.
  .ctLaplacePlacedRefine(fit, module, cont, est, record, control,
    tol = .ctBackendGapTolerance(fit), stay = isTRUE(control$stay), chunks = chunks,
    verbose = verbose, sink = sink, method = method,
    uncertaintycontrol = uncertaintycontrol, report = report_quadrature,
    timed = timed)
}
