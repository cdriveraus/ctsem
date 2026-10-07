# Function to compute the Hessian using bootstrap resampling
bootstrapHessian <- function(standata, sm, est, finishsamples, cores, scores=NULL) {
  # `scores` may be supplied by a backend that computes per-subject gradients
  # itself (see .ctBackendScoreMatrix); only the Stan path has to reconstruct
  # them a subject at a time.
  if(is.null(scores)) scores <- scorecalc(standata = standata,est = est,stanmodel = sm,
    subjectsonly = ctOptimNSubjects(standata) >= 2,
    returnsubjectlist = F,cores=cores)
  num_bootstrap_samples <- max(c(finishsamples,1000))
  alpha_max = 100 # Maximum bootstrap sample size factor
  alpha_min = 1 # Minimum bootstrap sample size factor
  n_threshold=1000 # Threshold n for alpha correction
  alpha <- alpha_max - (alpha_max - alpha_min) * (min(1000,nrow(scores))  / n_threshold)  # Bootstrap sample size factor
  num_bootstrap_samples  # Total number of bootstrap samples
  n <- nrow(scores)  # Number of observations
  p <- ncol(scores)  # Number of parameters
  
  # Create a bootstrap resampling matrix
  resample_matrix <- matrix(sample(1:n, size = round(alpha * n) * num_bootstrap_samples, replace = TRUE),
    nrow = num_bootstrap_samples, ncol = round(alpha * n))
  
  # Generate random weights for smoothing
  weights <- matrix(runif(length(resample_matrix), min = 0.1, max = 2),
    nrow = num_bootstrap_samples)
  
  # Aggregate gradients using matrix multiplication
  gradsamples <- matrix(0, nrow = num_bootstrap_samples, ncol = p)  # Initialize gradsamples
  
  # Compute the weighted sum of gradients for each bootstrap sample
  for (i in 1:num_bootstrap_samples) {
    gradsamples[i, ] <- colSums(scores[resample_matrix[i, ], , drop = FALSE] * weights[i, ])
  }
  
  # Trim outliers
  trim_percent <- 0#.05  # Proportion of outliers to trim
  if (trim_percent > 0) {
    lower <- apply(gradsamples, 2, quantile, probs = trim_percent, na.rm = TRUE)
    upper <- apply(gradsamples, 2, quantile, probs = 1 - trim_percent, na.rm = TRUE)
    gradsamples <- pmax(pmin(gradsamples, upper), lower)        }
  
  #  Compute the Hessian with alpha correction
  hess <- -corpcor::cov.shrink(gradsamples,verbose=FALSE) / alpha
  return(list(hess=hess,scores=scores))
}

# What the importance sampler actually delivered, said out loud.
#
# `ctLaplaceCorrect()` and `ctParticleCorrect()` already warn when their
# importance sampling ends short of its effective-sample-size target.
# `uncertainty='is'` did not, and it is the call that most needs it, because
# the condition that defeats it is ordinary rather than exotic: a likelihood
# that is flat in some raw direction -- an individually varying parameter with
# no individual differences behind it, and the raw correlation that goes with
# it -- makes the target improper in that direction. The weights then have
# infinite variance, so the effective sample size does not converge to
# anything, and the run burns every one of `imisMaxIter` iterations to arrive
# at an effective sample of a few dozen. Measured on a 400-subject, 8-wave
# linear fit whose `rawcor` coordinate carries no curvature: 51,000
# log-probability evaluations, an ESS oscillating between 1.1 and 64.8 against
# a target of 100, and standard errors of 873 and 1214 on the two population
# standard deviations the Hessian puts at 1.5 and 0.8.
#
# The second half matters as much. When the weighted covariance comes back
# non-finite -- which it does when proposal draws land where the model cannot
# be evaluated -- both backends fall back to the unweighted covariance of the
# resampled draws. That is a different estimator, and until now nothing on the
# fit or in the session said which of the two had produced the intervals.
# "These intervals rest on too few effective points to mean what they look like
# they mean."
#
# Four places asked that and each wrote its own answer: here, in
# `ctLaplaceCorrect()`, and twice in `ctParticleCorrect()`. Three of the four
# did not call this one, and the thresholds had already drifted -- so a change
# to the rule, or to the wording, had to be made in up to four places to stay
# consistent, and nothing made it.
#
# `floor` is an argument because the rule differs by caller: a run that had a
# target warns when it ends short of it (`.ctOptimImisReport()`, which passes
# the target and names it), while reweighting one fixed batch of draws has no
# target and `max(50, 0.1 * n)` is a rule about that batch. `remedy` is the
# other real per-caller part -- what to do instead differs by where you are,
# and naming the wrong alternative is worse than naming none.
#' @keywords internal
.ctOptimEffectiveSampleWarn <- function(ess, floor, remedy, ndraws = NULL,
  target = NULL){
  ess <- if(is.null(ess)) NA_real_ else as.numeric(ess)[1L]
  if(!is.finite(ess) || !is.finite(floor) || ess >= floor) return(invisible(ess))
  warning('Importance sampling reached an effective sample size of ',
    round(ess, 1),
    if(is.null(ndraws)) '' else paste0(' from ', ndraws, ' draws'),
    if(is.null(target)) '' else paste0(', short of its target of ', target),
    '. The intervals rest on that many points, not on the number of draws. ',
    remedy, call.=FALSE)
  invisible(ess)
}

# `tailremedy` is what to do when the Pareto k below says the weights cannot
# be trusted, which is a different failure from a short effective sample and
# can have a different answer: reweighting cannot recover what the proposal
# never visits, so the remedy is to sample the target itself where that exists.
#
# Returns the effective sample size and Pareto k, for the caller to record.
.ctOptimImisReport <- function(is_res, target, weighted,
  remedy = paste0('A direction the data does not identify cannot be importance ',
    'sampled at all -- check the identifiability report, and consider ',
    'uncertainty = "hessian".'), tailremedy = remedy){
  ess <- if(is.null(is_res$ess)) NA_real_ else as.numeric(is_res$ess)[1L]
  k <- .ctImisParetoK(is_res)
  if(!isTRUE(weighted)) warning(
    'The weighted importance-sampling covariance was not finite, so the ',
    'unweighted covariance of the resampled draws was used instead.',
    call.=FALSE)
  if(is.finite(k) && k > .ctImisParetoKBar()) {
    warning('Importance sampling weights have Pareto k ', round(k, 2),
      ', above ', .ctImisParetoKBar(), ': the proposal misses part of the ',
      'posterior, so these intervals are unreliable whatever the effective ',
      'sample size. ', tailremedy, call.=FALSE)
  } else if(is.na(k) && length(is_res$log_weights) &&
      !requireNamespace('loo', quietly = TRUE)) {
    message('Pareto k of the importance weights not checked: install the loo ',
      'package.')
  }
  # A profile path that had not fallen off at the last rung: see
  # `.ctImisPaths()`. Positional, since names are attached later; in the
  # identified subspace the paths run along its directions, not parameters.
  unbounded <- is_res$paths$reachesLimit
  if(length(unbounded)) {
    what <- if(is.null(attr(is_res, 'subspace'))) 'raw parameter' else
      'direction of the identified subspace'
    warning('The posterior had not fallen off ', is_res$paths$limit,
      ' standard errors of the ',
      'curvature out along ', what, if(length(unbounded) > 1) 's ' else ' ',
      paste0(abs(unbounded), ifelse(unbounded < 0, ' (below)', ' (above)'),
        collapse = ', '),
      ', so it may be improper there, which no reweighting can represent. ',
      tailremedy, call.=FALSE)
  }
  # Short of the target at all, not of half of it: the rounds stop on reaching
  # it, so ending below it means they ran out, and a run that aimed for 200
  # and stopped at 150 should say so rather than pass silently.
  .ctOptimEffectiveSampleWarn(ess,
    floor = if(is.finite(target)) target else NA_real_, remedy = remedy,
    target = if(is.finite(target)) target else NULL)
  invisible(list(ess = ess, k = k))
}

# Whether importance weights can be trusted at all, as the shape k of a
# generalized Pareto fitted to their upper tail (Vehtari, Simpson, Gelman, Yao
# and Gabry, "Pareto smoothed importance sampling"), through `.ctPsisWeights()`,
# the fit `ctLOO(method = 'psis')` already makes, rather than a second one.
#
# The effective sample size cannot answer this. It is computed from the same
# weights, so when their variance does not exist -- a proposal with lighter
# tails than its target, which a normal proposal built on the curvature at the
# mode is whenever the posterior is skewed or heavy-tailed -- it still returns
# a healthy-looking number while the draws that were never made carry the
# missing mass. k estimates the tail directly: below 0.5 the variance exists,
# above 0.7 the estimates are not usable at any practical draw count, which is
# loo's rule and the one `ctLOO()` already reports against.
#
# NA when loo is not installed (it is suggested, not imported), or when there
# are too few weights to fit a tail to.
.ctImisParetoKBar <- function() 0.7

.ctImisParetoK <- function(is_res) {
  lw <- as.numeric(is_res$log_weights)
  lw <- lw[is.finite(lw)]
  if(length(lw) < 20L || !requireNamespace('loo', quietly = TRUE)) return(NA_real_)
  k <- tryCatch(.ctPsisWeights(matrix(lw, ncol = 1L))$k,
    error = function(e) NA_real_)
  as.numeric(k)[1L]
}

# Importance sampling against a reference density, and the covariance and draws
# that come out of it.
#
# Three places do this, and they are the same six lines each time -- `imis_is`,
# then the weighted covariance with an unweighted fallback, then the
# effective-size check. What differs is only which density is handed in and what
# to suggest when the effective size is short:
#
#   the uncertainty stage's `uncertainty='is'`, against the model's own density;
#   `ctLaplaceCorrect(draws='imis')`, against the adaptive-quadrature posterior;
#   `ctParticleCorrect(draws='imis')`, against the particle-filter likelihood.
#
# The last two are corrections *to a different objective* -- the reference is
# more accurate than what was optimised -- where the first reweights the same
# one. That is a real difference in what the answer means, and none in how it is
# computed, which is why only this part is shared.
#
# The IMIS proposal-inflation scale, named once. Before this it was written
# four times: `.ctOptimDrawSamples()`'s own formal defaults (the stan path
# reaches them by not overriding), `.ctOptimImisDraws()`'s own formal defaults
# (also stan-shaped, and also unused -- both its callers pre-scale their own
# proposal and pass `scaleInit = 1`), `.ctBackendUncertainty()`'s julia call
# site, and `imis_is()`'s own formal defaults -- the last two matching each
# other by coincidence rather than by reading from one place, which is exactly
# the shape a later edit to one number and not the other drifts through
# unnoticed.
#
# The two backend values are both deliberate and both measured, on different
# regimes: 1.1 on a 400-subject fit's identified subspace beats 1.5 at every
# evaluation count (2,000 evaluations to reach ESS 149.6 against 4,000 to reach
# 112.9, `IS-importance-sampling-2026-09-06.md` s5), which is why stan, whose
# models in the pass-3 benchmark are of that shape, uses it. 1.5/1.2 was raised
# for a 40-subject model with a nonlinear or variance parameter, where the
# posterior is genuinely wider than the Hessian curvature and a narrower
# proposal cannot see the extra width at all (`.ctBackendUncertainty()`'s own
# comment has the fuller account); whitening the null directions out does not
# change that argument; it only stops a different, unrelated failure
# (`.ctImisSubspace()`) from also being blamed on the scale. One number still
# cannot serve both sample sizes, so this stays two named constants rather than
# collapsing to one -- the point of naming them once is that this file is the
# only place a change to either has to be made.
#
# A bare call to `imis_is()` (a dev script, or the reproduction in the IS note
# itself) has no backend to ask, so its own formal defaults use the julia
# value: the wider, more conservative proposal, on the reasoning that costs
# more evaluations rather than the one that can quietly under-cover.
#
# `df`, the components' t degrees of freedom, is one value for both: the
# heavier tail is what lets a proposal built on the curvature reach draws its
# normal would not, and it does not depend on which engine evaluates them.
.ctImisProposalDefaults <- function(backend = c('julia', 'stan')) {
  backend <- match.arg(backend)
  if (backend == 'stan') list(scaleInit = 1.1, tailScale = 1.1, df = 5) else
    list(scaleInit = 1.5, tailScale = 1.2, df = 5)
}

# The effective sample size every route that produces draws aims for unless
# told otherwise, named once: `uncertainty = 'sample'` (`minESS`, of the worst
# parameter), `uncertainty = 'is'` (`isESS`), and the importance-sampling
# corrections `ctLaplaceCorrect()` and `ctParticleCorrect()` (`target_ess`,
# whose exported signatures write the number out; test-ess-target.R holds them
# to this). They differed -- 100 on the corrections, 200 elsewhere, and
# 'sample' stopping at its draw budget below its own 200 -- and one default
# across them was Charles's decision (2026-09-30).
#
# 200 because the 2.5% and 97.5% quantiles `summary()` reports are what needs
# the draws, for every parameter; and at 100 importance sampling stopped with a
# variance's tail still resting on a handful of heavy draws, its posterior sd
# moving by a third between seeds (job M2, `imis_is()`). A route that ends
# short of it says so (`.ctOptimImisReport()`, `.ctSampleWarn()`).
.ctEssTarget <- 200

# Which directions of a proposal covariance carry no information at all, by
# the same test `.ctOptimIdentifiedInverse()` applies to an information matrix
# and the same default tolerance, `.ctFlatDirectionRtol()` -- so the two agree
# on what "rank deficient" means, even though one looks at a covariance's small
# eigenvalues and the other at an information matrix's. That is the right
# correspondence rather than a coincidence: a direction `.ctOptimIdentifiedInverse()`
# projects out of the information contributes exactly zero variance to the
# covariance it builds (`vectors %*% (t(vectors) / values[keep])`, summed over
# kept directions only), so the same direction shows up here as a covariance
# eigenvalue at or near zero, not as a large one -- there is no ridge-floored
# covariance reaching this code any more (see `IS-importance-sampling-2026-09-06.md`
# s2 for what one looked like before that was fixed).
#
# Returns NULL when nothing is dropped -- the ordinary fit, and the common
# case -- so a caller pays for one extra `eigen()` and nothing else. That is
# negligible next to what it is about to spend on `imis_is`.
.ctImisSubspace <- function(cov, rtol = .ctFlatDirectionRtol()) {
  cov <- as.matrix(cov)
  cov <- (cov + t(cov)) / 2
  eig <- try(eigen(cov, symmetric = TRUE), silent = TRUE)
  if (inherits(eig, 'try-error')) return(NULL)
  values <- eig$values
  scale <- max(values)
  if (!is.finite(scale) || scale <= 0) return(NULL)
  keep <- values > rtol * scale
  if (!any(keep) || all(keep)) return(NULL)
  V <- eig$vectors[, keep, drop = FALSE]
  d <- values[keep]
  nullvectors <- eig$vectors[, !keep, drop = FALSE]
  # Same convention as `.ctOptimIdentifiedInverse()`'s `nullParameters`: the
  # coordinates with a share of the dropped directions at `.ctNullMassBar()`.
  loaded <- if (ncol(nullvectors))
    which(rowSums(nullvectors^2) >= .ctNullMassBar()) else integer()
  list(V = V, d = d, k = sum(keep), n = nrow(cov), nnull = sum(!keep),
    nullEigenvalues = values[!keep], nullParameters = loaded)
}

# The affine map from a `subspace$k`-dimensional whitened coordinate `z` back
# onto the raw parameters `subspace` was built from: `theta = centre + A %*% z`
# with `A = V %*% diag(sqrt(d))`, so `z` has covariance `I` exactly when `theta`
# has the covariance `.ctImisSubspace()` was handed (restricted to the kept
# directions; a null direction gets no column in `A` at all, so it can never
# move away from `centre`). Takes and returns a matrix with draws as ROWS,
# `imis_is`'s own convention for `x_new`.
.ctImisUnwhitenMatrix <- function(Z, centre, subspace) {
  A <- sweep(subspace$V, 2, sqrt(subspace$d), '*')
  sweep(as.matrix(Z) %*% t(A), 2, as.numeric(centre), '+')
}

# `lpg`, wrapped to take a whitened `z` and evaluate the original density at
# the raw point it maps to. Carries `lpg`'s own `'batch'` attribute through the
# same map -- see `.ctBackendLpgFunc()` -- so whitening does not silently
# defeat the one-bridge-call route: without this, `imis_is` would find no
# `'batch'` attribute on the wrapped closure and fall back to its per-draw
# loop, which is the cost this whole repair exists to remove.
#
# The `'gradbatch'` attribute goes through the map too, as `A' g`, for the
# profile-path search in `imis_is()`; without it the search would silently
# not run in the identified subspace.
.ctImisWhitenDensity <- function(lpg, centre, subspace) {
  A <- sweep(subspace$V, 2, sqrt(subspace$d), '*')
  wrapped <- function(z) lpg(as.numeric(
    .ctImisUnwhitenMatrix(matrix(z, nrow = 1L), centre, subspace)))
  batchbase <- attr(lpg, 'batch')
  if (!is.null(batchbase)) {
    attr(wrapped, 'batch') <- function(Z) batchbase(.ctImisUnwhitenMatrix(Z, centre, subspace))
  }
  gradbase <- attr(lpg, 'gradbatch')
  if (!is.null(gradbase)) {
    attr(wrapped, 'gradbatch') <- function(Z) {
      r <- gradbase(.ctImisUnwhitenMatrix(Z, centre, subspace))
      r$gradient <- as.matrix(r$gradient) %*% A
      r
    }
  }
  wrapped
}

# `imis_is()`'s own result, computed in whitened coordinates, mapped back onto
# the raw parameters. `mean` and `covariance` are recomputed from the mapped
# draws and the run's own final weights with the same `diagis` functions
# `imis_is` used internally -- not transformed analytically -- so this does not
# depend on separately re-deriving how a covariance moves under an affine map;
# it just asks the same weighted moment for the same weights at the mapped
# points. A weighted mean and variance are both affine-equivariant, so the two
# routes agree exactly; recomputing is the one that cannot get the direction of
# a transpose wrong.
.ctImisUnwhitenResult <- function(result, centre, subspace) {
  n <- length(centre)
  raw <- function(Z) if (length(Z)) .ctImisUnwhitenMatrix(Z, centre, subspace) else
    matrix(numeric(0), 0, n)
  result$theta <- raw(result$theta)
  result$full_theta <- raw(result$full_theta)
  haveweights <- length(result$full_weights) > 0
  result$mean <- if (haveweights)
    as.numeric(diagis::weighted_mean(result$full_theta, result$full_weights)) else
    rep(NA_real_, n)
  result$covariance <- if (haveweights)
    diagis::weighted_var(result$full_theta, result$full_weights) else
    matrix(NA_real_, n, n)
  result
}

# `imis_is()`, run in the identified subspace of `cov` when it is rank
# deficient (see `.ctImisSubspace()`) and passed straight through otherwise.
# Returns exactly what `imis_is()` returns -- a drop-in replacement at both of
# its call sites -- with the subspace, if any was used, carried as the
# `'subspace'` attribute rather than a new list field: nothing reads
# `attr(is_res, ...)` today, so this cannot collide with an existing `$name`
# read the way a new list field risks doing under R's partial matching (the
# `$ctstanmodel`/`$ctstanmodelbase` incident is the reason to say this
# explicitly rather than assume it).
#
# Sampling only the identified subspace is what removes the infinite-variance
# weights a flat raw direction produces: measured on the 400-subject fit
# `IS-importance-sampling-2026-09-06.md` is written against, 51,000 evaluations
# that never reached an effective sample of 100 became 2,000 to 4,000 that did,
# and the flat coordinate's reported spread went from 397 to exactly 0 rather
# than from 0 to 397 by accident -- see `.ctImisUnwhitenMatrix()`: a null
# direction has no column in the map at all.
.ctImisRun <- function(lpg, centre, cov, rtol = .ctFlatDirectionRtol(), ...) {
  centre <- as.numeric(centre)
  cov <- as.matrix(cov)
  subspace <- .ctImisSubspace(cov, rtol = rtol)
  if (is.null(subspace)) return(imis_is(lpg, mu_hat = centre, Sigma_hat = cov, ...))
  wrapped <- .ctImisWhitenDensity(lpg, centre, subspace)
  whitened <- imis_is(wrapped, mu_hat = rep(0, subspace$k),
    Sigma_hat = diag(subspace$k), ...)
  result <- .ctImisUnwhitenResult(whitened, centre, subspace)
  attr(result, 'subspace') <- subspace
  result
}

# The gradient a density can give `imis_is()`'s profile-path search: its
# `'gradbatch'` attribute, a function of a draws matrix (rows as draws)
# returning list(value, gradient) with one gradient row per draw, or NULL.
# Read from the attribute and never probed for: a density says what it can do
# by its attributes, as `'batch'` does, so the julia route's lpg carries one
# (one bridge call per batch), stan's gets one from `.ctImisPointGradbatch()`,
# and the quadrature and particle-filter densities of `ctLaplaceCorrect()` and
# `ctParticleCorrect()` carry none and skip the search. Probing would cost a
# particle filter per call there, and put a row in its evaluation record.
.ctImisGradient <- function(parlp, centre) {
  gb <- attr(parlp, 'gradbatch')
  if (is.function(gb)) gb else NULL
}

# A per-point density whose value carries a `'gradient'` attribute -- stan's
# `ctOptimFitLpgFunc()` -- given the `'gradbatch'` attribute `.ctImisGradient()`
# reads, asking one point at a time. So both backends run the same search.
.ctImisPointGradbatch <- function(lpg) {
  attr(lpg, 'gradbatch') <- function(X) {
    X <- as.matrix(X)
    out <- lapply(seq_len(nrow(X)), function(i)
      tryCatch(lpg(X[i, ]), error = function(e) NA_real_))
    grads <- lapply(out, function(o) {
      g <- attr(o, 'gradient')
      if (is.null(g) || length(g) != ncol(X)) rep(NA_real_, ncol(X)) else as.numeric(g)
    })
    list(value = vapply(out, function(o) as.numeric(o)[1L], numeric(1)),
      gradient = matrix(unlist(grads), nrow = nrow(X), byrow = TRUE))
  }
  lpg
}

# Each parameter's profile path, walked out on both sides of the mode: at
# `first` and then `rungs` standard errors of the curvature, parameter j is
# held and the others moved to their conditional mode, by Newton steps in the
# complement with the mode's precision block as the metric, from the previous
# rung's point carried on linearly (the regression direction at the first).
# See `imis_is()` for why, and for what is built from the points.
#
# A path stops when its log density has fallen `maxdrop` below the mode's,
# when it cannot be evaluated, or after the first rung if it fell there by
# more than `heavy` of what the curvature predicts (`first^2 / 2`): such a
# path has no tail to follow, and on the models measured that is most of
# them, which is what keeps the search cheap. A path still above `maxdrop` at
# the last rung is reported in `reachesLimit`, as parameter index times side
# (negative for the lower), with the last rung as `limit`: the density there
# has not fallen off that many standard errors out, which is what an
# improper direction looks like. The ladder runs to 96 because a variance
# held up by nothing but its N(0, 1) raw prior was measured 24 standard
# errors out still under the drop limit, proper but long.
#
# Every active path moves one Newton step per iteration, so an iteration is
# one call of `gradfun` for all of them and one or more of `evaluate` for the
# step halving.
.ctImisPaths <- function(gradfun, evaluate, mu, S, first = 3,
  rungs = c(6, 12, 24, 48, 96), maxdrop = 8, steps = 12, tol = 1e-3, heavy = 0.6) {
  d <- length(mu)
  ladder <- c(first, rungs)
  P <- solve(S)
  se <- sqrt(diag(S))
  # The paths start from the density's own mode, which the point handed in
  # need not be: a Laplace fit whose estimate the quadrature correction moved
  # reports that corrected point, up to a standard error from the mode of the
  # Laplace objective this samples (measured: 0.98 on one parameter of a
  # 40-subject model). Newton steps with the curvature handed in, kept only
  # while they gain.
  lp0 <- evaluate(matrix(mu, 1L))
  ngrad <- 0L
  for (st in 1:8) {
    g <- gradfun(matrix(mu, 1L))
    ngrad <- ngrad + 1L
    gr <- as.numeric(g$gradient)
    if (!is.finite(g$value[1L]) || any(!is.finite(gr))) break
    step <- as.numeric(S %*% gr)
    # the step and three halvings of it, in one batch
    cand <- t(vapply(2^-(0:3), function(a) mu + a * step, numeric(d)))
    if (d == 1L) cand <- matrix(cand, ncol = 1L)
    v <- evaluate(cand)
    best <- which.max(v)
    if (!is.finite(v[best]) || v[best] <= lp0 + tol) break
    mu <- cand[best, ]
    lp0 <- v[best]
  }
  J <- rep(seq_len(d), each = 2L)
  side <- rep(c(1, -1), d)
  np <- length(J)
  chols <- if (d > 1) lapply(seq_len(d), function(j) chol(P[-j, -j, drop = FALSE])) else NULL
  rung <- rep(1L, np)
  iter <- rep(0L, np)
  active <- rep(TRUE, np)
  history <- replicate(np, list(mu), simplify = FALSE)
  # Where a path's next rung starts: parameter j moved to the rung, and the
  # others carried along the regression direction (the curvature's own
  # guess), extrapolated through the last two points of the path, or left
  # where the last point had them. All are evaluated and the best kept. The
  # curvature's guess alone is not enough: on a 40-subject Laplace model the
  # Hessian at the reported estimate put the start of a variance's first rung
  # 101 log units down, from which no Newton step recovered, where holding
  # the others still was 20 down and the path then found its ridge.
  starts <- function(i) {
    j <- J[i]
    v <- mu[j] + side[i] * ladder[rung[i]] * se[j]
    h <- history[[i]]
    n <- length(h)
    last <- h[[n]]
    cand <- list(last + (v - last[j]) * S[, j] / S[j, j], last)
    if (n >= 2L) {
      b <- h[[n - 1L]]
      cand[[3L]] <- last + (last - b) * (v - last[j]) / (last[j] - b[j])
    }
    out <- do.call(rbind, cand)
    out[, j] <- v
    out
  }
  place <- function(ids) {
    cl <- lapply(ids, starts)
    allc <- do.call(rbind, cl)
    v <- evaluate(allc)
    ncalls <<- ncalls + 1L
    at <- 0L
    for (m in seq_along(ids)) {
      nc <- nrow(cl[[m]])
      X[ids[m], ] <<- cl[[m]][which.max(v[at + seq_len(nc)]), ]
      at <- at + nc
    }
  }
  X <- matrix(rep(mu, each = np), np, d)
  points <- list()
  reaches <- integer(0)
  followed <- integer(0)
  ncalls <- 2L * ngrad
  place(seq_len(np))
  while (any(active)) {
    ids <- which(active)
    g <- gradfun(X[ids, , drop = FALSE])
    ngrad <- ngrad + length(ids)
    ncalls <- ncalls + 1L
    gv <- as.numeric(g$value)
    G <- as.matrix(g$gradient)
    ok <- is.finite(gv) & gv > -1e99 & apply(is.finite(G), 1L, all)
    dy <- matrix(0, length(ids), d)
    if (d > 1L) for (m in which(ok)) {
      j <- J[ids[m]]
      R <- chols[[j]]
      dy[m, -j] <- backsolve(R, forwardsolve(t(R), G[m, -j]))
    }
    newv <- gv
    accepted <- rep(FALSE, length(ids))
    pending <- ok & d > 1L
    a <- rep(1, length(ids))
    for (h in 0:6) {
      if (!any(pending)) break
      pm <- which(pending)
      cand <- X[ids[pm], , drop = FALSE] + a[pm] * dy[pm, , drop = FALSE]
      v <- evaluate(cand)
      ncalls <- ncalls + 1L
      better <- is.finite(v) & v > -1e99 & v >= gv[pm]
      acc <- pm[better]
      if (length(acc)) {
        X[ids[acc], ] <- cand[better, , drop = FALSE]
        newv[acc] <- v[better]
        accepted[acc] <- TRUE
        pending[acc] <- FALSE
      }
      a[pending] <- a[pending] / 2
    }
    advance <- integer(0)
    for (m in seq_along(ids)) {
      i <- ids[m]
      iter[i] <- iter[i] + 1L
      gain <- if (accepted[m]) newv[m] - gv[m] else 0
      if (ok[m] && accepted[m] && gain >= tol && iter[i] < steps) next
      value <- if (accepted[m]) newv[m] else gv[m]
      drop <- lp0 - value
      usable <- is.finite(value) && value > -1e99
      j <- J[i]
      if (usable && drop <= maxdrop) points[[length(points) + 1L]] <- list(j = j,
        side = side[i], rung = ladder[rung[i]], x = X[i, ],
        from = history[[i]][[length(history[[i]])]], drop = drop)
      last <- rung[i] == length(ladder)
      if (usable && drop <= maxdrop && last) reaches <- c(reaches, as.integer(j * side[i]))
      if (!usable || drop > maxdrop || last ||
          (rung[i] == 1L && drop > heavy * first^2 / 2)) {
        active[i] <- FALSE
        next
      }
      followed <- c(followed, j)
      history[[i]][[length(history[[i]]) + 1L]] <- X[i, ]
      rung[i] <- rung[i] + 1L
      iter[i] <- 0L
      advance <- c(advance, i)
    }
    if (length(advance)) place(advance)
  }
  list(points = points, followed = sort(unique(followed)), gradients = ngrad,
    calls = ncalls, reachesLimit = reaches, limit = max(ladder), centre = mu)
}

# `cov` is the proposal covariance as the caller wants it used. A caller that
# has already widened it passes `scaleInit = 1` rather than compounding two
# scalings, which is what `ctLaplaceCorrect()` does.
#' @keywords internal
.ctOptimImisDraws <- function(lpg, centre, cov, finishsamples, remedy,
  nbatch = 1000, target_ess = .ctEssTarget, maxiter = 50,
  scaleInit = .ctImisProposalDefaults('stan')$scaleInit,
  tailScale = .ctImisProposalDefaults('stan')$tailScale,
  df = Inf, verbose = 0, diagPlots = TRUE){

  is_res <- .ctImisRun(lpg, centre = centre, cov = cov,
    cl = NA, n_batch = as.integer(nbatch), target_ess = target_ess,
    max_iter = as.integer(maxiter), scale_init = scaleInit,
    tail_scale = tailScale, df = df,
    finishsamples = as.integer(finishsamples), diag_plots = diagPlots,
    # `verbose > 0`, not TRUE: this printed IMIS iteration progress at
    # `verbose = 0`, so the one argument meant two things across the backends --
    # silence on julia, a page of output on stan.
    verbose = verbose > 0)

  samples <- is_res$theta
  weighted <- !is.null(is_res$covariance) && all(is.finite(is_res$covariance))
  cov_out <- if(weighted) ctOptimSafeCov(is_res$covariance) else
    if(!is.null(samples) && nrow(samples) > 1) ctOptimSafeCov(stats::cov(samples)) else cov
  report <- .ctOptimImisReport(is_res, target_ess, weighted, remedy = remedy)

  list(samples = samples, cov = cov_out, ess = report$ess, k = report$k,
    weighted = weighted, is_res = is_res, subspace = attr(is_res, 'subspace'))
}

ctOptimSafeCov <- function(cov, ridge=1e-8){
  cov <- as.matrix(cov)
  cov <- (cov + t(cov)) / 2
  if(any(!is.finite(cov))) stop('Non-finite covariance values')
  diagnostics <- list(nearPD=FALSE, ridgeApplied=FALSE,
    minEigenOriginal=NA_real_, minEigenFinal=NA_real_, ridge=ridge)
  eig <- try(eigen(cov, symmetric=TRUE), silent=TRUE)
  if('try-error' %in% class(eig)) {
    cov <- as.matrix(Matrix::nearPD(cov, conv.norm.type='F')$mat)
    eig <- eigen(cov, symmetric=TRUE)
    diagnostics$nearPD <- TRUE
  }
  mineig <- min(eig$values)
  diagnostics$minEigenOriginal <- mineig
  if(mineig <= ridge){
    eig$values <- pmax(eig$values, ridge)
    cov <- eig$vectors %*% diag(eig$values, length(eig$values)) %*%
      t(eig$vectors)
    cov <- (cov + t(cov)) / 2
    diagnostics$ridgeApplied <- TRUE
  }
  diagnostics$minEigenFinal <- min(eigen(cov, symmetric=TRUE,
    only.values=TRUE)$values)
  attr(cov, 'ctOptimSafeCov') <- diagnostics
  cov
}

# Invert an information matrix over the subspace the data actually determines,
# leaving the rest alone.
#
# The alternative -- flooring every eigenvalue at a small `ridge` and inverting
# the result -- manufactures a variance of `1/ridge` along each direction that
# carries no information, and that number is a property of the ridge rather
# than of the data. With the default ridge of 1e-8 it is 1e8, a standard error
# of 1e4. Two things then go wrong, and both were measured on a benchmark fit
# (60 subjects x 150 occasions, a population SD with no individual differences
# behind it):
#
#   * The orientation of a null eigenvector is set by rounding error, because
#     the block it spans is numerically zero. It therefore has an arbitrary
#     small component on the *identified* parameters, and 1e8 multiplies that
#     component. Perturbing the information matrix by 1e-12 of its own scale --
#     less than the difference between two runs that reached the same optimum
#     to eight decimal places -- moved the reported standard error of a
#     well-determined drift parameter from 0.019 to 0.24, and a second such
#     perturbation to 0.087. A quantity that moves by an order of magnitude
#     under rounding is not a statement about the data.
#   * Nothing downstream can tell 1e4 from a real standard error, so it flows
#     into the draws, the quantiles and the mean over draws that `popmeans`
#     reports, which is how a fit at the right optimum came to report
#     intervals a hundred times too wide and point estimates that had wandered.
#
# So a direction whose curvature is negligible against the sharpest one is
# projected out instead of floored: the identified subspace is inverted
# exactly, and the null subspace contributes zero rather than 1/ridge. That
# leaves the identified parameters stable under rounding, which is the whole
# point, and it leaves the flat coordinates reported with no spread at all --
# which is why the null directions are recorded and named here, and warned
# about by the caller.
#
# `rtol` is deliberately *not* `.ctBackendIdentifiability()`'s 1e-8, and the
# difference is the point. That one is a judgement about identification and it
# only warns, so a false positive costs a warning. This one acts: a direction
# it drops comes back with no spread at all, which would be a new wrong answer
# for a direction that is weak but real. So it is set by numerical resolution
# instead. A symmetric eigendecomposition resolves eigenvalues to about
# `eps * largest`, so below ~1e-14 of the largest an eigenvalue is inside the
# error bar of zero; 1e-12 leaves two orders of margin above that and ten below
# the identification judgement. Every flat direction measured here sat between
# 1e-18 and 1e-25 of its matrix's largest eigenvalue, so nothing real is near
# this line. A direction between the two tolerances is inverted as usual --
# its variance is genuinely enormous, which is the truth about it -- and
# `.ctBackendIntervalCheck()` is what says so.
# Which flagged directions the likelihood is actually flat along.
#
# The curvature at the estimate is a local quadratic approximation, and on a
# direction the data does not determine it is measuring the wrong thing: the
# eigenvalue there is not a property of the model and the data, it is a
# residue of the transform's own derivative at whatever raw value the
# optimiser happened to stop at. Measured on a two-latent model fitted to
# noise, walking one diffusion correlation out along its flat ray: the log
# likelihood is -207.01897 at every raw value from -6 to -20, while the
# smallest relative eigenvalue falls from 1.6e-08 to 7.1e-16 and then turns
# negative from rounding. Whether that direction was reported as determined
# was therefore decided by where the optimiser stopped on a ray along which
# the likelihood is constant -- two runs of the same fit, one with the
# predicted-gain stopping rule on and one off, gave opposite diagnoses.
#
# So the eigenvalue selects candidates and the likelihood decides. This is the
# profile-likelihood criterion (Raue et al. 2009): a direction is not
# identified when the likelihood does not change along it, and the scale for
# "does not change" is the likelihood-ratio bound, `qchisq(1-alpha, 1) / 2` --
# 1.92 at 95% -- rather than a tolerance on a differentiated approximation.
# That bar is a statistical quantity, invariant to reparameterisation, and
# comparable across models, which no `rtol` on an eigenvalue is. The
# sloppy-model literature (Gutenkunst et al. 2007) is the general reason to
# expect no threshold to work: eigenvalue spectra are typically spread over
# many orders with no gap to cut at.
#
# Deliberately one-sided, and that is the whole of what makes it safe to act
# on. Walking a direction without re-optimising the other parameters is a
# slice, not a profile: if the likelihood stays flat we have *exhibited* a
# curve along which it is constant, which is non-identification and needs no
# further argument; if it rises we have learned nothing, because the flat
# manifold may be curved and a profile would have followed it. So a candidate
# that fails this test is left exactly as it was.
#
# ## Following a curved ridge, and why one side settles it
#
# The straight slice is where the answer used to depend on where the optimiser
# stopped. Measured on `test-stan-julia-parity.R`'s fixture -- six subjects, ten
# population correlations along one ridge -- at five stopping points on the
# same ridge, 4.4e-04 nats apart from first to last: the straight walk's worst
# drop at four raw units was 17.8, 19.9, 5.6, 2.1 and 0.68 nats. The ridge is
# flat throughout -- the optimiser walked six raw units along it for those
# 4.4e-04 -- but it is curved in raw coordinates, so a straight line leaves it,
# and how fast depends on where on it you start.
#
# So a rung that has dropped past the bar is followed back to the ridge before
# it is judged: up to `corrections` Newton steps in the directions the
# curvature does trust, with the Hessian already decomposed and the gradient at
# the rung, each kept only if it improves. The corrected point is still at the
# same displacement along the candidate direction -- the step is orthogonal to
# it -- and it is a point the model evaluated, so the profile at that
# displacement is at least as high. A corrected rung within the bar is
# therefore evidence of the same one-sided kind, and a rung the correction does
# not bring back is refused as before. On the fixture's earliest stopping
# point the four-unit rung's 2.51-nat drop came back to 1.39 in one step.
#
# And a side at a time: one side flat through the whole ladder exhibits a
# curve of constant likelihood running four raw units from the estimate, which
# is conclusive whatever the other side does -- the same reading
# `ctFitProfile()` gives a profile, one flat side being structural
# non-identification. Requiring both sides refused the fixture's earliest
# stopping point, whose ridge runs one way from where it stopped.
#
# ## Candidates
#
#   * curvature below `rtol = 1e-7` of the sharpest direction, ten times the
#     identifiability report's own eigenvalue rule. Because curvature along a
#     flat ray decays with where the optimiser stopped, that rule saw the
#     fixture's ridge at 1.7e-08, 4.1e-09, 1.2e-09, 6.8e-11 and 4.5e-12 over the
#     five stopping points -- above its 1e-8 at the first, so nothing was named
#     there -- while the fixture's first identified direction held at 1.2e-06
#     at every one. 1e-7 sits in that gap. Healthy fits cost nothing: on four
#     fits built as the test fixtures are (test-julia-convergence.R's,
#     test-julia-laplace.R's exact one, a 40-subject nonlinear Laplace fit,
#     and test-backend-summary.R's default one) the smallest relative
#     curvature was 9.6e-05 or more, apart from the last one's own flat
#     direction.
#     What this changes is that the report can now name a direction the
#     eigenvalue rule missed -- it names what the likelihood confirms, see
#     `.ctBackendIdentifiability()` -- and that such a direction leaves the
#     covariance as the rest do.
#   * with no candidate it returns `NULL` before evaluating anything, so a fit
#     with no flat direction -- the usual fit -- pays nothing and is unchanged.
#
# `lengths` are in raw parameter units, where ctsem's coordinates are
# standardised by construction, and are the same ladder
# engine's flat probe (`_ctsem_flat_probe`) walks for the sibling question ("does anything
# *improve* along here"). `maxdirections` caps the cost on a model with many
# flat directions, where the ones with the least curvature are the ones worth
# asking about.
#
# The change is measured as `abs`, not as a drop. A direction along which the
# likelihood *rises* is not flat -- it is a direction the optimiser has not
# finished with -- and on a far-out flat ray the smallest eigenvalue goes
# negative from rounding, so those arrive here as candidates and must not be
# confirmed.
.ctOptimFlatDirectionScreen <- function(info, lpgFunc, est,
  rtol=1e-7, bar=stats::qchisq(0.95, 1) / 2, lengths=c(0.25, 1, 4),
  maxdirections=20L, tolerance=1e-6, corrections=2L){
  if(is.null(info) || is.null(lpgFunc) || is.null(est)) return(NULL)
  if(!is.function(lpgFunc)) return(NULL)
  info <- as.matrix(info)
  n <- nrow(info)
  if(!n || n != ncol(info) || length(est) != n) return(NULL)
  if(!all(is.finite(info))) return(NULL)
  info <- (info + t(info)) / 2
  eig <- try(eigen(info, symmetric=TRUE), silent=TRUE)
  if('try-error' %in% class(eig)) return(NULL)
  values <- eig$values
  scale <- max(values)
  if(!is.finite(scale) || scale <= 0) return(NULL)
  candidates <- which(values <= rtol * scale)
  if(!length(candidates)) return(NULL)
  # Smallest curvature first, so the cap keeps the directions the question is
  # really about.
  candidates <- candidates[order(values[candidates])]
  if(length(candidates) > maxdirections) candidates <- candidates[seq_len(maxdirections)]
  # The value and, where the caller's function carries one, the gradient: the
  # correction below needs it and the evaluation has already paid for it.
  evaluate <- function(x){
    out <- try(lpgFunc(x), silent=TRUE)
    if('try-error' %in% class(out)) return(list(value=NA_real_, gradient=NULL))
    list(value=suppressWarnings(as.numeric(out)[1L]),
      gradient=attr(out, 'gradient'))
  }
  base <- evaluate(est)$value
  if(!isTRUE(is.finite(base))) return(NULL)
  flat <- rep(FALSE, n)
  change <- rep(NA_real_, n)
  # Counted rather than timed. What this costs is a number of likelihood
  # evaluations, which is the same on any machine and under any load; a wall
  # clock here would measure the box. One for the base point, then up to
  # `length(lengths)` per side per candidate plus the corrections, fewer for
  # every side that leaves the ladder early or settles the question.
  evaluations <- 1L
  # The directions a correction may move in: every one the curvature trusts,
  # which excludes the candidate being walked and every other candidate, whose
  # near-zero curvature would turn a gradient into an enormous step.
  trusted <- values > rtol * scale
  basis <- eig$vectors[, trusted, drop=FALSE]
  curvature <- values[trusted]
  # Newton steps from a rung back towards the ridge, each accepted only if it
  # improves, and no further once the drop is inside the bar.
  correct <- function(point, value, gradient){
    used <- 0L
    for(i in seq_len(corrections)){
      if(!ncol(basis) || base - value < bar) break
      g <- suppressWarnings(as.numeric(gradient))
      if(length(g) != n || !all(is.finite(g))) break
      step <- as.numeric(basis %*% (crossprod(basis, g) / curvature))
      if(!all(is.finite(step)) || !any(step != 0)) break
      moved <- FALSE
      for(alpha in c(1, 0.5, 0.25, 0.125)){
        trial <- evaluate(point + alpha * step)
        used <- used + 1L
        if(isTRUE(is.finite(trial$value)) && trial$value > value){
          point <- point + alpha * step
          value <- trial$value
          gradient <- trial$gradient
          moved <- TRUE
          break
        }
      }
      if(!moved) break
    }
    list(point=point, value=value, evaluations=used)
  }
  # The other half of what these evaluations are worth, and it used to be
  # thrown away. A candidate direction along which the likelihood *rises* says
  # the estimate is not a maximum, and the point that proved it is already paid
  # for -- so it is kept, with its sign, rather than collapsed into the
  # magnitude that decides flatness. `_ctsem_overshot` in the engine asks the
  # same question of a different set of directions (magnitude-ordered
  # coordinate prefixes), so this one can find what that one does not look at.
  gain <- 0
  gainpoint <- NULL
  gaindirection <- NA_integer_
  better <- function(value, point, k){
    if(value - base > gain){
      gain <<- value - base
      gainpoint <<- point
      gaindirection <<- k
    }
  }
  for(k in candidates){
    v <- eig$vectors[, k]
    sides <- c(NA_real_, NA_real_)
    for(side in 1:2){
      direction <- c(1, -1)[side]
      worst <- 0
      usable <- TRUE
      for(len in lengths){
        trialpoint <- est + direction * len * v
        trial <- evaluate(trialpoint)
        evaluations <- evaluations + 1L
        # A point the model cannot evaluate is not evidence of flatness, so
        # this side says nothing.
        if(!isTRUE(is.finite(trial$value))){
          usable <- FALSE
          break
        }
        value <- trial$value
        better(value, trialpoint, k)
        if(base - value >= bar){
          followed <- correct(trialpoint, value, trial$gradient)
          evaluations <- evaluations + followed$evaluations
          value <- followed$value
          better(value, followed$point, k)
        }
        worst <- max(worst, abs(value - base))
        # Past the bar this side is refused, and no further displacement can
        # un-refuse it. Worth the early exit rather than completing the
        # ladder: on the laplace route every one of these is an inner mode
        # solve per subject, and the candidates that are *not* flat are exactly
        # the ones a longer walk would spend the most on.
        if(worst >= bar) break
      }
      if(!usable) next
      sides[side] <- worst
      # One flat side settles it; see above.
      if(worst < bar) break
    }
    if(all(is.na(sides))) next
    change[k] <- min(sides, na.rm=TRUE)
    flat[k] <- change[k] < bar
  }
  # `gain` is only reported when it is larger than the optimiser's own
  # convergence tolerance: a rise of 1e-12 along a flat direction is the
  # arithmetic, not a better point.
  found <- is.finite(gain) && gain > tolerance && !is.null(gainpoint)
  list(eig=eig, flat=flat, change=change, bar=bar, rtol=rtol,
    candidates=candidates, lengths=lengths, base=base,
    evaluations=evaluations,
    gain=if(found) gain else 0, point=if(found) gainpoint else NULL,
    direction=if(found) gaindirection else NA_integer_)
}

#
# `eig` and `flat` come from `.ctOptimFlatDirectionScreen()` when the caller
# ran it: the decomposition so it is not taken twice, and a mask of directions
# the likelihood was measured to be flat along. `flat` can only *remove*
# directions from the kept subspace, never add one, so with it absent or all
# FALSE this is exactly the eigenvalue rule it has always been.
.ctOptimIdentifiedInverse <- function(info, rtol=.ctFlatDirectionRtol(),
  eig=NULL, flat=NULL){
  info <- (info + t(info)) / 2
  if(is.null(eig)) eig <- try(eigen(info, symmetric=TRUE), silent=TRUE)
  if('try-error' %in% class(eig)) return(NULL)
  values <- eig$values
  scale <- max(values)
  if(!is.finite(scale) || scale <= 0) return(NULL)
  threshold <- rtol * scale
  keep <- values > threshold
  if(!is.null(flat) && length(flat) == length(keep)) keep <- keep & !flat
  if(!any(keep)) return(NULL)
  vectors <- eig$vectors[, keep, drop=FALSE]
  cov <- vectors %*% (t(vectors) / values[keep])
  cov <- (cov + t(cov)) / 2
  if(any(!is.finite(cov))) return(NULL)
  # How much of each coordinate the projection took away, which is the number
  # that says whether its reported spread means anything: a coordinate with a
  # large share of the dropped subspace has an infinite asymptotic variance,
  # and the covariance built here reports it as almost none. Free from the
  # decomposition already taken, basis-invariant where any single
  # eigenvector's loading is not, and the same quantity
  # `.ctBackendIntervalCheck()` reads -- see the comment there for the fit this
  # was measured on.
  nullvectors <- eig$vectors[, !keep, drop=FALSE]
  mass <- if(ncol(nullvectors)) rowSums(nullvectors^2) else rep(0, nrow(info))
  # Which parameters the dropped directions involve: those with a share of
  # them at `.ctNullMassBar()`, the rule `.ctBackendIdentifiability()` and the
  # interval check name parameters by. It was a loading of 0.25 on any one
  # dropped direction, which named fewer and depended on the basis.
  loaded <- which(mass >= .ctNullMassBar())
  list(cov=cov, nnull=sum(!keep), nullEigenvalues=values[!keep],
    nullParameters=loaded, nullMass=mass, threshold=threshold)
}

ctOptimCovFromHessian <- function(hess, ridge=1e-8, rtol=.ctFlatDirectionRtol(), warn=TRUE,
  context='Hessian', screen=NULL){
  hess <- (hess + t(hess)) / 2
  info <- -hess
  infoEig <- try(eigen(info, symmetric=TRUE, only.values=TRUE), silent=TRUE)
  minInfoEig <- if('try-error' %in% class(infoEig)) NA_real_ else
    min(infoEig$values)
  # The largest as well, because the smallest on its own says nothing. A
  # minimum of -7e-11 is rounding when the largest is 2.6e4 and a real rank
  # deficiency when the largest is 1e-9, and the repair warning below could not
  # tell those apart: it reported a magnitude with no scale to read it against.
  # A warning that fires identically on both is one people learn to ignore,
  # which is the worst outcome, because the text is the same when it mattered.
  maxInfoEig <- if('try-error' %in% class(infoEig)) NA_real_ else
    max(infoEig$values)
  infoEigenRatio <- if(is.finite(maxInfoEig) && maxInfoEig > 0 &&
      is.finite(minInfoEig)) abs(minInfoEig) / maxInfoEig else NA_real_
  covOk <- function(x){
    if('try-error' %in% class(x) || any(!is.finite(x))) return(FALSE)
    cholcheck <- try(suppressWarnings(chol((x + t(x)) / 2)), silent=TRUE)
    !'try-error' %in% class(cholcheck)
  }
  nearPDCov <- function(x){
    out <- try(Matrix::nearPD((x + t(x)) / 2, conv.norm.type='F',
        base.matrix=TRUE)$mat, silent=TRUE)
    if('try-error' %in% class(out)) return(out)
    (out + t(out)) / 2
  }
  repairSteps <- character()
  rawSolveSucceeded <- FALSE
  rawCholSucceeded <- FALSE
  usedNearPD <- FALSE
  usedNullProjection <- FALSE
  nullDirections <- 0L
  nullEigenvalues <- numeric()
  nullParameters <- integer()
  nullMass <- numeric()
  usedGinv <- FALSE
  infoNearPD <- FALSE
  covNearPD <- FALSE
  covRidgeApplied <- FALSE
  covReady <- FALSE
  minInfoEigenFinal <- minInfoEig
  minCovEigenOriginal <- NA_real_
  minCovEigenFinal <- NA_real_

  # Decided before `solve()` is tried, not after it fails.
  #
  # Whether a matrix with a numerically zero eigenvalue makes LAPACK's `solve`
  # give up or merely return an enormous inverse is settled by rounding -- the
  # same rounding that orients the null eigenvector -- so a repair reached only
  # on failure is reached only some of the time. That is the coin flip behind
  # the two runs this was found from: same data, same starting values, the same
  # optimum to eight decimal places, and intervals differing by a factor of a
  # hundred, because one of them fell into this branch and the other did not.
  nullPresent <- is.finite(minInfoEig) && is.finite(maxInfoEig) &&
    maxInfoEig > 0 && minInfoEig <= rtol * maxInfoEig
  # And the same decision on measured rather than approximated evidence. A
  # direction the likelihood is flat along has to be projected out whether or
  # not its eigenvalue has underflowed yet, for the reason the paragraph above
  # gives about `solve()`: otherwise which answer a reader gets is settled by
  # where the optimiser stopped along that direction rather than by the data.
  # See `.ctOptimFlatDirectionScreen()`.
  profileFlat <- if(is.null(screen$flat)) 0L else sum(screen$flat)
  if(profileFlat > 0L) {
    nullPresent <- TRUE
    repairSteps <- c(repairSteps, paste0(profileFlat,
      ' direction(s) measured flat in the likelihood, within ',
      signif(screen$bar, 3), ' log units over displacements of ',
      paste(signif(screen$lengths, 3), collapse='/')))
  }

  rawcov <- if(nullPresent) {
    repairSteps <- c(repairSteps, paste0(
      'solve(-hessian) not attempted: smallest eigenvalue is ',
      signif(infoEigenRatio, 3), ' of the largest'))
    structure('skipped', class='try-error')
  } else try(suppressWarnings(solve(info)), silent=TRUE)
  rawSolveSucceeded <- !'try-error' %in% class(rawcov) &&
    all(is.finite(rawcov))
  if(rawSolveSucceeded) {
    rawcov <- (rawcov + t(rawcov)) / 2
    minCovEigenOriginal <- min(eigen(rawcov, symmetric=TRUE,
      only.values=TRUE)$values)
    rawCholSucceeded <- covOk(rawcov)
    if(rawCholSucceeded) {
      cov <- rawcov
      minCovEigenFinal <- minCovEigenOriginal
      diagnostics <- list(context=context, ridge=ridge,
        method='solve', minInfoEigenOriginal=minInfoEig,
        rawSolveSucceeded=rawSolveSucceeded,
        rawCholSucceeded=rawCholSucceeded,
        infoNearPD=FALSE, usedNullProjection=FALSE,
        nullDirections=0L, nullEigenvalues=numeric(),
        nullParameters=integer(), nullMass=numeric(),
        minInfoEigenFinal=minInfoEigenFinal,
        usedNearPD=FALSE, usedGinv=FALSE,
        covNearPD=FALSE, covRidgeApplied=FALSE,
        minCovEigenOriginal=minCovEigenOriginal,
        minCovEigenFinal=minCovEigenFinal, repairSteps=repairSteps)
      attr(cov, 'ctOptimCovFromHessian') <- diagnostics
      return(cov)
    }
    repairSteps <- c(repairSteps,
      'solve(-hessian) succeeded but covariance was not positive definite')
    npdcov <- nearPDCov(rawcov)
    if(covOk(npdcov)) {
      cov <- npdcov
      covReady <- TRUE
      usedNearPD <- TRUE
      covNearPD <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps, 'nearPD applied to solved covariance')
    }
  } else if(!nullPresent) {
    repairSteps <- c(repairSteps, 'solve(-hessian) failed')
  }

  if(!covReady) {
    # See `.ctOptimIdentifiedInverse()`. This used to floor the information
    # eigenvalues at `ridge` and invert, which put 1/ridge along every
    # direction the data does not determine and leaked it into the ones it
    # does. `covOk()` is deliberately not the test here: the projected
    # covariance is singular by construction -- that is what it is for -- so
    # `chol()` cannot succeed on it and asking would send every such matrix to
    # the generalized inverse below.
    projected <- .ctOptimIdentifiedInverse(info, rtol=rtol,
      eig=screen$eig, flat=screen$flat)
    if(!is.null(projected)) {
      cov <- projected$cov
      covReady <- TRUE
      usedNullProjection <- TRUE
      nullDirections <- projected$nnull
      nullEigenvalues <- projected$nullEigenvalues
      nullParameters <- projected$nullParameters
      nullMass <- projected$nullMass
      minInfoEigenFinal <- projected$threshold
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps, paste0(nullDirections,
        ' direction(s) with no curvature were projected out before inversion'))
    } else {
      repairSteps <- c(repairSteps,
        'the information matrix has no direction with positive curvature')
    }
  }
  
  if(!covReady) {
    usedGinv <- TRUE
    ginvcov <- try(MASS::ginv(info), silent=TRUE)
    if(!'try-error' %in% class(ginvcov)) ginvcov <- (ginvcov + t(ginvcov)) / 2
    if(covOk(ginvcov)) {
      cov <- ginvcov
      covReady <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps, 'MASS::ginv(-hessian) used')
    } else {
      npdcov <- if('try-error' %in% class(ginvcov)) ginvcov else
        nearPDCov(ginvcov)
      if(!covOk(npdcov)) {
        if('try-error' %in% class(ginvcov)) {
          stop('Could not construct covariance from Hessian using solve, ridge repair, or generalized inverse.',
            call.=FALSE)
        }
        npdcov <- ctOptimSafeCov(ginvcov, ridge=ridge)
        covRidgeApplied <- isTRUE(attr(npdcov, 'ctOptimSafeCov')$ridgeApplied)
      }
      cov <- npdcov
      covReady <- TRUE
      covNearPD <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps,
        'MASS::ginv(-hessian) used with covariance positive-definite cleanup')
    }
  }
  
  if(is.na(minCovEigenOriginal) && covReady) {
    minCovEigenOriginal <- min(eigen((cov + t(cov)) / 2, symmetric=TRUE,
      only.values=TRUE)$values)
  }
  diagnostics <- list(context=context, ridge=ridge,
    minInfoEigenOriginal=minInfoEig,
    maxInfoEigenOriginal=maxInfoEig,
    infoEigenRatio=infoEigenRatio,
    # `sqrt(eps)` is where a symmetric eigendecomposition stops being able to
    # tell a small eigenvalue from zero, so below it the repair is arithmetic
    # rather than a statement about the model.
    #
    # The sign is half the test, and leaving it out was wrong. Rounding error
    # takes an eigenvalue that should be positive and makes it slightly
    # *negative*; it does not make it exactly zero. So a zero -- or a positive
    # value small enough to have needed the ridge -- is a genuinely singular
    # direction, which is the most serious case rather than the least, and the
    # ratio test alone classified it as negligible. `-diag(c(1, 0))` is the
    # minimal example, and it stopped warning.
    infoRepairNegligible=is.finite(infoEigenRatio) && minInfoEig < 0 &&
      infoEigenRatio < sqrt(.Machine$double.eps),
    rawSolveSucceeded=rawSolveSucceeded,
    rawCholSucceeded=rawCholSucceeded,
    infoNearPD=infoNearPD,
    usedNullProjection=usedNullProjection,
    # The directions the data does not determine: how many, how flat, and which
    # parameters carry them. Reported rather than repaired away, because a
    # covariance that is silently missing a dimension is the thing this used to
    # hide behind a fabricated 1e4 standard error.
    nullDirections=nullDirections,
    nullEigenvalues=nullEigenvalues,
    nullParameters=nullParameters,
    # Per coordinate, so a caller can ask which reported spreads the
    # projection removed rather than only how many directions it dropped.
    nullMass=nullMass,
    # How many of those directions were dropped because the likelihood was
    # measured flat along them rather than because their eigenvalue had
    # underflowed, and by how much the likelihood moved when they were walked.
    # Separated because they are different evidence: the first is a statement
    # about the data, the second about arithmetic.
    profileFlatDirections=profileFlat,
    profileChange=if(is.null(screen$flat)) numeric() else
      screen$change[screen$flat],
    profileBar=if(is.null(screen)) NA_real_ else screen$bar,
    minInfoEigenFinal=minInfoEigenFinal,
    usedNearPD=usedNearPD,
    usedGinv=usedGinv,
    covNearPD=covNearPD,
    covRidgeApplied=covRidgeApplied,
    minCovEigenOriginal=minCovEigenOriginal,
    minCovEigenFinal=minCovEigenFinal,
    method=if(usedGinv) 'ginv' else if(usedNullProjection) 'nullprojection'
      else if(usedNearPD) 'nearPD_cov' else 'solve',
    repairSteps=repairSteps)
  attr(cov, 'ctOptimCovFromHessian') <- diagnostics
  # `repairSteps` is the audit trail and goes out whole on the diagnostics; the
  # warning gets one sentence per fact. Three of those steps describe the null
  # projection from three angles -- solve was skipped, directions were dropped,
  # they have no spread -- and printing all three spent most of a warning
  # saying one thing. R truncates a warning at `getOption('warning.length')`,
  # 1000 bytes by default, so the length was not merely untidy: the tail of
  # this warning and of the identifiability one was being cut off mid-word.
  issues <- if(isTRUE(diagnostics$usedNullProjection)) {
    repairSteps[!grepl('^solve\\(-hessian\\) not attempted', repairSteps) &
      !grepl('projected out before inversion$', repairSteps)]
  } else repairSteps
  if(isTRUE(diagnostics$infoNearPD)) issues <- c(issues,
    'nearPD was needed for the information matrix')
  if(isTRUE(diagnostics$usedNullProjection)) issues <- c(issues,
    paste0(nullDirections, ' direction(s) with no curvature left out of the',
      ' inversion',
      if(is.finite(infoEigenRatio))
        paste0(' (smallest eigenvalue ', signif(infoEigenRatio, 3),
          ' of the largest)') else '',
      # The count of *parameters* as well as of directions, because that is the
      # number a reader is about to be misled by. A parameter lying along a
      # dropped direction has an infinite variance and is given a small
      # reported one, which reads as precision rather than as a gap. See
      # `.ctBackendIntervalCheck()`, which names them on a julia fit.
      if(sum(nullMass >= .ctNullMassBar()) > 0) paste0(', and ',
        sum(nullMass >= .ctNullMassBar()),
        ' parameter(s) along them whose reported sd is that projection rather',
        ' than a small width') else ''))
  if(isTRUE(diagnostics$usedGinv)) issues <- c(issues,
    'MASS::ginv() was used')
  if(isTRUE(diagnostics$covNearPD) || isTRUE(diagnostics$covRidgeApplied)) {
    issues <- c(issues,
      'the resulting covariance required positive-definite cleanup')
  }
  if(warn && length(issues) > 0) {
    # Graded. A repair that floored an eigenvalue indistinguishable from zero is
    # arithmetic, and is reported as such; one that floored a substantively
    # negative or tiny eigenvalue is a statement about what the data can
    # determine, and keeps the warning. Making that distinction is the point --
    # the same text for both is what taught people to ignore it, and the cases
    # it currently conflates are genuinely different. A near-integrated trend
    # process *should* warn here.
    #
    # `usedNullProjection` is excluded from the quiet branch on purpose. A
    # direction that had to be left out of the inversion is a statement about
    # what the data determines whatever the sign of the smallest eigenvalue
    # was, and the parameters carrying it have no reported spread at all.
    if(isTRUE(diagnostics$infoRepairNegligible) &&
        !isTRUE(diagnostics$usedNullProjection) &&
        !isTRUE(diagnostics$usedGinv) && !isTRUE(diagnostics$infoNearPD)) {
      # Deliberately not "arithmetic, not a statement about the model", which
      # is what this said and could not support. A smallest eigenvalue of
      # -1e-12 relative to the largest is what rounding does to a
      # positive-definite matrix, and it is *also* what rounding does to a
      # genuinely singular one -- the two are indistinguishable from this
      # number alone. Usually the first, so a message rather than a warning;
      # never certainly the first, so it says which is which is checkable
      # elsewhere rather than pronouncing.
      message(context, ' covariance: the information matrix needed a numerical ',
        'nudge before inversion (smallest eigenvalue ',
        signif(minInfoEig, 3), ', ', signif(infoEigenRatio, 3),
        ' of the largest, which is indistinguishable from zero at machine ',
        'precision). Usually rounding rather than a flat direction, but the ',
        'two look the same at this magnitude; fit$identifiability names the ',
        'parameters if it is the latter.')
    } else {
      # Classed, so that a julia fit can hold it to its closing summary
      # (`.ctBackendFitWarnings()`), where the identifiability warning usually
      # says the same thing in terms a reader can act on.
      warning(warningCondition(paste0(context,
        ' covariance from Hessian required numerical repair: ',
        paste(issues, collapse='; '),
        if(is.finite(infoEigenRatio) &&
            infoEigenRatio >= sqrt(.Machine$double.eps))
          paste0('. The smallest eigenvalue is ', signif(infoEigenRatio, 3),
            ' of the largest, too large to be rounding: some direction of this ',
            'model is close to unidentified and the standard errors along it ',
            'are not trustworthy') else ''),
        class = 'ctsemCovarianceRepair'))
    }
  }
  cov
}

ctOptimNormalDraws <- function(mean, cov, n, df=Inf){
  cov <- ctOptimSafeCov(cov)
  z <- matrix(stats::rnorm(n * length(mean)), nrow=n)
  draws <- z %*% chol(cov)
  if(is.finite(df)){
    draws <- draws / sqrt(stats::rchisq(n, df=df) / df)
  }
  sweep(draws, 2, mean, '+')
}

ctOptimScoreMatrix <- function(standata, sm, est, cores=1, scores=NULL){
  if(!is.null(scores)) return(scores)
  scorecalc(standata=standata, est=est, stanmodel=sm,
    subjectsonly=ctOptimNSubjects(standata) >= 2,
    returnsubjectlist=FALSE,
    cores=cores)
}

ctOptimNSubjects <- function(standata){
  nsubjects <- suppressWarnings(as.integer(standata$nsubjects[1]))
  if(length(nsubjects) < 1 || is.na(nsubjects) || !is.finite(nsubjects)) {
    if(!is.null(standata$subject)) nsubjects <- length(unique(standata$subject))
  }
  if(length(nsubjects) < 1 || is.na(nsubjects) || !is.finite(nsubjects)) {
    nsubjects <- NA_integer_
  }
  nsubjects
}

ctOptimNDataPoints <- function(standata){
  ndatapoints <- suppressWarnings(as.integer(standata$ndatapoints[1]))
  if(length(ndatapoints) < 1 || is.na(ndatapoints) || !is.finite(ndatapoints)) {
    if(!is.null(standata$subject)) ndatapoints <- length(standata$subject)
  }
  if(length(ndatapoints) < 1 || is.na(ndatapoints) || !is.finite(ndatapoints)) {
    ndatapoints <- NA_integer_
  }
  ndatapoints
}

ctOptimCheckUncertaintyData <- function(standata, uncertainty, finishsamples,
  npars=NULL){
  nsubjects <- ctOptimNSubjects(standata)
  ndatapoints <- ctOptimNDataPoints(standata)
  if(is.null(npars)) npars <- NA_integer_
  npars <- suppressWarnings(as.integer(npars[1]))
  if(length(npars) < 1 || is.na(npars) || !is.finite(npars)) npars <- NA_integer_
  if(uncertainty %in% c('bootstrap','fullbootstrap') &&
      finishsamples < 2) {
    stop(uncertainty, ' uncertainty requires at least two samples / refits ',
      'to estimate a covariance; increase finishsamples.', call.=FALSE)
  }
  if(uncertainty == 'fullbootstrap'){
    if(is.na(nsubjects) || nsubjects < 2) {
      stop('fullbootstrap uncertainty requires at least two subjects.',
        call.=FALSE)
    }
    if(nsubjects < 10) {
      warning('fullbootstrap uncertainty requested with fewer than ten ',
        'independent subjects; the bootstrap distribution may be unstable.',
        call.=FALSE)
    }
    if(!is.na(npars) && finishsamples <= npars) {
      warning('fullbootstrap requested with finishsamples <= number of raw ',
        'parameters; the empirical covariance is rank limited and will be ',
        'regularised.', call.=FALSE)
    }
  }
  if(uncertainty %in% c('bootstrap','sandwich','opg')){
    subjectScores <- !is.na(nsubjects) && nsubjects >= 2
    nscore <- if(subjectScores) nsubjects else ndatapoints
    if(is.na(nscore) || nscore < 2) {
      stop(uncertainty, ' uncertainty requires at least two ',
        if(subjectScores) 'subject-level' else 'case-level',
        ' score contribution rows.', call.=FALSE)
    }
    if(subjectScores && nscore < 10) {
      warning(uncertainty, ' uncertainty is based on fewer than ten ',
        'independent subject-level score contributions; the covariance ',
        'estimate may be unstable.', call.=FALSE)
    }
    if(!is.na(npars) && nscore <= npars) {
      warning(uncertainty, ' uncertainty has no more score contribution rows ',
        'than raw parameters; the covariance estimate is rank limited and ',
        'will be regularised.', call.=FALSE)
    }
    if(!subjectScores) {
      warning(uncertainty, ' uncertainty for a single-subject model uses ',
        'case-level score contributions. This can be unreliable when ',
        'observations are serially dependent.',
        call.=FALSE)
    }
  }
  invisible(list(nsubjects=nsubjects, ndatapoints=ndatapoints))
}

ctOptimBootstrapDraws <- function(est, cov, scores, n=1000){
  scores <- as.matrix(scores)
  scores <- scale(scores, center=TRUE, scale=FALSE)
  nscore <- nrow(scores)
  if(nscore < 2) stop('At least two score rows are required for bootstrap uncertainty')
  draws <- matrix(NA_real_, nrow=n, ncol=length(est))
  for(i in seq_len(n)){
    idx <- sample.int(nscore, nscore, replace=TRUE)
    score_sum <- colSums(scores[idx,,drop=FALSE])
    draws[i,] <- est + as.numeric(cov %*% score_sum)
  }
  draws
}

ctOptimBootstrapStandata <- function(standata, subjects){
  subjects <- as.integer(subjects)
  if(length(subjects) < 2) stop('At least two subjects are required for full bootstrap uncertainty')
  long <- standatatolong(standata)
  longlist <- vector('list', length(subjects))
  for(i in seq_along(subjects)){
    longi <- long[long$subject %in% subjects[i], , drop=FALSE]
    if(nrow(longi) < 1) stop('Subject ', subjects[i], ' not found in standata')
    longi$subject <- i
    longlist[[i]] <- longi
  }
  longboot <- do.call(rbind, longlist)
  row.names(longboot) <- NULL
  standataboot <- standatalongremerge(long=longboot, standata=standata)
  standataboot$ndatapoints <- as.integer(nrow(longboot))
  standataboot$nsubjects <- as.integer(length(subjects))
  standataboot$subject <- array(as.integer(longboot$subject))
  if(standata$ntipred > 0) {
    standataboot$tipredsdata <- standata$tipredsdata[subjects, , drop=FALSE]
  }
  standataboot$idmap <- data.frame(
    original=paste0('boot', seq_along(subjects), '_subject', subjects),
    new=seq_along(subjects))
  standataboot
}

ctOptimFullBootstrapOne <- function(i, est, standata, sm, fitCores, tol,
  verbose=0){
  subjects <- unique(standata$subject)
  sampledSubjects <- sample(subjects, length(subjects), replace=TRUE)
  standataboot <- ctOptimBootstrapStandata(standata=standata,
    subjects=sampledSubjects)
  standataboot$savesubjectmatrices <- 0L
  standataboot$nsubsets <- 1L
  lpgsetup <- ctOptimDataLpgFunc(sm=sm, standata=standataboot,
    cores=fitCores)
  on.exit({
    if(!is.null(lpgsetup$cl)) try(parallel::stopCluster(lpgsetup$cl),
      silent=TRUE)
    if(!is.null(lpgsetup$smfile) && nzchar(lpgsetup$smfile)) {
      try(file.remove(lpgsetup$smfile), silent=TRUE)
    }
  }, add=TRUE)
  if(verbose > 0) message('Full bootstrap sample ', i)
  opt <- try(ctOptim(init=est, lpgFunc=lpgsetup$lpg, tol=tol,
    nsubsets=1L, stochastic=FALSE, stochasticTolAdjust=1,
    bfgsType='mize'), silent=TRUE)
  if('try-error' %in% class(opt) || is.null(opt$par) ||
      length(opt$par) != length(est) || any(!is.finite(opt$par))) {
    msg <- if('try-error' %in% class(opt) && !is.null(attr(opt,
          'condition'))) {
      conditionMessage(attr(opt, 'condition'))
    } else 'non-finite optimized parameters'
    return(list(ok=FALSE, par=rep(NA_real_, length(est)),
      sampledSubjects=sampledSubjects, message=msg))
  }
  list(ok=TRUE, par=opt$par, sampledSubjects=sampledSubjects,
    value=opt$value, message=opt$message)
}

ctOptimFullBootstrapDraws <- function(est, standata, sm, n=1000, cores=1,
  control=list(), verbose=0){
  if(standata$nsubjects < 2) {
    stop('Full bootstrap uncertainty requires at least two subjects')
  }
  if(is.null(control$bootstrapFitCores)) control$bootstrapFitCores <- 1L
  if(is.null(control$bootstrapTol)) control$bootstrapTol <- 1e-5
  fitCores <- suppressWarnings(as.integer(control$bootstrapFitCores[1]))
  if(!is.finite(fitCores) || is.na(fitCores) || fitCores < 1) fitCores <- 1L
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  outerCores <- min(n, max(1L, floor(cores / fitCores)))
  fitCores <- min(fitCores, standata$nsubjects)
  if(verbose > 0 || outerCores > 1) {
    message('Fitting ', n, ' full bootstrap samples using ', outerCores,
      ' bootstrap worker(s) and ', fitCores, ' core(s) per refit')
  }
  
  if(outerCores > 1){
    cl <- makeClusterID(outerCores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
    bootHelpers <- c('ctOptimFullBootstrapOne',
      'ctOptimBootstrapStandata', 'ctOptimDataLpgFunc', 'ctOptim',
      'standatatolong', 'standatalongremerge', 'standatalongobjects',
      'stan_reinitsf', 'getcxxfun', 'suppressOutput', 'makeClusterID',
      'parallelStanSetup', 'clusterIDexport', 'clusterIDeval',
      'singlecoreStanSetup', 'parlptext')
    parallel::clusterExport(cl, c('bootHelpers', bootHelpers),
      envir=environment(ctOptimFullBootstrapDraws))
    parallel::clusterEvalQ(cl, {
      bootEnv <- new.env(parent=.GlobalEnv)
      for(fn in bootHelpers){
        obj <- get(fn, envir=.GlobalEnv)
        if(is.function(obj)){
          environment(obj) <- bootEnv
        }
        assign(fn, obj, envir=bootEnv)
      }
      options(ctsem.bootstrap.env=bootEnv)
      rm(list=c(bootHelpers, 'bootHelpers'), envir=.GlobalEnv)
      rm(bootHelpers, bootEnv)
      NULL
    })
    out <- parallel::parLapplyLB(cl, seq_len(n), function(i){
      bootEnv <- getOption('ctsem.bootstrap.env')
      if(is.null(bootEnv)) stop('Missing ctsem bootstrap worker environment')
      get('ctOptimFullBootstrapOne', envir=bootEnv)(i=i, est=est, standata=standata,
        sm=sm, fitCores=fitCores, tol=control$bootstrapTol,
        verbose=verbose)
    })
  } else {
    out <- lapply(seq_len(n), function(i){
      if(verbose == 0) message('\rFull bootstrap sample ', i, '/', n,
        appendLF=FALSE)
      ctOptimFullBootstrapOne(i=i, est=est, standata=standata, sm=sm,
        fitCores=fitCores, tol=control$bootstrapTol, verbose=verbose)
    })
    if(verbose == 0) message('')
  }
  ok <- vapply(out, `[[`, logical(1), 'ok')
  if(!any(ok)) {
    msgs <- unique(vapply(out, function(x) x$message, character(1)))
    stop('All full bootstrap refits failed. First errors: ',
      paste(utils::head(msgs, 3), collapse='; '))
  }
  if(any(!ok)) {
    warning(sum(!ok), ' full bootstrap refits failed and were omitted.',
      call.=FALSE)
  }
  draws <- do.call(rbind, lapply(out[ok], `[[`, 'par'))
  list(draws=draws,
    sampledSubjects=lapply(out[ok], `[[`, 'sampledSubjects'),
    failures=out[!ok], outerCores=outerCores, fitCores=fitCores,
    tol=control$bootstrapTol)
}

ctOptimSurrogateDirections <- function(p, n){
  dirs <- rbind(diag(p), -diag(p))
  while(nrow(dirs) < n){
    addn <- ceiling((n - nrow(dirs)) / 2)
    z <- matrix(stats::rnorm(addn * p), nrow=addn)
    z <- z / sqrt(rowSums(z^2))
    dirs <- rbind(dirs, z, -z)
  }
  dirs[seq_len(n),,drop=FALSE]
}

ctOptimSurrogateEvalPoint <- function(est, lpgFunc, cholcov, z, baseValue){
  rawstep <- as.numeric(z %*% cholcov)
  lp <- try(suppressMessages(suppressWarnings(lpgFunc(est + rawstep))),
    silent=TRUE)
  value <- NA_real_
  grad <- rep(NA_real_, length(est))
  if(!'try-error' %in% class(lp)) {
    value <- lp[1]
    lpgrad <- attributes(lp)$gradient
    if(length(lpgrad) == length(est)) grad <- lpgrad
  }
  drop <- baseValue - value
  finite <- is.finite(value) && is.finite(drop) && all(is.finite(grad))
  list(white=as.numeric(z), raw=rawstep, value=value, gradient=grad,
    drop=drop, finite=finite)
}

ctOptimSurrogateBestEval <- function(evals, targetDrop, preferRange=NULL){
  drops <- vapply(evals, `[[`, numeric(1), 'drop')
  finite <- vapply(evals, `[[`, logical(1), 'finite') &
    is.finite(drops) & drops > 0
  if(!is.null(preferRange)) {
    inrange <- finite & drops >= preferRange[1] & drops <= preferRange[2]
    if(any(inrange)) finite <- inrange
  }
  if(!any(finite)) {
    finite <- vapply(evals, `[[`, logical(1), 'finite')
  }
  if(!any(finite)) return(evals[[length(evals)]])
  ii <- which(finite)
  ii <- ii[which.min(abs(log(pmax(drops[ii], .Machine$double.eps) /
      targetDrop)))]
  evals[[ii]]
}

ctOptimSurrogateRadiusForDrop <- function(radius, drop, targetDrop,
  minFactor=.2, maxFactor=8, safety=1.1){
  if(!is.finite(drop) || drop <= 0) return(radius * maxFactor)
  factor <- sqrt(targetDrop / drop) * safety
  radius * min(maxFactor, max(minFactor, factor))
}

ctOptimSurrogateTargetPoint <- function(est, lpgFunc, cholcov, direction,
  targetDrop, dropRange, initialRadius, maxRadius=64, maxEval=6,
  targetFactor=1.35, baseValue){
  direction <- as.numeric(direction)
  direction <- direction / sqrt(sum(direction^2))
  evals <- list()
  evalAt <- function(radius){
    ev <- ctOptimSurrogateEvalPoint(est=est, lpgFunc=lpgFunc,
      cholcov=cholcov, z=radius * direction, baseValue=baseValue)
    ev$radius <- radius
    ev
  }
  addEval <- function(radius){
    evals[[length(evals) + 1L]] <<- evalAt(radius)
    evals[[length(evals)]]
  }
  goodTarget <- function(ev){
    ev$finite && is.finite(ev$drop) && ev$drop > 0 &&
      abs(log(ev$drop / targetDrop)) <= log(targetFactor)
  }
  initialRadius <- max(.Machine$double.eps, initialRadius)
  ev <- addEval(initialRadius)
  if(goodTarget(ev)) {
    ev$targeted <- TRUE
    ev$neval <- length(evals)
    return(ev)
  }
  
  low <- 0
  high <- initialRadius
  if(ev$finite && is.finite(ev$drop) && ev$drop > 0 &&
      ev$drop < targetDrop) {
    low <- initialRadius
    while(length(evals) < maxEval && high < maxRadius) {
      high <- min(maxRadius, ctOptimSurrogateRadiusForDrop(high,
        ev$drop, targetDrop, minFactor=1.5))
      ev <- addEval(high)
      if(goodTarget(ev)) {
        ev$targeted <- TRUE
        ev$neval <- length(evals)
        return(ev)
      }
      if(!ev$finite || !is.finite(ev$drop) || ev$drop >= targetDrop) break
      low <- high
    }
  }
  
  while(length(evals) < maxEval && high > 0 && high > low) {
    if(is.finite(ev$drop) && ev$drop > targetDrop && low == 0) {
      mid <- ctOptimSurrogateRadiusForDrop(high, ev$drop, targetDrop,
        minFactor=.1, maxFactor=.8, safety=.9)
      mid <- min(high * .95, max(.Machine$double.eps, mid))
    } else {
      mid <- (low + high) / 2
    }
    ev <- addEval(mid)
    if(goodTarget(ev)) {
      ev$targeted <- TRUE
      ev$neval <- length(evals)
      return(ev)
    }
    if(!ev$finite || !is.finite(ev$drop) || ev$drop >= targetDrop) {
      high <- mid
    } else {
      low <- mid
    }
  }
  
  ev <- ctOptimSurrogateBestEval(evals=evals, targetDrop=targetDrop,
    preferRange=dropRange)
  ev$targeted <- ev$finite && is.finite(ev$drop) && ev$drop >= dropRange[1] &&
    ev$drop <= dropRange[2]
  ev$neval <- length(evals)
  ev
}

ctOptimSurrogateBacktransformHessian <- function(hessWhite, cholcov){
  invchol <- backsolve(cholcov, diag(ncol(cholcov)))
  hess <- invchol %*% hessWhite %*% t(invchol)
  (hess + t(hess)) / 2
}

ctOptimSurrogateProfileDirections <- function(est, lpgFunc, cholcov,
  directions, targetDrop=2, maxStep=64, tol=.02, maxIter=25,
  initialStep=1, maxExpand=4, baseValue=NULL, verbose=0){
  if(is.null(baseValue)) {
    baseValue <- suppressMessages(suppressWarnings(lpgFunc(est)))[1]
  }
  directions <- as.matrix(directions)
  if(nrow(directions) < 1) {
    return(data.frame())
  }
  evalDrop <- function(z){
    rawstep <- as.numeric(z %*% cholcov)
    lp <- try(suppressMessages(suppressWarnings(lpgFunc(est + rawstep))),
      silent=TRUE)
    if('try-error' %in% class(lp)) return(NA_real_)
    baseValue - lp[1]
  }
  out <- vector('list', nrow(directions) * 2L)
  oi <- 0L
  for(di in seq_len(nrow(directions))){
    diri <- directions[di,]
    diri <- diri / sqrt(sum(diri^2))
    maxStepi <- maxStep[min(length(maxStep), di)]
    initialStepi <- initialStep[min(length(initialStep), di)]
    if(!is.finite(maxStepi) || maxStepi <= 0) maxStepi <- 64
    if(!is.finite(initialStepi) || initialStepi <= 0) initialStepi <- 1
    initialStepi <- min(initialStepi, maxStepi)
    for(sgn in c(-1, 1)){
      oi <- oi + 1L
      low <- 0
      high <- initialStepi
      dhigh <- evalDrop(sgn * high * diri)
      while(is.finite(dhigh) && dhigh < targetDrop && high < maxStepi){
        low <- high
        high <- min(maxStepi, ctOptimSurrogateRadiusForDrop(high,
          dhigh, targetDrop, minFactor=1.5))
        dhigh <- evalDrop(sgn * high * diri)
      }
      expand <- 0L
      while(is.finite(dhigh) && dhigh < targetDrop && expand < maxExpand) {
        expand <- expand + 1L
        low <- high
        high2 <- ctOptimSurrogateRadiusForDrop(high, dhigh, targetDrop,
          minFactor=1.5, maxFactor=16)
        if(!is.finite(high2) || high2 <= high) high2 <- high * 2
        high <- high2
        dhigh <- evalDrop(sgn * high * diri)
      }
      reached <- is.finite(dhigh) && dhigh >= targetDrop
      if(reached) {
        for(iter in seq_len(maxIter)){
          if(low == 0 && is.finite(dhigh) && dhigh > targetDrop) {
            mid <- ctOptimSurrogateRadiusForDrop(high, dhigh, targetDrop,
              minFactor=.1, maxFactor=.8, safety=.9)
            mid <- min(high * .95, max(.Machine$double.eps, mid))
          } else {
            mid <- (low + high) / 2
          }
          dmid <- evalDrop(sgn * mid * diri)
          if(!is.finite(dmid) || dmid >= targetDrop) high <- mid else low <- mid
          if(is.finite(dmid) && abs(dmid - targetDrop) <= tol * targetDrop) break
        }
        step <- high
        drop <- evalDrop(sgn * step * diri)
        curvature <- 2 * targetDrop / step^2
      } else {
        step <- high
        drop <- dhigh
        curvature <- if(is.finite(drop) && drop > 0) 2 * drop / step^2
          else NA_real_
      }
      status <- if(reached) 'reached' else if(!is.finite(drop)) {
        'nonfinite'
      } else if(drop < targetDrop) {
        'belowTargetAtMaxStep'
      } else {
        'notReached'
      }
      out[[oi]] <- data.frame(direction=di, sign=sgn, step=step,
        drop=drop, curvature=curvature, reached=reached, status=status,
        expansions=expand)
    }
  }
  profiles <- do.call(rbind, out)
  if(verbose > 0) {
    message('Surrogate profiled ', nrow(directions), ' direction(s); ',
      sum(profiles$reached), '/', nrow(profiles),
      ' one-sided targets reached')
    if(any(!profiles$reached)) {
      message('Surrogate profile misses: ',
        paste(names(table(profiles$status[!profiles$reached])),
          as.integer(table(profiles$status[!profiles$reached])),
          sep='=', collapse=', '))
    }
  }
  profiles
}

ctOptimSurrogateProfileCurvature <- function(hessWhite, est, lpgFunc,
  cholcov, targetDrop=2, maxStep=64, baseValue=NULL, verbose=0){
  infoWhite <- -((hessWhite + t(hessWhite)) / 2)
  eig <- eigen(infoWhite, symmetric=TRUE)
  directions <- t(eig$vectors)
  expectedStep <- rep(maxStep, length(eig$values))
  positive <- is.finite(eig$values) & eig$values > 0
  expectedStep[positive] <- sqrt(2 * targetDrop / eig$values[positive])
  expectedStep[!is.finite(expectedStep) | expectedStep <= 0] <- maxStep
  profileMaxStep <- pmax(maxStep, expectedStep * 2)
  profileInitialStep <- pmax(.Machine$double.eps, expectedStep)
  profiles <- ctOptimSurrogateProfileDirections(est=est, lpgFunc=lpgFunc,
    cholcov=cholcov, directions=directions, targetDrop=targetDrop,
    maxStep=profileMaxStep, initialStep=profileInitialStep,
    baseValue=baseValue, verbose=verbose)
  adjusted <- 0L
  newvals <- eig$values
  for(i in seq_along(newvals)){
    curv <- profiles$curvature[profiles$direction == i &
        is.finite(profiles$curvature)]
    if(length(curv) > 0) {
      profileCurv <- max(curv)
      if(is.finite(profileCurv) && profileCurv > newvals[i]) {
        newvals[i] <- profileCurv
        adjusted <- adjusted + 1L
      }
    }
  }
  infoWhite <- eig$vectors %*% diag(newvals, length(newvals)) %*%
    t(eig$vectors)
  infoWhite <- (infoWhite + t(infoWhite)) / 2
  list(hessWhite=-infoWhite, profiles=profiles, nProfiled=length(newvals),
    nAdjusted=adjusted)
}

ctOptimSurrogateHessian <- function(est, lpgFunc, cov, npoints=NULL,
  scale=.5, ridge=1e-6, profile=TRUE, profileTargetDrop=NULL,
  profileMaxStep=64, verbose=0){
  p <- length(est)
  if(is.null(npoints)) npoints <- max(4 * p, 50)
  cov <- ctOptimSafeCov(cov, ridge=ridge)
  cholcov <- chol(cov)
  targetDrop <- 2
  dropRange <- c(.25, 6)
  defaultScale <- .5
  initialRadius <- (scale / defaultScale) * sqrt(2 * targetDrop)
  base <- suppressMessages(suppressWarnings(lpgFunc(est)))
  baseValue <- base[1]
  basegrad <- attributes(base)$gradient
  if(is.null(basegrad) || any(!is.finite(basegrad))) basegrad <- rep(0, p)
  
  directions <- ctOptimSurrogateDirections(p=p, n=npoints)
  evals <- vector('list', npoints)
  for(i in seq_len(npoints)){
    if(verbose > 0) message('\rFitting local quadratic surrogate, point ',
      i, '/', npoints, appendLF=FALSE)
    evals[[i]] <- ctOptimSurrogateTargetPoint(est=est, lpgFunc=lpgFunc,
      cholcov=cholcov, direction=directions[i,], targetDrop=targetDrop,
      dropRange=dropRange, initialRadius=initialRadius,
      baseValue=baseValue)
  }
  if(verbose > 0) message('')
  values <- vapply(evals, `[[`, numeric(1), 'value')
  drops <- vapply(evals, `[[`, numeric(1), 'drop')
  targeted <- vapply(evals, `[[`, logical(1), 'targeted')
  neval <- vapply(evals, `[[`, numeric(1), 'neval')
  gradients <- do.call(rbind, lapply(evals, `[[`, 'gradient'))
  design <- do.call(rbind, lapply(evals, `[[`, 'raw'))
  whiteDesign <- do.call(rbind, lapply(evals, `[[`, 'white'))
  finite <- is.finite(values) & is.finite(drops) &
    apply(gradients, 1, function(x) all(is.finite(x)))
  inrange <- finite & drops >= dropRange[1] & drops <= dropRange[2]
  positive <- finite & drops > 0
  keep <- inrange
  if(sum(keep) < npoints && any(positive & !keep)) {
    add <- which(positive & !keep)
    add <- add[order(abs(log(pmax(drops[add], .Machine$double.eps) /
        targetDrop)))]
    keep[add[seq_len(min(length(add), npoints - sum(keep)))]] <- TRUE
  }
  if(sum(keep) <= p) stop('Too few finite surrogate evaluations')
  design <- design[keep,,drop=FALSE]
  whiteDesign <- whiteDesign[keep,,drop=FALSE]
  gradients <- gradients[keep,,drop=FALSE]
  values <- values[keep]
  drops <- drops[keep]
  diagnostics <- list(nRequested=npoints, nUsed=nrow(design),
    nTargeted=sum(targeted), nInRange=sum(inrange),
    nFinitePositive=sum(positive), nFinite=sum(finite),
    meanEvaluations=mean(neval), maxEvaluations=max(neval),
    targetDrop=targetDrop, dropRange=dropRange)
  rawgrad <- sweep(gradients, 2, basegrad, '-')
  y <- rawgrad %*% t(cholcov)
  pointWeights <- 1 / pmax(abs(log(pmax(drops, .Machine$double.eps) /
        targetDrop)), .25)
  pointWeights <- pointWeights / mean(pointWeights, na.rm=TRUE)
  sqrtw <- sqrt(pointWeights)
  xw <- whiteDesign * sqrtw
  yw <- y * sqrtw
  xtx <- crossprod(xw) + diag(ridge, p)
  coef <- solve(xtx, crossprod(xw, yw))
  hessWhite <- (t(coef) + coef) / 2
  if(isTRUE(profile)) {
    profiled <- ctOptimSurrogateProfileCurvature(hessWhite=hessWhite,
      est=est, lpgFunc=lpgFunc, cholcov=cholcov,
      targetDrop=if(is.null(profileTargetDrop)) targetDrop else profileTargetDrop,
      maxStep=profileMaxStep, baseValue=baseValue, verbose=verbose)
    hessWhite <- profiled$hessWhite
  } else {
    profiled <- list(profiles=data.frame(), nProfiled=0L, nAdjusted=0L,
      note='Directional profiling disabled')
  }
  hess <- ctOptimSurrogateBacktransformHessian(hessWhite, cholcov)
  list(hessian=hess, values=values, gradients=gradients, design=design,
    whiteDesign=whiteDesign, drops=drops, targetDrop=targetDrop, dropRange=dropRange,
    scale=initialRadius, parScale=rep(1, p), nfinite=sum(keep),
    rounds=1, diagnostics=diagnostics, profile=profiled)
}

ctOptimComputeUncertainty <- function(est, standata, sm, lpgFunc,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg'),
  finishsamples=1000, cores=1, matsetup=NA, control=list(), verbose=0,
  scores=NULL, hessian=NULL){
  
  uncertainty <- match.arg(uncertainty)
  ctOptimCheckUncertaintyData(standata=standata, uncertainty=uncertainty,
    finishsamples=finishsamples, npars=length(est))
  if(is.null(control$ridge)) control$ridge <- 1e-8
  if(is.null(control$hessianStep)) control$hessianStep <- 1e-3
  if(is.null(control$surrogateScale)) control$surrogateScale <- .5
  if(is.null(control$surrogateNpoints)) control$surrogateNpoints <- NULL
  if(is.null(control$surrogateProfile)) control$surrogateProfile <- TRUE
  if(is.null(control$surrogateProfileTargetDrop)) {
    control$surrogateProfileTargetDrop <- NULL
  }
  if(is.null(control$surrogateProfileMaxStep)) {
    control$surrogateProfileMaxStep <- 64
  }
  
  base <- suppressMessages(suppressWarnings(lpgFunc(est)))
  base_gradient <- attributes(base)$gradient
  if(is.null(base_gradient) || any(!is.finite(base_gradient))) {
    base_gradient <- rep(0, length(est))
  }
  
  hessian_result <- NULL
  scoremat <- NULL
  draws <- NULL
  method_details <- list()
  covavailable <- FALSE
  
  if(uncertainty == 'surrogate' && !is.null(control$initialCov)){
    cov <- ctOptimSafeCov(control$initialCov, ridge=control$ridge)
    covavailable <- TRUE
  }
  
  if(uncertainty %in% c('hessian','sandwich','bootstrap','is') ||
      (uncertainty == 'surrogate' && !covavailable)){
    if(!is.null(hessian)){
      # An exact Hessian supplied by the caller -- see .ctBackendHessian(),
      # which gets one from the julia engine by differentiating its own
      # reverse-mode gradient. There is nothing to difference, so there is no
      # step to choose and no pair of one-sided estimates to reconcile.
      hess <- hessian
    } else {
      message('Estimating Hessian')
      # Both directions on one progress line (see numericHessianFunc()).
      progress <- .ctBackendProgressSink(0)
      hess1 <- numericHessianFunc(pars=est, step=control$hessianStep,
        verbose=verbose, directions=1, lpgFunc=lpgFunc,
        base_value=base[1], base_gradient=base_gradient, progress=progress)
      hess2 <- numericHessianFunc(pars=est, step=control$hessianStep,
        verbose=verbose, directions=-1, lpgFunc=lpgFunc,
        base_value=base[1], base_gradient=base_gradient, progress=progress)
      if (!is.null(progress)) progress('', 'break')
      hessian_result <- processHessianMatrices(hess1, hess2, verbose, matsetup)
      hess <- hessian_result$hess
    }
    # Which of the flat directions in this curvature the likelihood is really
    # flat along, measured rather than inferred from the eigenvalue. Costs
    # nothing on a fit with no flat direction -- the screen returns before
    # evaluating anything -- and is what keeps the answer from depending on
    # where along such a direction the optimiser stopped. See
    # `.ctOptimFlatDirectionScreen()`.
    #
    # Only on the curvature routes. `opg`, `sandwich` and `bootstrap` build
    # their matrix from score contributions rather than from the likelihood's
    # own second derivatives, so walking the likelihood along an eigenvector of
    # *that* is not the question this answers.
    #
    # `control$flatScreen = FALSE` turns it off, the same shape as
    # `analyticHessian`: a caller who wants the eigenvalue rule on its own --
    # to reproduce an older result, or to measure what this is worth -- has
    # asked for it, and accepting the argument and ignoring it would be worse
    # than either answer.
    screen <- if(identical(control$flatScreen, FALSE)) NULL else
      .ctOptimFlatDirectionScreen(-(hess + t(hess)) / 2, lpgFunc, est)
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge, screen=screen)
    # A point better than the estimate, found while asking a different
    # question. Reported whether or not anything was confirmed flat: it says
    # the fit is not at a maximum, which is a more serious finding than
    # anything else this stage produces.
    if(!is.null(screen) && screen$gain > 0) {
      method_details$notmaximum <- list(gain = screen$gain,
        point = screen$point, direction = screen$direction)
    }
    # Kept whenever the screen asked anything, flat or not, so a reader can
    # tell "asked, and nothing was flat" from "never asked". The vectors are
    # the confirmed directions themselves -- at most `maxdirections` columns --
    # because the identifiability report names what this confirmed, and a
    # direction is only nameable by its loadings.
    if(!is.null(screen)) {
      method_details$flatdirections <- list(
        n = sum(screen$flat), bar = screen$bar, lengths = screen$lengths,
        rtol = screen$rtol, candidates = length(screen$candidates),
        evaluations = screen$evaluations,
        change = screen$change[screen$flat],
        eigenvalue = screen$eig$values[screen$flat] /
          max(screen$eig$values),
        vectors = screen$eig$vectors[, screen$flat, drop = FALSE])
    }
    covavailable <- TRUE
  }
  
  if(uncertainty == 'opg'){
    message('Estimating score / OPG covariance')
    score_hessian <- bootstrapHessian(standata=standata, sm=sm, est=est,
      finishsamples=finishsamples, cores=cores, scores=scores)
    hess <- score_hessian$hess
    scoremat <- score_hessian$scores
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge)
    method_details$opg <- list(
      note='OPG-style information estimate; local prior curvature is not represented unless it appears in score variability.'
    )
  }
  
  if(uncertainty == 'fullbootstrap'){
    message('Fitting full bootstrap refits')
    fullbootstrap <- ctOptimFullBootstrapDraws(est=est, standata=standata,
      sm=sm, n=finishsamples, cores=cores, control=control,
      verbose=verbose)
    draws <- fullbootstrap$draws
    cov <- ctOptimSafeCov(stats::cov(draws), ridge=control$ridge)
    method_details$fullbootstrap <- fullbootstrap
    method_details$fullbootstrap$draws <- NULL
  }
  
  if(uncertainty %in% c('sandwich','bootstrap')){
    message('Computing score contributions')
    scoremat <- ctOptimScoreMatrix(standata=standata, sm=sm, est=est,
      cores=cores, scores=scores)
    centered <- scale(scoremat, center=TRUE, scale=FALSE)
    meat <- crossprod(centered)
    if(uncertainty == 'sandwich'){
      cov <- cov %*% meat %*% cov
      cov <- ctOptimSafeCov(cov, ridge=control$ridge)
    } else {
      draws <- ctOptimBootstrapDraws(est=est, cov=cov, scores=scoremat,
        n=finishsamples)
      cov <- ctOptimSafeCov(stats::cov(draws), ridge=control$ridge)
    }
  }
  
  if(uncertainty == 'surrogate'){
    message('Estimating local quadratic surrogate')
    surrogate <- ctOptimSurrogateHessian(est=est, lpgFunc=lpgFunc, cov=cov,
      npoints=control$surrogateNpoints, scale=control$surrogateScale,
      ridge=control$ridge,
      profile=control$surrogateProfile,
      profileTargetDrop=control$surrogateProfileTargetDrop,
      profileMaxStep=control$surrogateProfileMaxStep,
      verbose=verbose)
    hess <- surrogate$hessian
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge)
    method_details$surrogate <- surrogate
  }
  covDiagnostics <- attr(cov, 'ctOptimCovFromHessian')
  if(!is.null(covDiagnostics)) method_details$covariance <- covDiagnostics
  
  if(!is.null(hessian)) method_details$hessian <- list(
    source='exact, by forward-mode differentiation of the reverse-mode gradient')
  list(method=uncertainty, cov=cov, hessian=if(exists('hess')) hess else NULL,
    scores=scoremat, draws=draws, base_value=base[1],
    base_gradient=base_gradient, details=method_details)
}

ctOptimUpdateTransformed <- function(fit, samples, cores=1){
  savesubjectmatrices <- fit$standata$savesubjectmatrices
  sdat <- fit$standata
  if(!as.logical(savesubjectmatrices)) sdat <- standatact_specificsubjects(sdat, 1)
  fit$stanfit$transformedpars <- stan_constrainsamples(sm=fit$stanmodel,
    standata=sdat, savesubjectmatrices=savesubjectmatrices,
    savescores=fit$standata$savescores,
    dokalman=as.logical(savesubjectmatrices), samples=samples,
    cores=cores, quiet=TRUE)
  sds <- try(suppressWarnings(sqrt(diag(fit$stanfit$cov))), silent=TRUE)
  if('try-error' %in% class(sds)) sds <- rep(NA_real_, length(fit$stanfit$rawest))
  smf <- stan_reinitsf(fit$stanmodel, fit$standata)
  fit$stanfit$transformedpars_old <- NA
  try(fit$stanfit$transformedpars_old <- cbind(
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest - 1.96 * sds)),
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest)),
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest + 1.96 * sds))),
    silent=TRUE)
  try(colnames(fit$stanfit$transformedpars_old) <- c('2.5%','mean','97.5%'),
    silent=TRUE)
  fit
}

ctOptimDataLpgFunc <- function(sm, standata, cores=1){
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  cores <- min(cores, length(unique(standata$subject)))
  
  if(cores <= 1){
    smuse <- sm
    if(!is.null(standata$recompile) && standata$recompile == 0) {
      smuse <- utils::getFromNamespace("stanmodels", "ctsem")$ctsm
    }
    smf <- stan_reinitsf(smuse, standata)
    lpg <- function(parm){
      out <- try(rstan::log_prob(smf, upars=parm, adjust_transform=TRUE,
        gradient=TRUE), silent=FALSE)
      if('try-error' %in% class(out) || is.nan(out)) {
        out <- -1e100
        attributes(out) <- list(gradient=rep(0, length(parm)))
      }
      out
    }
    return(list(lpg=lpg, cl=NULL, standata=standata, cores=1L))
  }
  
  smfile <- ''
  if(standata$recompile > 0){
    smfile <- file.path(tempdir(), paste0('ctsem_sm_',
      ceiling(stats::runif(1, 0, 100000)), '.rda'))
    save(sm, file=smfile, eval.promises=FALSE, precheck=FALSE)
  }
  cl <- makeClusterID(cores)
  parallelStanSetup(cl=cl, standata=standata, split=TRUE,
    smfile=if(standata$recompile > 0) smfile else '')
  lpg <- function(parm){
    clusterIDexport(cl, 'parm')
    out2 <- parallel::clusterEvalQ(cl=cl, parlp(parm))
    out <- try(sum(unlist(out2)), silent=TRUE)
    for(i in seq_along(out2)){
      if(i == 1) attributes(out)$gradient <- attributes(out2[[1]])$gradient
      if(i > 1) attributes(out)$gradient <-
          attributes(out)$gradient + attributes(out2[[i]])$gradient
    }
    if('try-error' %in% class(out) || is.nan(out)) {
      out <- -1e100
      attributes(out) <- list(gradient=rep(0, length(parm)))
    }
    out
  }
  list(lpg=lpg, cl=cl, standata=standata, cores=cores, smfile=smfile)
}

ctOptimFitLpgFunc <- function(fit, cores=1){
  standata <- fit$standata
  standata$savesubjectmatrices <- 0L
  ctOptimDataLpgFunc(sm=fit$stanmodel, standata=standata, cores=cores)
}

# Defined above that block, not between it and the function it documents.
# Roxygen attaches a block to whatever definition follows it, so sitting
# below it cost `ctOptimUncertainty` its export: a regenerated NAMESPACE
# carried `export(.ctResolveDraws)` in its place, un-exporting a documented
# user-facing function that three error messages tell people to call.
# How samples are produced, given how the covariance was estimated.
#
# `uncertainty` and `draws` are not independent, and pretending otherwise is
# what made this confusing: `uncertainty='is'` *is* importance sampling, so
# `draws='normal'` alongside it was accepted, warned about, and then overridden
# -- the argument appeared to be a choice and was not. It is now derived, and an
# incompatible request is an error rather than a warning about something the
# code went on to ignore.
.ctResolveDraws <- function(uncertainty, draws) {
  implied <- switch(uncertainty,
    is = 'imis',
    bootstrap = 'empirical',
    fullbootstrap = 'empirical',
    'normal')
  if(identical(draws, 'auto')) return(implied)
  if(!identical(draws, implied)) {
    stop("uncertainty='", uncertainty, "' produces draws by '", implied,
      "', so draws='", draws, "' cannot be honoured. Pass draws='auto' (the ",
      "default) or draws='", implied, "'.", call.=FALSE)
  }
  draws
}

# uncertainty='stored': redraw from the covariance the fit already carries.
#
# Every other method here builds a covariance out of model evaluations -- 2*npar
# log-probability/gradient calls for the finite-difference Hessian, more for the
# score and bootstrap methods -- and then draws from it. Asking for a different
# number of draws is not asking for that work again: the draws are iid from a
# covariance that has not changed, so this path costs no model evaluations at
# all on either backend.
#
# It replaces `ctFitAddSamples()`, which did the same thing on stan alone and
# appended rather than replaced. Appending is what makes it stan-only in spirit
# as well as in code: the rows already on the fit may have come from `is` or
# `bootstrap`, and adding normal draws to those leaves a posterior that is part
# one distribution and part another, with nothing recording the mixture.
# Replacing them cannot do that, and the warning below says when the previous
# draws were of a kind this cannot reproduce.
# Turn a computed covariance into the draws a fit reports, given which kind of
# draws were asked for.
#
# This was written twice -- once in `ctOptimUncertainty()`'s stan branch and
# once in `.ctBackendUncertainty()` -- around the same shared core
# (`ctOptimComputeUncertainty()`), in the same order, with the same five
# `imis*` constants filled in by hand on each side. The julia copy's comment
# records what that cost: `empirical` was missing there, so
# `uncertainty='bootstrap'` fell through to normal draws while recording
# `draws='empirical'`, and nothing looked wrong because the *covariance* was
# still the bootstrap's.
#
# The two genuine differences are arguments here rather than hardcoded
# constants, so that the divergence is visible in one place instead of being
# two similar-looking blocks:
#
#   `scaleInit`/`tailScale` -- stan uses 1.1/1.1, julia 1.5/1.2. Deliberate on
#   both sides and measured: a proposal narrower than its target cannot correct
#   it, which argues for the wider default, but a 400-subject model measures
#   better at 1.1 where a 40-subject one measures better at 1.5. One constant
#   does not serve both sample sizes; see the note in `.ctBackendUncertainty()`.
#
#   `lpg` -- where the log density comes from. Value-only on julia, because
#   `imis_is` reads the log probability and nothing else, so a `lpgFunc` that
#   also computes a reverse pass per draw has that work thrown away.
#
#   `tailremedy` -- how to sample the posterior itself when the weights' Pareto
#   k says reweighting cannot (see `.ctOptimImisReport()`), which is spelled
#   differently on each backend.
#
# @return list(samples, uncertaintyfit, control) -- `control` comes back
#   because the defaults filled in here are what gets recorded in `$settings`.
.ctOptimDrawSamples <- function(uncertaintyfit, draws, control, est,
  finishsamples, lpg, verbose = 0,
  scaleInit = .ctImisProposalDefaults()$scaleInit,
  tailScale = .ctImisProposalDefaults()$tailScale,
  df = .ctImisProposalDefaults()$df,
  tailremedy = 'Sample the posterior rather than reweighting an approximation to it.') {

  if (draws == 'empirical' && !is.null(uncertaintyfit$draws)) {
    return(list(samples = uncertaintyfit$draws, uncertaintyfit = uncertaintyfit,
      control = control))
  }

  if (draws != 'imis') {
    return(list(samples = ctOptimNormalDraws(est, uncertaintyfit$cov, finishsamples),
      uncertaintyfit = uncertaintyfit, control = control))
  }

  if (is.null(control$imisMaxIter)) control$imisMaxIter <- 50
  if (is.null(control$imisScaleInit)) control$imisScaleInit <- scaleInit
  if (is.null(control$imisTailScale)) control$imisTailScale <- tailScale
  # t with 5 degrees of freedom, not normal: see `.ctImisProposalDefaults()`.
  if (is.null(control$imisDf)) control$imisDf <- df
  # `.ctEssTarget`, the one default every draw-producing route shares.
  if (is.null(control$isESS)) control$isESS <- .ctEssTarget
  if (is.null(control$isitersize)) control$isitersize <- 1000

  # `.ctImisRun()` runs in the identified subspace of `uncertaintyfit$cov` when
  # it is rank deficient -- an individually varying parameter with no
  # individual differences behind it, and its raw correlations, are the usual
  # cause -- and passes straight through to `imis_is()` otherwise; see
  # `IS-importance-sampling-2026-09-06.md`. `lpg` value-only carries the batch
  # attribute `.ctBackendLpgFunc(fit, gradient=FALSE)` sets, whitened or not.
  is_res <- .ctImisRun(lpg, centre = est, cov = uncertaintyfit$cov,
    max_iter = control$imisMaxIter, scale_init = control$imisScaleInit,
    tail_scale = control$imisTailScale, df = control$imisDf,
    target_ess = control$isESS, n_batch = control$isitersize, cl = NA,
    finishsamples = finishsamples,
    # `verbose > 0`, not TRUE: this printed IMIS iteration progress at
    # `verbose = 0`, so the one argument meant two things across the backends --
    # silence on julia, a page of output on stan.
    verbose = verbose > 0)
  subspace <- attr(is_res, 'subspace')

  samples <- is_res$theta
  uncertaintyfit$proposal_cov <- uncertaintyfit$cov
  # `$details$covariance` diagnoses the covariance `ctOptimComputeUncertainty()`
  # produced, which after this point is the *proposal* rather than the reported
  # one. Renaming it was done on stan and not on julia -- so a julia fit
  # corrected by importance sampling carried diagnostics describing a matrix it
  # was no longer reporting. Unconditional here, which changes that julia field
  # and nothing else.
  if (!is.null(uncertaintyfit$details$covariance)) {
    uncertaintyfit$details$proposal_covariance <- uncertaintyfit$details$covariance
    uncertaintyfit$details$covariance <- NULL
  }
  weighted <- !is.null(is_res$covariance) && all(is.finite(is_res$covariance))
  if (weighted) {
    uncertaintyfit$cov <- ctOptimSafeCov(is_res$covariance)
  } else if (nrow(samples) > 1) {
    uncertaintyfit$cov <- ctOptimSafeCov(stats::cov(samples))
  }
  uncertaintyfit$imis <- is_res
  report <- .ctOptimImisReport(is_res, control$isESS, weighted,
    tailremedy = tailremedy)
  uncertaintyfit$details$importance_sampling <- list(ess = is_res$ess,
    # NA when loo is not installed; see `.ctImisParetoK()`.
    pareto_k = report$k,
    df_used = is_res$df_used, weighted = weighted,
    covariance = if (weighted) 'weighted importance-sampling covariance' else
      'unweighted covariance of the resampled draws',
    # NULL on the ordinary fit, where nothing was held. Positional in the raw
    # parameter vector, as `ctOptimCovFromHessian()`'s `nullParameters` is,
    # because names are not attached to this vector until
    # `.ctFitNameRawUncertainty()` runs, further down the caller.
    subspace = if (is.null(subspace)) NULL else list(
      nullDirections = subspace$nnull, heldParameters = subspace$nullParameters),
    # How much it cost and what the profile-path search found: the log
    # density's evaluations, the rounds of draws, and the parameters (or
    # subspace directions) whose paths had a tail to follow; see `imis_is()`.
    evaluations = is_res$evaluations, rounds = nrow(is_res$rounds),
    paths = is_res$paths)

  list(samples = samples, uncertaintyfit = uncertaintyfit, control = control)
}

.ctOptimStoredRedraw <- function(fit, finishsamples, cores, verbose=0){
  julia <- inherits(fit, 'ctJuliaFit')
  cov <- if(julia) fit$estimate$cov else fit$stanfit$cov
  est <- if(julia) as.numeric(fit$estimate$raw) else fit$stanfit$rawest
  if(is.null(cov) || !length(cov) || any(!is.finite(cov))) stop(
    "uncertainty='stored' redraws from the covariance already on the fit, and ",
    "this fit has no usable one. Run ctFitUncertainty() with a method that ",
    "computes one first -- 'hessian' is the default.", call.=FALSE)
  if(is.null(est) || length(est) != ncol(cov)) stop(
    "uncertainty='stored' needs the fit's raw estimate and stored covariance ",
    "to describe the same parameters; they do not.", call.=FALSE)

  uncertaintyfit <- if(julia) fit$uncertainty else fit$stanfit$uncertainty
  if(is.null(uncertaintyfit)) uncertaintyfit <- list(method='stored', cov=cov)
  previousdraws <- uncertaintyfit$settings$draws
  if(!is.null(previousdraws) && !identical(previousdraws, 'normal')) {
    warning("This fit's draws came from '", previousdraws, "'; redrawing from ",
      "the stored covariance gives normal draws instead. Rerun with ",
      "uncertainty='", uncertaintyfit$settings$method, "' to keep them.",
      call.=FALSE)
  }

  samples <- ctOptimNormalDraws(est, cov, finishsamples)

  # `method` keeps naming the method that produced the covariance, because that
  # is what every reader of it wants to know and it has not changed; `redrawn`
  # records that the draws were regenerated from it afterwards.
  uncertaintyfit$draws <- 'normal'
  uncertaintyfit$settings$method <- .ctJuliaOr(uncertaintyfit$settings$method,
    'stored')
  uncertaintyfit$settings$draws <- 'normal'
  uncertaintyfit$settings$finishsamples <- finishsamples
  uncertaintyfit$settings$cores <- cores
  uncertaintyfit$settings$redrawn <- TRUE

  if(julia){
    fit$estimate$rawposterior <- samples
    fit <- .ctFitNameRawUncertainty(fit)
    fit$uncertainty <- uncertaintyfit
    fit$transformedpars <- .ctBackendConstrain(fit)
    return(fit)
  }
  fit$stanfit$rawposterior <- samples
  fit <- .ctFitNameRawUncertainty(fit)
  fit$stanfit$uncertainty <- uncertaintyfit
  if(verbose > 0) message('Redrawing ', finishsamples,
    ' samples from the stored covariance')
  ctOptimUpdateTransformed(fit, samples=samples, cores=cores)
}

#' Compute or sample a fit's uncertainty
#'
#' Recomputes the approximate raw-parameter uncertainty for an optimized
#' \code{\link{ctFit}} object and refreshes the approximate raw-parameter
#' samples. This is the entry point for both backends; \code{ctFit} itself
#' calls it to finish an optimized fit. \code{\link{ctOptimUncertainty}} is
#' the previous name, kept for fits written against ctsem 3.11.1; the two are
#' otherwise identical except that only this one offers
#' \code{uncertainty = 'sample'}.
#'
#' Every method except \code{'sample'} writes \emph{pseudo-posterior} draws to
#' \code{$rawposterior}: a sample from a covariance fitted to the
#' log-posterior surface around the optimum, not a sample from the posterior
#' itself, and the point estimate stays at the optimum throughout.
#' \code{uncertainty = 'sample'} is the other thing, genuine posterior draws
#' by Markov chain Monte Carlo, and it is the one method that moves the point
#' estimate: \code{$estimate$raw} (or \code{$stanfit$rawest}) becomes the
#' posterior mean, the way \code{ctFit(optimize = FALSE)} already reports a
#' sampled fit. Everything downstream reads either kind of draw from the same
#' slot, so the difference is recorded rather than visible in the shape of the
#' result: a curvature-based fit carries \code{$uncertainty$settings}, a
#' sampled one carries \code{$sample}, and every method except
#' \code{'sample'} itself refuses a fit that already carries a posterior
#' rather than replacing its draws with an approximation -- see
#' \code{uncertainty = 'sample'} below for what re-sampling one does instead.
#'
#' To change only the number of draws, use \code{uncertainty='stored'}, which
#' redraws from the covariance the fit already carries and evaluates no model.
#'
#' @section Backend differences:
#' The methods are the same on both backends and the covariance they produce is
#' comparable, but three things differ and are worth knowing before comparing
#' output.
#'
#' \emph{Where the result is stored.} A \code{ctStanFit} comes back with
#' \code{fit$stanfit$cov}, \code{fit$stanfit$rawposterior} and
#' \code{fit$stanfit$uncertainty}; a \code{ctJuliaFit} with
#' \code{fit$estimate$cov}, \code{fit$estimate$se},
#' \code{fit$estimate$rawposterior} and \code{fit$uncertainty}. Both record the
#' resolved settings in \code{$uncertainty$settings}. On both backends the
#' draws, the covariance and the standard errors are labelled by raw parameter,
#' so \code{fit$estimate$se['drift']} or
#' \code{apply(fit$estimate$rawposterior, 2, quantile)} reads without counting
#' columns. The two backends spell the population-SD and correlation blocks
#' differently (\code{popsd_x}/\code{rawcor_x__y} on julia,
#' \code{x_SD}/\code{y_x_corr} on stan); the order is the same.
#'
#' \emph{The Hessian.} The stan path finite-differences its gradient with a
#' single global step (\code{control$hessianStep}). The julia engine
#' differentiates its own reverse-mode gradient in forward mode, so its
#' Hessian is exact; \code{control$analyticHessian = FALSE} falls back to the
#' shared finite difference for comparison. The two agree to about 1e-4 on a
#' well-conditioned model.
#'
#' \emph{Directions with no curvature.} A direction whose curvature is
#' negligible against the sharpest one is left out of the inversion rather
#' than floored at \code{control$ridge}, so the parameters that carry it come
#' back with no spread at all instead of a standard error of \code{1/ridge}.
#' Floored, that number leaks into the parameters the data \emph{does}
#' determine, and by an amount rounding decides: it can differ by a factor of a
#' hundred between two fits that reach the same optimum.
#' A direction with little curvature along which the likelihood itself is
#' measured flat -- within the likelihood-ratio bar of 1.92 over four raw
#' units on one side, following the ridge where it curves -- is treated the
#' same way whether or not its curvature has yet decayed to negligible, since
#' that depends on where along such a ridge the optimiser stopped.
#' \code{$uncertainty$cov} records which directions were dropped, and
#' \code{fit$identifiability} names the parameters.
#'
#' \emph{Whether the intervals are as wide as the curvature allows.} A julia
#' fit carries \code{fit$uncertainty$intervalcheck}: for each raw parameter,
#' the reported standard error against \code{1/sqrt(information[i,i])}, the
#' width that parameter's own curvature supports. Their ratio is 1 when the
#' parameter is separable from the rest and grows without bound as it stops
#' being; anything past about 100 means the reported width comes from the
#' entanglement rather than from the data, and will not repeat between runs.
#' The same object answers the opposite question, which the ratio cannot:
#' \code{$nullmass} is each coordinate's share of the directions left out of
#' the inversion, and \code{$unidentified} names those with enough of it that
#' their true asymptotic variance is infinite. Those get almost none of the
#' variance the projected inverse has to give, so they would otherwise be
#' reported as precisely estimated -- more precisely the more data there is.
#' \code{summary()} reports their sd, interval and z as \code{NA}.
#'
#' \emph{Backend-specific arguments.} \code{uncertainty='fullbootstrap'} and
#' its \code{control$bootstrapFitCores} / \code{control$bootstrapTol}, and
#' \code{control$parsteps}, are stan-only and are refused by name on a
#' \code{ctJuliaFit}. \code{control$analyticHessian} is julia-only.
#' \code{cores} means R worker processes on stan and engine threads on julia
#' (see below). The IMIS proposal defaults also differ:
#' \code{imisScaleInit = 1.1} and \code{imisTailScale = 1.1} on stan against
#' \code{1.5} and \code{1.2} on julia. Both are measured, on different
#' regimes: the narrower pair costs fewer evaluations to reach the target
#' effective sample size on a well-identified fit, and the wider pair is what
#' a small-sample posterior genuinely wider than the Hessian curvature needs
#' to be seen at all -- a proposal no wider than the curvature cannot correct
#' a posterior wider than it -- which is the setting the julia default was
#' raised for. Set them explicitly to compare the two. On either backend, a
#' proposal covariance with a raw direction the data does not identify -- an
#' individually varying parameter with no individual differences behind it is
#' the usual cause -- is sampled in the identified subspace only, holding that
#' direction at the estimate rather than manufacturing an importance weight
#' for a density that is not one; see \code{uncertainty='is'} below.
#'
#' Transformed-parameter summaries are refreshed on stan and not on julia,
#' which has no parameter-matrix reconstruction through
#' \code{rstan::constrain_pars}; the raw-scale covariance and draws are
#' complete on both.
#'
#' @param fit Optimized \code{ctStanFit} or \code{ctJuliaFit} object. For a
#' \code{ctJuliaFit}, every \code{uncertainty} method except
#' \code{'fullbootstrap'} is available; that one re-optimises each resample and
#' so needs the model rebuilt rather than re-evaluated. \code{uncertainty =
#' 'sample'} needs \code{backend = 'julia'} specifically -- it is refused by
#' name on a \code{ctStanFit}, which samples through \code{ctFit(backend =
#' 'stan', optimize = FALSE)} instead. A sampled fit of either backend is
#' refused by every \emph{other} method: it already carries a posterior, and
#' replacing it with a curvature-based approximation would discard it.
#' \code{uncertainty = 'sample'} may be run again on an already-sampled fit,
#' to draw more, or differently, from where it now stands.
#' @param uncertainty Uncertainty approximation. \code{'hessian'} uses the
#' Hessian at the estimate -- exact on julia, which differentiates its own
#' gradient, and by finite differences on stan -- \code{'surrogate'} fits a local quadratic
#' surrogate around the optimum, \code{'is'} runs adaptive importance sampling
#' against the fitted log posterior. Each raw parameter is first walked out on
#' both sides of the mode, from 3 to 96 standard errors of the curvature, with
#' the other parameters moved to their conditional mode, and a proposal
#' component is placed wherever that walk finds a tail the curvature does not
#' show -- a skewed one, or one bending away along a ridge the other
#' parameters follow, as a variance with few subjects often has. Rounds of
#' \code{isitersize} draws from these multivariate t components
#' (\code{imisDf}) then add components at the highest-weighted draws until the
#' effective sample size reaches \code{isESS} and the Pareto k of the weights
#' is below 0.7 (k needs the loo package), each draw weighted against the whole
#' mixture. The effective size, k, the evaluations used and which parameters
#' had a tail to follow are recorded in
#' \code{fit$uncertainty$details$importance_sampling}; k above 0.7 is warned
#' of, and so is a parameter whose walk has not fallen off 96 standard errors
#' out, which is what an improper posterior looks like. When the Hessian
#' covariance is rank deficient the sampling runs in the whitened
#' eigen-coordinates of the directions it has curvature in, holding the rest
#' at the estimate. It needs the log posterior's gradient for the walk, which
#' both backends' fits supply.
#'
#' The density \code{'is'} weights toward is the fit's own objective. On a
#' fit whose random effects are integrated -- \code{intoverpop = 'laplace'}
#' or \code{'augmented'} -- that is the approximate marginal posterior the route
#' maximised, so the draws correct the normal approximation's shape but keep
#' the route's own approximation: where the Laplace approximation is biased
#' (random effects on variances, in nonlinear dynamics, or with non-Gaussian
#' indicators), \code{'is'} reproduces the bias faithfully. Measured against
#' long reference samples on five bench models and two closed forms where
#' that approximation is accurate, its worst errors in a posterior sd and in
#' the 2.5\% and 97.5\% quantiles were at most those of \code{'sample'} at its
#' defaults, in a twentieth to a half of the time; on twenty models including
#' those where it is not, it cost 5 to 12 times the fit and was no more
#' accurate than the Hessian on balance. \code{'sample'} draws from the exact
#' posterior.
#' \code{'bootstrap'} uses one-step score bootstrap draws with
#' Hessian bread, \code{'fullbootstrap'} resamples subjects and fully
#' re-optimizes each sample from the original maximum likelihood or MAP
#' estimate using mize L-BFGS, \code{'sandwich'} uses Hessian bread with score
#' covariance meat, and \code{'opg'} uses an OPG-style score information
#' approximation. \code{'stored'} computes nothing: it redraws
#' \code{finishsamples} normal draws from the covariance already on the fit and
#' leaves that covariance, the Hessian and the recorded \code{method}
#' untouched, marking \code{$uncertainty$settings$redrawn = TRUE}. It is the
#' way to change the number of draws, or to reseed them, without paying for the
#' covariance again -- the other methods cost at least \code{2 * npar}
#' log-probability evaluations, this one costs none -- and it warns if the
#' fit's existing draws came from \code{'is'} or \code{'bootstrap'}, which
#' normal draws from that covariance do not reproduce.
#' \code{'sample'} draws from the genuine posterior by Markov chain Monte
#' Carlo -- SAEM's kernel on the joint posterior of parameters and random
#' effects, the No-U-Turn sampler otherwise -- through the same runner
#' \code{ctFit(backend = 'julia', optimize = FALSE)} uses to fit and sample
#' together -- \code{julia}-only, see \code{fit} above. Its settings are
#' entries of \code{control} rather than \code{draws}/\code{finishsamples},
#' which this method ignores: \code{chains} (default 4), \code{warmup} (500),
#' \code{draws} (500), \code{seed}, \code{saveEffects} (FALSE, whether to keep
#' every draw of every random effect rather than only their summary; either way
#' they are each parameter's deviation from its population value on the raw
#' scale, per subject or group, in \code{fit$sample$effects},
#' \code{effect_mean} and \code{effect_sd}), and
#' \code{processes} (TRUE, one R process per chain). A run stops once every
#' parameter's effective sample size reaches \code{minESS} (200) and the
#' chains' R-hat is below \code{rhatTarget} (1.01) -- \code{meanESS} adds a
#' target for the average -- within a budget of four times \code{draws} per
#' chain, or \code{maxDraws} when given; \code{minESS = 0} takes exactly
#' \code{draws}. A run that reaches its budget short of the target warns.
#'
#' Lower targets give a quicker, approximate posterior, and say how
#' approximate. A reported quantile's Monte Carlo error is about
#' \code{1.96 * 2.1 / sqrt(ESS)} posterior sds at the 5\% and 95\% points (the
#' tail ESS) and \code{1.96 * 1.25 / sqrt(ESS)} at the median, and chains
#' whose means disagree by \code{d} posterior sds give an R-hat of about
#' \code{sqrt(1 + d^2)}; measured on short runs, the error these imply was
#' about a third optimistic. So \code{minESS = 50, rhatTarget = 1.07} aims at
#' about half a posterior sd at the centre and more in the tails, and the
#' defaults at about a fifth at the centre and two fifths in the tails. On the models compared, runs of four chains
#' and 100 to 200 draws each overtook the Laplace normal approximation's
#' accuracy at two to three times a Laplace fit's time, and were faster and
#' more accurate outright where the Laplace fit was slow or biased; short runs
#' are limited by whether the chains agree, so R-hat is usually the binding
#' check. With \code{processes = TRUE} the chains run in up to \code{cores}
#' worker processes, each starting its own Julia session -- a fixed cost of
#' some tens of seconds that a run much shorter than a minute does not recover
#' -- and a worker runs its share of the chains in turn, so four chains at the
#' default of two cores run two at a time. \code{processes = FALSE}, or one
#' core, keeps the chains in this session, as threads when it has them.
#' \code{settleTol} ends warmup early once the metric stops moving (off by
#' default, and slower when measured). \code{control$sampler} chooses the
#' kernel for the joint posterior: \code{'saem'} (the default there) draws the
#' random effects by SAEM's Metropolis-within-Gibbs sweeps and the parameters
#' given them by NUTS, with moves that re-express the effects for each level's
#' means and scales; \code{'nuts'} runs NUTS on the parameters and every
#' random effect at once. Both draw the same posterior, from the same
#' placement, and stop by the same rule; SAEM's kernel was the more reliable
#' and the faster of the two in comparisons on variance and nested models,
#' where NUTS on the joint vector mixes poorly. \code{adapt_metric},
#' \code{adapt_effects} and \code{settleTol} belong to NUTS and are refused
#' with SAEM's kernel. A marginal target is sampled by NUTS. \code{control$placement} says where the chains
#' start: \code{'saem'}, the default on the joint posterior, runs SAEM from the
#' fit's estimate -- on the exact marginal posterior, where the fit's optimum is
#' the Laplace approximation's -- and starts each chain from SAEM's estimate and
#' its draws of the random effects, for either sampler; \code{'fit'} starts them
#' around the fit's own estimate, and is the only placement for a marginal
#' target, which has no random effects to start. \code{control$target}
#' says which posterior: \code{'auto'} (the default) is the one the fit's
#' \code{intoverpop} names -- the filter's marginal for \code{'augmented'}, the
#' Laplace marginal for a fit that asked for \code{'laplace'} (an approximate
#' posterior whose dimension does not grow with the subject count), and the
#' joint posterior over population parameters \emph{and} every subject's
#' random effects for \code{'none'} or for a Laplace route that
#' \code{intoverpop = 'auto'} chose, where the Laplace fit only places the
#' chains. \code{'marginal'}/\code{'joint'} ask for one explicitly regardless
#' of route; \code{'joint'} is refused by name on an \code{intoverpop =
#' 'augmented'} fit, which has no separate random effect to sample jointly
#' with the parameters. See
#' \code{\link{ctJuliaSetup}} for the thread count that decides whether
#' chains run concurrently, and \code{\link{ctFit}}'s \code{intoverpop} for
#' what each route means.
#'
#' Every method that produces draws aims for the same effective sample size
#' unless told otherwise, 200: \code{minESS} for \code{'sample'} (of the worst
#' parameter, since the 2.5\% and 97.5\% quantiles need it for each),
#' \code{isESS} for \code{'is'}, and \code{target_ess} for
#' \code{\link{ctLaplaceCorrect}} and \code{\link{ctParticleCorrect}} with
#' \code{draws = 'imis'}. Each warns when it ends short of its target.
#' @param draws Approximate raw-parameter draw method. \code{'auto'} uses
#' empirical draws for \code{uncertainty='bootstrap'} and
#' \code{uncertainty='fullbootstrap'} and normal draws otherwise.
#' \code{'normal'} draws from a multivariate normal using the selected
#' covariance, \code{'empirical'} uses empirical draws when available, and
#' \code{'imis'} runs the importance sampler using the selected covariance as
#' proposal. For \code{uncertainty='is'}, \code{draws} is set to \code{'imis'}.
#' @param finishsamples Number of approximate raw-parameter samples. If
#' \code{NULL}, the existing number of rows in the fit's raw posterior
#' (\code{fit$stanfit$rawposterior} for stan, \code{fit$estimate$rawposterior}
#' for julia) is reused when available; otherwise 1000 samples are used.
#' @param cores The most CPU cores the call uses at once, counting every
#' process it starts. If \code{NULL}, one core is used, and nothing
#' is parallelised unless a value above one is asked for -- except for
#' \code{'sample'}, which takes \code{getOption("mc.cores", 2)} as
#' \code{ctFit} does, and shares it between at most that many worker
#' processes, which run the chains. On a
#' \code{ctStanFit} these are R worker processes: each
#' log-probability/gradient evaluation is split across subjects and reassembled,
#' and score contributions and transformed quantities use them too. On a
#' \code{ctJuliaFit} there is no R cluster; the value becomes the engine's
#' subject-loop chunk ceiling for the duration of the call, and the engine caps
#' it at its own thread count. Neither route changes what is computed, though
#' both change the order things are summed in, so results are reproducible at
#' \code{cores = 1} and agree to rounding above it.
#' @param control List of method-specific options. For \code{uncertainty =
#' 'sample'} these are the sampler settings described under \code{uncertainty}
#' above (\code{chains}, \code{warmup}, \code{draws}, \code{seed},
#' \code{saveEffects}, \code{processes}, \code{target}, \code{sampler},
#' \code{placement}, \code{minESS}, \code{meanESS}, \code{rhatTarget},
#' \code{maxDraws}, \code{settleTol}, \code{init_scale}, and NUTS's
#' \code{maxdepth}, \code{target_accept}, \code{adapt_metric} and
#' \code{adapt_effects}); none of the entries below apply to it. For every other method, useful entries include
#' \code{ridge}, \code{hessianStep}, \code{surrogateNpoints},
#' \code{surrogateScale}, \code{surrogateProfile},
#' \code{surrogateProfileTargetDrop}, \code{surrogateProfileMaxStep},
#' \code{bootstrapFitCores}, \code{bootstrapTol}, \code{imisMaxIter},
#' \code{imisScaleInit}, \code{imisTailScale}, \code{imisDf}, \code{isESS},
#' and \code{isitersize}. Omitted entries use
#' \code{ridge = 1e-8}, \code{hessianStep = 1e-3},
#' \code{surrogateScale = .5}, \code{surrogateNpoints = NULL},
#' \code{surrogateProfile = TRUE},
#' \code{surrogateProfileTargetDrop = NULL},
#' \code{surrogateProfileMaxStep = 64},
#' \code{bootstrapFitCores = 1}, \code{bootstrapTol = 1e-5},
#' \code{imisMaxIter = 50}, \code{imisScaleInit = 1.1},
#' \code{imisTailScale = 1.1} (1.5 and 1.2 on julia), \code{imisDf = 5} (the
#' proposal components' t degrees of freedom; \code{Inf} for normal),
#' \code{isESS = 200}, and \code{isitersize = 1000}. When
#' \code{surrogateNpoints} is \code{NULL}, the
#' surrogate uses at least \code{max(4 * npars, 50)} local directions. The
#' surrogate is fit in whitened coordinates relative to the proposal covariance.
#' Each direction is radially adjusted with a small evaluation budget so that
#' the retained points are close to an informative local log-probability drop,
#' rather than relying on a random cloud to land in the desired range. With
#' \code{surrogateProfile = TRUE}, all fitted surrogate curvature directions
#' are then profiled until they reach \code{surrogateProfileTargetDrop}, or the
#' surrogate target drop when \code{NULL}. The profile search uses the fitted
#' surrogate curvature to choose direction-specific starting distances, so very
#' flat directions can be checked beyond \code{surrogateProfileMaxStep} when
#' the surrogate itself predicts that a larger distance is needed. If the
#' observed profile is still flatter than the fitted surrogate predicted, a
#' small magnitude-adjusted expansion budget is used before reporting a missed
#' target. \code{parsteps} may be supplied internally to keep stepwise-fixed
#' raw parameters fixed while estimating uncertainty for the remaining
#' parameters; existing fixed indices from a previous uncertainty calculation
#' are retained when no new \code{parsteps} are supplied. It is stan-only,
#' since only \code{\link{stanoptimis}} has a stepwise phase to inherit fixed
#' parameters from, and is refused rather than ignored on a
#' \code{ctJuliaFit}. \code{analyticHessian = FALSE} is julia-only and asks for
#' the shared finite-difference Hessian instead of the engine's exact one.
#' Hessian-based covariance construction first attempts the unmodified
#' \code{solve(-hessian)} covariance and a Cholesky check. It warns when
#' numerical repair is needed, such as positive-definite projection, ridge
#' flooring of information eigenvalues, or fallback to \code{MASS::ginv};
#' diagnostics are stored in
#' \code{fit$stanfit$uncertainty$details$covariance}.
#' Score-based methods use subject-level score contributions when there are
#' at least two subjects; single-subject models warn and use case-level
#' contributions. Score-based methods warn when there are fewer than ten
#' independent subjects or no more score rows than raw parameters. Full
#' bootstrap requires at least two subjects and warns below ten independent
#' subjects. Bootstrap-style methods require at least two returned samples /
#' refits.
#' @param verbose Integer controlling progress detail.
#' @param ... Not used. Anything passed here is an error rather than ignored:
#' a method's settings are entries of \code{control}.
#'
#' @return The fit, of the class it came in as. For every method except
#' \code{'sample'}: the resolved method, draw strategy, sample count, cores,
#' and non-internal controls are recorded in
#' \code{fit$stanfit$uncertainty$settings} for a \code{ctStanFit} and in
#' \code{fit$uncertainty$settings} for a \code{ctJuliaFit}; see the backend
#' differences above for the other slots each writes. For
#' \code{uncertainty = 'sample'}: a \code{ctJuliaFit} with
#' \code{estimate$rawposterior} holding the draws and \code{$sample} holding
#' the chain diagnostics (split R-hat and effective sample size per
#' parameter, divergences, tree depths, step sizes, E-BFMI (NUTS only),
#' \code{sampler} naming the kernel, \code{converged}
#' and \code{diagnosis}, and \code{target} naming which posterior was
#' sampled), exactly as \code{ctFit(optimize = FALSE)} returns.
#' @seealso \code{\link{ctFitAddSamples}} is the deprecated stan-only
#' predecessor of \code{uncertainty='stored'}. \code{\link{ctFit}} for
#' \code{optimize = FALSE}, which reaches \code{uncertainty = 'sample'}
#' through the same pipeline that placed the fit being sampled.
#' @export
ctFitUncertainty <- function(fit,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg','stored','sample'),
  draws=c('auto','normal','empirical','imis'), finishsamples=NULL,
  cores=NULL, control=list(), verbose=0, ...){
  
  # `...` used to swallow whatever reached it: `chains = 4` here, meant for
  # `control`, ran the method with its own defaults and said nothing.
  if(...length()) {
    dotnames <- ...names()
    if(is.null(dotnames)) dotnames <- rep('', ...length())
    dotnames <- ifelse(nzchar(dotnames), dotnames, '(unnamed)')
    stop("ctFitUncertainty() does not take ", paste(dotnames, collapse=', '),
      ". A method's settings are entries of control, as in ",
      "control = list(", dotnames[1L], " = ...).", call.=FALSE)
  }
  uncertainty <- match.arg(uncertainty)

  # `'sample'` is a different kind of thing from the other eight methods --
  # genuine posterior draws by Hamiltonian Monte Carlo, through the runner
  # `ctFit(backend = 'julia', optimize = FALSE)` uses to fit and sample
  # together -- so it is dispatched here, before any of the curvature-based
  # machinery below, rather than threaded through it. Julia-only, refused by
  # name rather than left to fail inside `ctOptimComputeUncertainty()`, which
  # has no julia branch at all.
  if(identical(uncertainty, 'sample')) {
    if(!.ctFitIsJulia(fit)) {
      stop("uncertainty='sample' draws by Hamiltonian Monte Carlo through ",
        "the julia sampler; it is not available for a ctStanFit. Refit with ",
        "backend='julia', or sample a stan fit with ctFit(backend='stan', ",
        "optimize=FALSE).", call.=FALSE)
    }
    # As ctFit(optimize = FALSE): its chains run at most `cores` at a time.
    if(is.null(cores)) cores <- getOption("mc.cores", 2L)
    cores <- max(1L, suppressWarnings(as.integer(cores[1])))
    if(is.na(cores)) cores <- 1L
    # The fit's own prior is what is sampled under. An ML fit has none on its
    # population sds (or none at all), and where the data allow one near zero
    # the posterior is improper there -- said before the run, not after it.
    scope <- fit$args$input$priorscope
    varying <- isTRUE(tryCatch(.ctAnyVarying(.ctFitBaseModel(fit)),
      error = function(e) FALSE))
    if(identical(scope, 'none') || (identical(scope, 'randomCorr') && varying))
      warning("This fit has no prior on ", if(identical(scope, 'none'))
        "its parameters" else "its population sds", ", so the sampled ",
        "posterior is improper wherever the data leave one flat. ",
        "Refit with priors = TRUE for a prior on every parameter.", call. = FALSE)
    return(.ctBackendUncertaintySample(fit, control=control, cores=cores,
      verbose=verbose))
  }

  draws <- match.arg(draws)
  # backend='julia' fits reach the same ctOptimComputeUncertainty() below,
  # through a log-probability/gradient function built from their own engine;
  # see R/ctBackendUncertainty.R.
  draws <- .ctResolveDraws(uncertainty, draws)
  if(inherits(fit, 'ctJuliaFit')) {
    # `fit$sample` (class "ctSampleDiagnostics") is set only by
    # `.ctBackendSampleAssemble()`, the routine shared by `ctFit(optimize =
    # FALSE)` and `uncertainty = 'sample'` above -- so it marks a julia fit
    # built from real draws either way. Nothing below knows that:
    # `.ctBackendUncertainty()`
    # treats `fit$estimate$raw` as a point estimate, builds a Hessian or score
    # matrix around it, and overwrites `fit$estimate$rawposterior` with fresh
    # curvature-based draws. Run on a sampled fit that would silently discard
    # the actual posterior draws in favour of a Gaussian approximation around
    # their mean -- wrong in a way nothing downstream would notice, since the
    # replacement is the same shape and a plausible size. The stan branch
    # below already refuses the equivalent case (`fit$stanfit$stanfit@sim`
    # populated); this mirrors it for julia.
    if(!is.null(fit$sample)) {
      stop("ctFitUncertainty() applies to an optimized ctJuliaFit; this fit ",
        "was sampled (ctFit(optimize = FALSE) or uncertainty = 'sample') and ",
        "already carries its posterior in fit$estimate$rawposterior. Read ",
        "that directly, resample it with uncertainty = 'sample', or refit ",
        "with optimize = TRUE if a curvature-based approximation is what you ",
        "want.", call.=FALSE)
    }
    # `control$parsteps` is a `stanoptimis()` concept and nothing else: the
    # stan optimiser can hold a block of raw parameters at zero for an early
    # step, and the branch below then estimates uncertainty for the remainder
    # and pads the fixed entries back in. The julia optimiser has no such
    # phase, so there is nothing for this to name -- and it reached
    # `ctOptimComputeUncertainty()`, which never reads it, as an inert list
    # element. Same standard errors as without it, and the request recorded in
    # `$uncertainty$settings$control` as though it had been honoured. Refused
    # rather than translated, because there is no julia-side meaning to
    # translate it to.
    if(!is.null(control$parsteps)) {
      stop("control$parsteps is only available for backend='stan' fits: it ",
        "holds parameters that stanoptimis() fixed during a stepwise ",
        "optimisation, and the julia optimiser has no such step. Drop it, or ",
        "refit with backend='stan'.", call.=FALSE)
    }
    # Same rule as the stan branch below, reading the julia fit's own slot.
    # Hardcoding 1000 here meant `ctOptimUncertainty(fit)` after
    # `ctOptimUncertainty(fit, finishsamples=200)` silently resampled to 1000
    # on julia and stayed at 200 on stan, against one documented default.
    if(is.null(finishsamples)) {
      finishsamples <- if(!is.null(fit$estimate$rawposterior))
        nrow(fit$estimate$rawposterior) else 1000
    }
    if(is.null(cores)) cores <- 1L
    cores <- max(1L, suppressWarnings(as.integer(cores[1])))
    if(is.na(cores)) cores <- 1L
    if(uncertainty == 'stored') return(.ctOptimStoredRedraw(fit,
      finishsamples=finishsamples, cores=cores, verbose=verbose))
    return(.ctBackendUncertainty(fit=fit, uncertainty=uncertainty, draws=draws,
      finishsamples=finishsamples, cores=cores, control=control,
      verbose=verbose))
  }
  # Named rather than asserted. "fit must be a ctStanFit object" told a caller
  # holding some other object about a class; what it needs to say is which
  # objects this does work on, since it works on both backends' fits.
  if(!'ctStanFit' %in% class(fit)) stop(
    'ctFitUncertainty() takes an optimized ctsem fit: a ctStanFit from ',
    "ctFit(..., backend='stan') or a ctJuliaFit from ctFit(..., ",
    "backend='julia'). This object is neither.", call.=FALSE)
  if(length(fit$stanfit$stanfit@sim) > 0) {
    stop('ctFitUncertainty currently applies to optimized ctStanFit objects')
  }
  # The mirror of the julia branch's `parsteps` refusal above. `analyticHessian`
  # selects the julia engine's exact Hessian, differentiated out of its own
  # gradient; the stan path has no such thing. It was neither read nor stripped
  # here, so it reached `$uncertainty$settings$control` unchanged, and a caller
  # who passed `analyticHessian=FALSE` got a settings record indistinguishable
  # from one where the request had been honoured. Refused by name rather than
  # dropped quietly, because there is no stan-side meaning to translate it to.
  if(!is.null(control$analyticHessian)) {
    stop("control$analyticHessian is only available for backend='julia' fits: ",
      "it selects the engine's exact Hessian, which the stan path does not ",
      "have. Drop it, or refit with backend='julia'.", call.=FALSE)
  }
  if(is.null(finishsamples)) {
    finishsamples <- if(!is.null(fit$stanfit$rawposterior))
      nrow(fit$stanfit$rawposterior) else 1000
  }
  if(is.null(cores)) cores <- 1
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  # Before `ctOptimFitLpgFunc()`, which reinitialises the stan model object,
  # and before any of it is needed: a stored redraw evaluates no model.
  if(uncertainty == 'stored') return(.ctOptimStoredRedraw(fit,
    finishsamples=finishsamples, cores=cores, verbose=verbose))
  lpg_cores <- if(uncertainty %in% c('opg','fullbootstrap') &&
      draws != 'imis') 1L else cores
  lpgsetup <- ctOptimFitLpgFunc(fit, cores=lpg_cores)
  on.exit({
    if(!is.null(lpgsetup$cl)) try(parallel::stopCluster(lpgsetup$cl), silent=TRUE)
    if(!is.null(lpgsetup$smfile) && nzchar(lpgsetup$smfile)) {
      try(file.remove(lpgsetup$smfile), silent=TRUE)
    }
  }, add=TRUE)
  if(uncertainty == 'surrogate' && is.null(control$initialCov) &&
      !is.null(fit$stanfit$cov)) {
    control$initialCov <- fit$stanfit$cov
  }
  if(is.null(control$parsteps) &&
      !is.null(fit$stanfit$uncertainty$fixedpars)) {
    control$parsteps <- fit$stanfit$uncertainty$fixedpars
  }
  parsteps <- integer()
  if(!is.null(control$parsteps)) {
    parsteps <- sort(unique(as.integer(unlist(control$parsteps))))
    parsteps <- parsteps[is.finite(parsteps) & parsteps >= 1 &
        parsteps <= length(fit$stanfit$rawest)]
    control$parsteps <- NULL
  }
  freepars <- setdiff(seq_along(fit$stanfit$rawest), parsteps)
  if(length(freepars) < 1) stop('No free parameters for uncertainty estimation')
  estuse <- fit$stanfit$rawest
  lpguse <- lpgsetup$lpg
  if(length(parsteps) > 0) {
    estuse <- fit$stanfit$rawest[freepars]
    if(!is.null(control$initialCov)) {
      control$initialCov <- control$initialCov[freepars, freepars, drop=FALSE]
    }
    lpguse <- function(parm){
      fullparm <- fit$stanfit$rawest
      fullparm[freepars] <- parm
      out <- lpgsetup$lpg(fullparm)
      grad <- attributes(out)$gradient
      if(!is.null(grad) && length(grad) == length(fullparm)) {
        attributes(out)$gradient <- grad[freepars]
      }
      out
    }
  }
  matsetup <- if(!is.null(fit$setup$matsetup)) fit$setup$matsetup else NA
  uncertaintyfit <- ctOptimComputeUncertainty(est=estuse,
    standata=lpgsetup$standata, sm=fit$stanmodel, lpgFunc=lpguse,
    uncertainty=uncertainty, finishsamples=finishsamples, cores=cores,
    matsetup=matsetup, control=control, verbose=verbose)
  if(length(parsteps) > 0) {
    freecov <- uncertaintyfit$cov
    fullcov <- diag(1e-10, length(fit$stanfit$rawest))
    fullcov[freepars, freepars] <- freecov
    uncertaintyfit$free_cov <- freecov
    uncertaintyfit$freepars <- freepars
    uncertaintyfit$fixedpars <- parsteps
    uncertaintyfit$cov <- fullcov
    if(!is.null(uncertaintyfit$hessian)) {
      fullhess <- matrix(0, length(fit$stanfit$rawest),
        length(fit$stanfit$rawest))
      fullhess[freepars, freepars] <- uncertaintyfit$hessian
      uncertaintyfit$free_hessian <- uncertaintyfit$hessian
      uncertaintyfit$hessian <- fullhess
    }
    if(!is.null(uncertaintyfit$draws)) {
      fulldraws <- matrix(rep(fit$stanfit$rawest,
          each=nrow(uncertaintyfit$draws)),
        nrow=nrow(uncertaintyfit$draws), byrow=FALSE)
      fulldraws[, freepars] <- uncertaintyfit$draws
      uncertaintyfit$free_draws <- uncertaintyfit$draws
      uncertaintyfit$draws <- fulldraws
    }
  }
  
  # `scaleInit`/`tailScale` at stan's own value rather than julia's -- see
  # `.ctImisProposalDefaults()` for why one constant does not serve both.
  stanImis <- .ctImisProposalDefaults('stan')
  drawn <- .ctOptimDrawSamples(uncertaintyfit, draws = draws, control = control,
    est = fit$stanfit$rawest, finishsamples = finishsamples,
    lpg = .ctImisPointGradbatch(lpgsetup$lpg), verbose = verbose,
    scaleInit = stanImis$scaleInit, tailScale = stanImis$tailScale,
    tailremedy = paste0("ctFit(..., optimize = FALSE) samples the posterior ",
      "itself, as does ctFitUncertainty(fit, 'sample') on a backend = 'julia' fit."))
  samples <- drawn$samples
  uncertaintyfit <- drawn$uncertaintyfit
  control <- drawn$control
  
  fit$stanfit$cov <- uncertaintyfit$cov
  fit$stanfit$rawposterior <- samples
  # As on julia. A no-op when this is called from inside `stanoptimis()`, where
  # the fit is still a stub with no model attached to read names from; `ctFit()`
  # calls it again on the assembled object.
  fit <- .ctFitNameRawUncertainty(fit)
  fit$stanfit$uncertainty <- uncertaintyfit
  fit$stanfit$uncertainty$draws <- draws
  storedControl <- control
  storedControl$initialCov <- NULL
  storedControl$parsteps <- NULL
  fit$stanfit$uncertainty$settings <- list(
    method=uncertainty,
    draws=draws,
    finishsamples=finishsamples,
    cores=cores,
    control=storedControl
  )
  if(!is.null(uncertaintyfit$scores)) fit$stanfit$subjectscores <- uncertaintyfit$scores
  message('Computing posterior approximation with ', nrow(samples), ' samples')
  fit <- ctOptimUpdateTransformed(fit, samples=samples, cores=cores)
  fit
}

#' @describeIn ctFitUncertainty The name this function shipped under in ctsem
#' 3.11.1, kept so that a call written against that release keeps meaning
#' exactly what it did. Identical to \code{ctFitUncertainty} in every other
#' respect; new code should prefer \code{ctFitUncertainty}, which is also
#' where \code{uncertainty = 'sample'} is documented.
#' @export
ctOptimUncertainty <- function(fit,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg','stored'),
  draws=c('auto','normal','empirical','imis'), finishsamples=NULL,
  cores=NULL, control=list(), verbose=0, ...){
  uncertainty <- match.arg(uncertainty)
  ctFitUncertainty(fit, uncertainty=uncertainty, draws=draws,
    finishsamples=finishsamples, cores=cores, control=control,
    verbose=verbose, ...)
}
