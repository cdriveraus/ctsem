# Cross validation for ctsem fits: K fold by refitting, or PSIS for Laplace fits

Estimates out-of-sample predictive performance, by refitting the model
with each fold of the data held out in turn or, for a julia Laplace fit,
by Pareto smoothed importance sampling without refitting.

## Usage

``` r
ctLOO(
  fit,
  folds = 10,
  cores = getOption("mc.cores", 2L),
  parallelFolds = FALSE,
  tol = NULL,
  subjectwise = ifelse(length(unique(.ctFitRowSubject(fit))) >= folds, TRUE, FALSE),
  keepfirstobs = FALSE,
  leaveOutN = NA,
  refit = TRUE,
  casewiseApproximation = FALSE,
  method = c("kfold", "psis"),
  nsamples = 1000
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
  either backend.

- folds:

  Number of cross validation splits to use – 10 folds implies that the
  model is re-fit 10 times, each time to a data set with 1/10 of the
  observations randomly removed.

- cores:

  Number of processor cores to use.

- parallelFolds:

  compute folds in parallel or use cores to finish single folds faster.
  parallelFolds will use folds times as much memory.

- tol:

  Objective tolerance for each fold's refit, by default 1e-5 on a stan
  fit. On a julia fit, by default a fold stops once a step is predicted
  to gain less than 1e-4 in log likelihood, which keeps the error in
  each held-out subject's log likelihood to about 0.014 times the square
  root of its influence (see `casewiseApproximation`).

- subjectwise:

  drop random subjects instead of data rows? With `method = 'psis'`,
  `TRUE` (the default there) leaves out one subject or group at a time,
  and `FALSE` one row at a time.

- keepfirstobs:

  do not drop first observation (more stable estimates)

- leaveOutN:

  if a positive integer is given, the folds argument is ignored and
  instead the folds are calculated by leaving out every Nth row from the
  data when fitting. Leaving 2 out would result in 3 folds (starting at
  rows 1,2,3), each containing one third of the data.

- refit:

  if FALSE, do not optimise parameters for the new data set, just
  compute the likelihoods etc from the original parameters

- casewiseApproximation:

  `TRUE` replaces each fold's refit by one Newton step from the full
  data estimate, which is much faster and close to a refit when no fold
  moves the estimate far. On a julia fit each fold's influence, the
  squared length of the step in posterior standard deviations, is
  returned in `foldinfluence`. A number instead of `TRUE` steps only the
  folds whose influence is below it, and refits the rest. Needs an
  optimised fit with a covariance. On a stan fit only subject folds and
  `TRUE` or `FALSE` are possible.

- method:

  `'kfold'` (the default) refits the model once per fold. `'psis'` is
  Pareto smoothed importance sampling leave-one-out without refitting,
  for `backend = 'julia'` fits (see Details). It needs the loo package,
  and refuses the K fold arguments `folds`, `parallelFolds`, `tol`,
  `keepfirstobs`, `leaveOutN`, `refit` and `casewiseApproximation`.

- nsamples:

  Number of importance sampling draws for `method = 'psis'`. For a
  sampled fit, its posterior draws are used instead, all of them, or
  thinned to `nsamples` if it is given. With too few draws the Pareto k
  of a heavy tailed correction is estimated too low.

## Value

For `method = 'kfold'`, a list with `foldrows`, `foldpars`,
`insampleLogLikRow`, `LogLikRowFolds`, `outsampleLogLikRow`,
`insampleLogLik`, `outsampleLogLik`, entropy and standard deviation
summaries, and `processes`, whether the folds ran in worker processes.
When a julia fit's folds are refitted, `foldoptim` has a row per fold,
with `method` (`'refit'` or `'step'`), and for a refit the optimiser's
`iterations`, `gradients`, `converged`, `gain` (the predicted gain in
log likelihood still left when it stopped) and `seconds`. A warning
counts the refits that stopped without converging. With
`casewiseApproximation`, `foldinfluence` and `foldstepped` give each
fold's influence and whether it was stepped. `scoring` is added whenever
folds were refitted, saying what was scored (see Details). An
`intoverpop = 'laplace'` fit also adds `insampleLogLikSubject` and
`outsampleLogLikSubject`. For `method = 'psis'`, an object of class
`ctLOOpsis`, with `elpd_loo`, `se_elpd_loo`, `p_loo`, `looic`,
`pointwise` (a data frame of `elpd_loo`, `lpd`, `p_loo` and `pareto_k`
per subject, group or row, and `refined` when leaving out subjects or
groups), `n_high_k`, `ndraws`, `ndropped`, `proposal_k` (the Pareto k of
the correction from the normal approximation to the posterior: one value
when leaving out subjects or groups, `NA` for a sampled fit, and one per
subject when leaving out rows), `level` and `note`.

## Details

Works for both backends. On julia, a held-out row is treated as missing,
so the filter passes through it without an update or a likelihood
contribution. For julia fits, `parallelFolds = TRUE` runs up to `cores`
folds at once, each in a background R session with its own Julia
session, started before the first fold. This needs the future package.
The sessions stay open for the next call, and
[`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)`('workers')`
ends them. Without it the folds run one at a time on `cores` threads.
Separate processes are usually faster, except when each fold takes only
seconds. `$processes` on the result says which was used.

Each fold is a separate optimisation, so on a weakly identified model,
or one where leaving out a fold leaves a parameter poorly determined,
two folds (or the two backends) can end at noticeably different raw
parameters with similar likelihoods.

**K fold on an `intoverpop = 'laplace'` fit.** The random effects are
integrated out, so what is scored depends on the type of fold. With row
folds, a held-out row is scored with each subject's random effects at
their mode given the fold's training rows, at the fold's parameters.
With subject folds, the out-of-sample quantity is each held-out
subject's Laplace marginal likelihood at the fold's parameters, returned
as `outsampleLogLikSubject`, and `outsampleLogLikRow` is `NA`, because
the marginal likelihood does not split into rows. In a multilevel model,
a subject's marginal likelihood is conditional on the rest of its
outermost group, and subjects of one group held out together share it
evenly. `insampleLogLik` is the Laplace marginal likelihood at the
estimate for subject folds, and the sum over rows at the full data modes
for row folds.

**`casewiseApproximation`.** Each fold's refit is replaced by one Newton
step from the full data estimate, `est - V %*% s`, where `s` is the
score of the held-out data at the estimate and `V` the inverse of the
negative Hessian there. This is the one-step (infinitesimal jackknife)
approximation to the refit, and its error grows faster than the step. A
fold's influence, `t(s) %*% V %*% s`, is the squared length of the step
in posterior standard deviations. Given a number, folds with an
influence at or above it are refitted, so the approximation is used only
where the step is short. `foldstepped` says which folds were stepped. On
a julia fit a stepped fold scores only its held-out rows, so its column
of `LogLikRowFolds` is `NA` elsewhere.

**`method = 'psis'`.** Pareto smoothed importance sampling (Vehtari,
Gelman and Gabry, 2017), without refitting. It estimates the Bayesian
leave-one-out predictive density, a subject's likelihood averaged over
the posterior without it. K fold instead scores each fold at its
refitted estimate, so the two can differ, most for a subject that moves
the estimate. There are two levels:

- `subjectwise = TRUE`: leave out one subject at a time, or in a
  multilevel model one group at a time, since subjects that share a
  group effect are not independent. Each subject's or group's likelihood
  has its random effects integrated out, by the Kalman filter under
  `intoverpop = 'augmented'` and by the Laplace approximation otherwise.
  On a sampled fit (`optimize = FALSE`, or
  `ctFitUncertainty(fit, 'sample')`) the population parameters are its
  posterior draws. On an optimised fit they are drawn from its normal
  approximation (`fit$estimate$cov`, so the fit must have one), widened
  by a factor of 1.5, and importance weighted to the posterior. A
  subject or group with a Pareto k above 0.7 is estimated again with
  `nsamples` more draws, centred where leaving it out moves the
  estimate, and the estimate with the lower k is kept
  (`pointwise$refined`). This is not done for a fit sampled jointly over
  parameters and random effects (`intoverpop = 'none'`). When the normal
  approximation itself is poor (`proposal_k` above 0.7), every estimate
  is unreliable, and a warning says to sample the fit and use its draws.
  Even without a warning, k is estimated too low from few draws, so a
  value near 0.7, or a posterior you expect to be skewed (an outlying
  subject, a variance near zero, binary indicators), is worth rerunning
  with larger `nsamples` or on a sampled fit.

- `subjectwise = FALSE`, on an optimised `intoverpop = 'laplace'` fit
  only: leave out one row at a time, conditional on the estimated
  population parameters, with each subject's random effects drawn from
  the normal distribution the Laplace approximation fits at their mode,
  widened by a factor of 1.5. This is leave-one-out of each row's
  density given the rows before it, which is not exact time series
  leave-one-out, and it ignores uncertainty in the population
  parameters.

Draws for which the random effects' mode cannot be found are dropped and
counted. A Pareto k above 0.7 marks an unreliable estimate for that
subject, group or row, for which K fold is the better tool.

**Laplace fits refined by quadrature.** Both methods cross-validate the
Laplace approximation, so for a fit refined with
`intoverpop = 'quadrature'` (in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) or
[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md))
they start from the Laplace estimate,
`fit$laplace$correction$laplace_estimate`, and the folds are not
refined. `scoring` says so.

**Fold refits are not checked.** Each refitted fold skips the final
convergence check and Hessian of a full fit
(`optimcontrol$certify = FALSE`), which would cost each fold a Hessian.
`scoring` notes this whenever folds are refitted.

## References

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27, 1413–1432.

## Examples

``` r
# \donttest{
ctLOO(ctstantestfit)
# }
```
