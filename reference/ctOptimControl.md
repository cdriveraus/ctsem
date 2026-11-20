# Optimiser settings for ctFit

The settings
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) reads
from `optimcontrol`. Each has the same meaning on both backends. A
setting the backend does not have is an error, unless its value asks for
what that backend does anyway (`gradient = 'adjoint'` on stan,
`nsubsets = 1` on julia), and so is a misspelled name. Settings left
unset take the backend's default.

## Stopping

- `tol`:

  Both. Stop when the objective changes by less than this. Default
  `1e-8` on stan, off on julia.

- `g_tol`:

  Both. Stop when the l2 norm of the gradient falls below this. Default
  `1e-8` on julia, off on stan.

- `x_tol`:

  Both. Stop when the parameter step falls below this. Off by default.

- `maxiter`:

  Both. Iteration cap. Default 1000 on julia, mize's on stan.

- `lbfgs_scores`:

  Julia. Default TRUE. Keeps a second model of the curvature beside
  L-BFGS, built from the outer product of the subjects' scores, and
  takes each step with whichever model predicted the last step better.

- `lbfgs_local`:

  Julia. Default TRUE. With `intoverpop = 'laplace'`, parameters that
  only a few subjects or groups inform (a study's, a group's, or an item
  few answered) start L-BFGS from their exact curvature, computed from
  those subjects alone.

- `lbfgs_memory`:

  Both. How many curvature pairs L-BFGS keeps. Default on julia one per
  free parameter, between 20 and 1000, and mize's default on stan.

- `stallretries`:

  Both. Default 2. How many times a fit that stopped short of a maximum
  is continued. Stan restarts from new random values when the gradient
  exceeds `stalltol` per data point. Julia resumes from a better point
  found by its end-of-fit checks, or from where it stopped, and keeps
  the resumed run only if it improves the objective. 0 turns it off.

## Starting values

- `initsd`:

  Both. Default 0.01. The sd of the random start, on the raw scale.

- `datastart`:

  Julia. Default TRUE. Starting values for the diagonals of DRIFT,
  DIFFUSION, T0VAR and MANIFESTVAR are set from the scale of the data,
  using `log(1 + y)` for counts, and without using binary or ordinal
  indicators. If the objective cannot be evaluated there, the fit starts
  from the default values and says so. Supplied `inits` are always used
  as given.

- `carefulfit`:

  Both. Default TRUE. A first pass with ctsem's `normal(0,1)` priors on
  every raw parameter, to find starting values, after which the fit's
  own objective is maximised. Skipped with `priors = TRUE`, which has
  those priors already, and when `inits` are supplied. On julia it runs
  until it converges, for at most 200 iterations, and a number sets that
  limit. It keeps a fit from starting at a boundary, such as a
  correlation of one or a variance of zero. `fit$optim$carefulfit`
  records whether it ran.

## The optimiser

- `stochastic`:

  Both. Start with ctsem's stochastic-gradient optimiser and finish with
  L-BFGS. `'auto'` does so above 50 parameters. Default TRUE on stan,
  FALSE on julia. It can be much faster than L-BFGS alone on large,
  poorly conditioned models, and can also take a fit further toward a
  degenerate limit, such as a variance going to zero.
  `fit$optim$sgd_iterations` records it.

- `gradient`:

  Julia. `'adjoint'` (reverse mode, the default) or `'forward'`. Both
  give the same gradient, and the cost of the adjoint does not grow with
  the number of parameters.

- `batch`:

  Julia. Default TRUE, or FALSE when `inits` are supplied. Start on a
  random subset of the subjects and add more as the fit needs them.
  Needs at least 80 subjects or top-level groups. `TRUE` also batches a
  fit started from supplied `inits`.

- `coreset`:

  Julia. Default `'auto'`. Run parts of the optimisation on a weighted
  subset of the subjects (or of the top-level groups, with
  `intoverpop = 'laplace'`), chosen to reproduce the full gradient, and
  check each part against the full gradient. The fit finishes on all the
  data, so the estimates and standard errors use every subject. `'auto'`
  uses subsets for data with at least 200 subjects or top-level groups
  and 10,000 data cells (rows times manifest variables). On smaller data
  they save no time, and can change which of several maxima a fit
  reaches. `TRUE` uses them whenever they are expected to save
  gradients, `FALSE` never. `fit$optim$coreset_sizes` records them.

- `explosive_forward`:

  Julia. Default FALSE. For a subject whose predictions grow more than
  100-fold between observations (explosive dynamics over a long gap),
  compute its gradient in forward mode as a cross-check. The default
  reverse mode is accurate there too. The first such subject in a
  session compiles the filter again, which takes minutes for a large
  model. `fit$optim$explosive_passes` counts such subjects, and the fit
  warns when any is explosive at the estimate.

- `lbfgs_diagonal`:

  Julia. Default TRUE. L-BFGS learns a separate scale for each parameter
  from the curvature (Gilbert and Lemarechal). FALSE uses one scale for
  all, which can be much slower on large models. Like `stochastic`, it
  can take a fit further toward a degenerate limit.

- `lbfgs_gll`:

  Julia. Default 0. A window `W`: a step is accepted if it improves on
  the worst of the last `W` objective values (Grippo, Lampariello and
  Lucidi), which lets the fit cross a curved valley.

- `lbfgs_nonmonotone`:

  Julia. Between 0 (the default) and 1. A step is accepted against a
  running average of recent objective values (Zhang and Hager).

- `initial_alpha`:

  Julia. Default 0.1. The length of the first trial step, in raw units.

- `precondition`:

  Julia. Default FALSE. TRUE sets L-BFGS's initial scaling of each
  parameter from the model's transforms, instead of learning it from the
  curvature. This can make fits whose intercepts start far from the data
  end in worse optima.

- `saem`:

  Julia, `intoverpop = 'laplace'`. Default FALSE. Start with SAEM, which
  samples the random effects from their conditional distribution instead
  of integrating them with the Laplace approximation. The usual
  optimiser then finishes on the Laplace objective, so the estimate is
  usually the same as without SAEM, but SAEM can lead to a different
  maximum when the objective has more than one. Because SAEM's result
  has Monte Carlo error, the fit also runs once from its usual start
  without SAEM, and keeps the maximum with the higher log posterior
  after the Laplace approximation is corrected by quadrature.
  `fit$optim$saem_guard` records both. For a more accurate estimate than
  the Laplace approximation gives, sample with `optimize = FALSE`. A
  number sets the maximum number of SAEM iterations (TRUE is 10000).
  `fit$optim$saem_iterations`, `saem_settled`, `saem_trend`,
  `saem_chains`, `saem_acceptance` and `saem_trace` record the SAEM run.

- `saem_proposal`:

  Julia. SAEM's proposal, `'rw'` or `'laplace'`.

- `restarts`:

  Julia. Default 0. The number of random restarts for a fit that has not
  converged after its other corrections, typically because the
  likelihood has more than one maximum. Each starts from the fit's own
  start plus normal noise with sd `restartsd` (default 1) on the raw
  scale, and the best result is kept. Restarts run in worker processes
  when `cores > 1` and the fit took more than 30 seconds. Not used when
  `inits` or `maxiter` were supplied. `fit$optim$restarts` records each
  start.

## Newton steps and convergence checks (julia)

- `newton`:

  Default TRUE. Take Newton steps near the maximum, with the exact
  Hessian, updated after each step and computed again when progress
  slows. The Hessian at the end is used for the convergence check and
  the standard errors. The steps also move away from a saddle point
  along the direction of negative curvature, and test directions with no
  curvature. With `intoverpop = 'laplace'` they run only when the fit is
  certified (`certify`). The value may name the curvature the steps use:
  `'chord'` (the default except as below), `'exact'` (a new Hessian
  whenever a step makes slow progress), `'subset'`, or `'krylov'` (the
  outer product of the subjects' scores, corrected by Hessian-vector
  products at each step, with the exact Hessian computed at the end).
  `'krylov'` is the default with `intoverpop = 'laplace'` when the exact
  Hessian would cost at least 200 gradients. `fit$optim$newton_steps`
  records the steps, and `newton_krylov_products` the gradients the
  Krylov steps used.

- `certify`:

  Default TRUE. Check from the curvature that the fit has converged.
  `estonly` also skips it.

- `gaptol`:

  Default `5e-5`. The fit continues while the curvature at the estimate
  predicts more than this gain in the objective. The default corresponds
  to 0.01 standard errors from the maximum, close enough that the
  Hessian at the estimate is the Hessian at the maximum. A smaller
  value, such as `1e-6`, holds the estimate closer still.

- `hessianreuse`:

  Default 0.3. How far, in standard errors, the estimate may be from
  where a Hessian was computed for that Hessian still to give the
  standard errors. A Hessian from this far away adds roughly this
  fraction to the error of the standard errors. Computing a new one
  costs a gradient per parameter, or two with `intoverpop = 'laplace'`.
  0 always computes it at the estimate.

- `innergaptol`:

  The gain that a less reliable prediction (from L-BFGS, or from a
  Newton step on a Hessian computed elsewhere) must fall below for the
  fit to stop. By default a hundredth of `gaptol`. A gain predicted by
  the exact Hessian at the estimate is compared with `gaptol` itself.

- `overshoot`:

  `'magnitude'` (the default), `'saturation'` or `'off'`. Which
  directions are pulled back after a fit, to test whether it stopped at
  a maximum.

- `stallwindow`, `stallratio`, `stallcooldown`:

  The stall check (defaults 80, 0.001 and 30). The number of iterations
  over which a fit's gain is measured (a gain of half a unit of log
  likelihood or less starts a search for a transform that has gone
  flat), how much of its responsiveness a transform must have lost to
  count as the reason the fit stopped, and how many iterations the check
  waits after finding a fit slow but not stuck.

- `escapesaturated`:

  Default FALSE. Refit a fit that stalled on a saturated transform, with
  that parameter reset to zero, and keep the better result.

- `escapepin`:

  Default FALSE. Hold the reset parameters fixed while the others are
  optimised, then free them.

## Laplace approximation (julia, `intoverpop = 'laplace'`)

- `laplace_inner_maxiter`, `laplace_inner_tol`:

  Defaults 200 and `1e-10`. Iteration limit and tolerance for finding
  each subject's (or group's) random-effect mode. Raise the first when
  the fit reports modes that did not converge. The tolerance applies to
  the gradient, and is at least `1e-10 * (1 + abs(value))`.

- `laplace_onestep`:

  `FALSE` (the default) or `TRUE`. Where the Kalman filter can integrate
  a subject's random effects on initial means, intercepts or predictor
  effects exactly, given the subject's other random effects, it does,
  and the Laplace approximation covers only the others. `TRUE` uses one
  Laplace approximation over all of them instead. With random dynamics
  (a random drift or diffusion) the default is more accurate. Models
  whose dynamics depend on the state always use the single
  approximation, and
  [`ctModelImprove()`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
  and
  [`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
  use it too when comparing such a model with one that does not.

- `laplace_soft_directions`:

  `'quadrature'` (the default) or `'logdet_floor'`. How to integrate a
  subject's random effects in a direction where their posterior is
  flatter than the prior, where the plain Laplace approximation
  overstates the likelihood. `'quadrature'` integrates the flattest
  direction with a three-point quadrature when the curvature has an
  eigenvalue below 0.7, blending into plain Laplace between 0.2 and 0.7.
  `'logdet_floor'` floors the log determinant at zero. Recorded in
  `fit$laplace$soft_directions` and used by every function that rebuilds
  the objective, and `fit$laplace$conditioning` counts such subjects.
  Also applies when sampling with `intoverpop = FALSE` and
  `sampleControl$placement = 'fit'`.

## Uncertainty

- `uncertainty`:

  Both. Default `'hessian'` (`'sandwich'` with subject weights on
  julia). The uncertainty method computed when the optimiser ends, as in
  [`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md).
  The default is meant for building a model, with the checks below. For
  intervals to report, sample with `ctFitUncertainty(fit, 'sample')`.
  `'sample'` is not accepted here. Use `optimize = FALSE`, or
  `ctFitUncertainty(fit, 'sample')` after the fit.

- `checks`:

  Julia. Default `TRUE`, which runs both checks below, as the `checks`
  argument of
  [`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
  does. `FALSE` runs neither, and `'randomeffects'` or `'intervals'`
  runs one. Neither runs with `estonly`. The random-effects check (with
  `intoverpop = 'laplace'`) scores the Laplace approximation for every
  subject at draws of its random effects, and warns if integrating them
  exactly would move an estimate by 0.1 standard errors or more. The
  warning gives the size of the move when it is below one standard
  error, and how far the Laplace log likelihood is from the integrated
  one, and points to `ctFitUncertainty(fit, intoverpop = 'quadrature')`
  and to sampling. When the check runs, `fit$estimate$loglik` is the
  Laplace log likelihood corrected by the check, with `loglik_se` its
  standard error, `loglik_laplace` the uncorrected value, and
  `loglik_method` `'checked'`. Models whose random effects enter
  linearly pass at once. If the optimiser stopped where the Laplace
  approximation for some subjects jumps, the fit warns instead and names
  them. `fit$laplace$correction` records the check. The interval check
  (with `uncertainty = 'hessian'`) tests whether the normal
  approximation from the Hessian matches the posterior, and says so in a
  message and in [`summary()`](https://rdrr.io/r/base/summary.html)
  where it does not. With `intoverpop = 'quadrature'` it uses the
  refined fit, and can then take longer than the refinement.
  `checks = 'randomeffects'` leaves it out.

- `uncertaintyDraws`, `finishsamples`, `uncertaintyControl`:

  Both. The `draws` argument of `ctOptimUncertainty` (default `'auto'`),
  and the `finishsamples` (1000) and `control` arguments of
  `ctFitUncertainty`.

- `accuracy`, `sampleControl`:

  Julia. A target accuracy in posterior standard deviations, off by
  default. After the uncertainty method has run, its estimate of its own
  error is compared with this, and if it falls short, or gives none, the
  fit is sampled to that accuracy (`control$accuracy` in
  `ctFitUncertainty`). `sampleControl` holds the other sampler settings
  for that run.

- `estonly`:

  Both. Default FALSE. Point estimates only, without a Hessian, so there
  are no standard errors and no convergence check. With
  `intoverpop = 'laplace'` there are also no Newton steps and no
  random-effects check.

## Progress and output (julia)

- `callback`:

  A function called while the fit runs, with arguments
  `(iteration, total, objective, gradient_norm, parameters)`, where
  `parameters` is the raw parameter vector. A function with only the
  first four arguments is called with four. It is called at regular time
  intervals and once at the end, for example to show progress in a front
  end. An error inside it disables it, with a warning. `fit$optim$trace`
  holds every iteration after the fit, and
  [`ctTracePlot`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md)
  plots it. To save a checkpoint to resume from, use `checkpoint`.

- `checkpoint`:

  A file path. While the fit optimises, its current point is saved there
  every `checkpoint_every` seconds (default 60) and at the end of each
  stage, together with the random-effect modes (with
  `intoverpop = 'laplace'`) and the smoothing parameters of any smooth
  terms. Each exact Hessian the fit computes is saved beside it
  (`<checkpoint>.hessian`). If the file already exists, the fit starts
  from it, so running the same call again resumes an interrupted fit,
  including its modes, smoothing parameters and Hessian. A file written
  for a different model stops the fit instead of being overwritten.

- `progress`:

  Force the progress line on or off.

- `tipredMissing`:

  How missing values of time independent predictors are handled. Each
  missing value gets a normal distribution from a regression on the
  other variables. `'impute'` sets it to the mean of that distribution,
  `'sample'` treats it as a parameter with that distribution as its
  prior, and `'auto'`, the default, imputes when optimizing and samples
  when sampling. `ctFitUncertainty(fit, 'sample')` uses the fit's
  setting, and its `control$tipredMissing` overrides it.

- `tipredMissingIncludeOutcome`:

  Default TRUE. Whether the outcome informs the imputation of missing TI
  predictors.

- `saveEffects`:

  A sampling setting, read here for older scripts. Use
  `sampleControl$saveEffects`.

## Stan only

- `nsubsets`, `subsamplesize`, `lproughnesstarget`,
  `stochasticTolAdjust`, `parsteps`:

  Tuning of stan's stochastic-gradient optimiser and its stepwise phase.

- `stalltol`:

  Default 0.01. The gradient per data point above which stan calls a fit
  stalled and restarts it (see `stallretries`).

## See also

[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md),
[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md),
[`ctTracePlot`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md).
