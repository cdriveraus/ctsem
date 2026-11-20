# Optimize / importance sample a stan or ctStan model.

Optimize / importance sample a stan or ctStan model.

## Usage

``` r
stanoptimis(
  standata,
  sm,
  init = "random",
  initsd = 0.01,
  estonly = FALSE,
  tol = 1e-08,
  stochastic = TRUE,
  priors = TRUE,
  carefulfit = TRUE,
  stallretries = 2,
  stalltol = 0.01,
  uncertainty = c("hessian", "checked", "is", "bootstrap", "fullbootstrap", "sandwich",
    "opg"),
  uncertaintyDraws = "auto",
  uncertaintyControl = list(),
  subsamplesize = 1,
  parsteps = c(),
  plot = FALSE,
  finishsamples = 1000,
  lproughnesstarget = 0.2,
  verbose = 0,
  cores = 2,
  matsetup = NA,
  nsubsets = 1,
  stochasticTolAdjust = 1000,
  g_tol = NULL,
  x_tol = NULL,
  maxiter = NULL,
  lbfgs_memory = NULL
)
```

## Arguments

- standata:

  list object conforming to rstan data standards.

- sm:

  compiled stan model object.

- init:

  vector of unconstrained parameter values, or character string 'random'
  to initialise with random values very close to zero.

- initsd:

  positive numeric specifying sd of normal distribution governing random
  sample of init parameters, if init='random' .

- estonly:

  if TRUE,just return point estimates under \$rawest subobject.

- tol:

  objective tolerance. The optimizer stops when the objective changes by
  less than this between evaluations.

- stochastic:

  Logical. Use stochastic gradient descent as main optimizer. Always
  finishes (double checks) with mize (bfgs) optimizer.

- priors:

  logical. If TRUE, a priors integer is set to 1 (TRUE) in the standata
  object – only has an effect if the stan model uses this value.

- carefulfit:

  Logical. If TRUE, priors are always used for a rough first pass to
  obtain starting values when priors=FALSE

- stallretries:

  Integer. Number of times to restart the optimizer from new random
  values when it stops somewhere that is not a maximum (see `stalltol`).
  Only applies when `init='random'`. With supplied inits a failed fit is
  reported but not retried. Set to 0 to disable retrying.

- stalltol:

  Gradient per data point above which the optimizer is taken to have
  stopped short of convergence. Rough likelihoods, such as with binary
  indicators, can leave the optimizer unable to find an improving step,
  which it reports as convergence at the starting values. Converged fits
  are usually near 1e-5 per data point, and stalled ones at 1 or more.

- uncertainty:

  Character string selecting the optimized-fit uncertainty
  approximation. Options are `'hessian'`, `'checked'`, `'is'`,
  `'bootstrap'`, `'fullbootstrap'`, `'sandwich'`, and `'opg'`.

- uncertaintyDraws:

  Character string controlling approximate raw-parameter draws from the
  approximate uncertainty. `'auto'` uses empirical draws for
  `uncertainty='bootstrap'` or `uncertainty='fullbootstrap'` and normal
  draws otherwise. `'normal'` draws from a multivariate normal using the
  selected covariance, `'empirical'` uses empirical draws when
  available, and `'imis'` runs importance sampling.

- uncertaintyControl:

  List of method-specific options passed to
  [`ctOptimUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
  internals. Score-based methods use subject-level score contributions
  when there are at least two subjects; single-subject models warn and
  use case-level contributions. Score-based methods warn when there are
  fewer than ten independent subjects or no more score rows than raw
  parameters. Full bootstrap requires at least two subjects and warns
  below ten independent subjects. Bootstrap-style methods require at
  least two returned samples / refits.

- subsamplesize:

  value between 0 and 1 representing proportion of subjects to include
  in first pass fit.

- parsteps:

  ordered list of vectors of integers denoting which parameters should
  begin fixed at zero, and freed sequentially (by list order). Useful
  for complex models, e.g. keep all cross couplings fixed to zero as a
  first step, free them in second step.

- plot:

  Logical. If TRUE, plot iteration details. Probably slower.

- finishsamples:

  Number of samples to draw (either from hessian based covariance or
  posterior distribution) for final results computation.

- lproughnesstarget:

  target log posterior roughness for stochastic optimizer (suggest
  between .05 and .4).

- verbose:

  Integer from 0 to 2. Higher values print more information during model
  fit – for debugging.

- cores:

  Number of cpu cores to use, should be at least 2.

- matsetup:

  subobject of ctStanFit output. If provided, parameter names instead of
  numbers are output for any problem indications.

- nsubsets:

  number of subsets for stochastic optimizer. Subsets are further split
  across cores, but each subjects data remains whole – processed by one
  core in one subset.

- stochasticTolAdjust:

  Multiplier for stochastic optimizer tolerance.

- g_tol:

  gradient tolerance. Stop when the l2 norm of the gradient falls below
  this. `NULL` (the default) leaves it off, so `tol` is the stopping
  rule.

- x_tol:

  parameter-step tolerance. Stop when the update to the parameter vector
  is smaller than this. `NULL` leaves it off.

- maxiter:

  maximum optimizer iterations. `NULL` leaves the existing caps in place
  (99999 for the bfgs optimizer, 5000 for the stochastic one).

- lbfgs_memory:

  number of curvature pairs L-BFGS keeps. `NULL` uses 100.

## Value

list containing fit elements
