# Improve a model step by step, from score tests and predictive checks

Starting from a fitted model, repeatedly screens it with
[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
refits it with the most strongly indicated relaxations one at a time,
starting from the current estimates, and keeps the best of them if the
refit's likelihood confirms it. The search stops when no relaxation is
indicated and
[`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
finds the model adequate. When nothing is indicated but the model is not
adequate, it tries structural changes in the same way: a slow process
driving a latent
([`ctModelSlowComponent`](https://cdriveraus.github.io/ctsem/reference/ctModelSlowComponent.md)),
a cycle such as a daily rhythm
([`ctModelPeriodicComponent`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md)),
a process of its own for an item
([`ctModelUniqueProcess`](https://cdriveraus.github.io/ctsem/reference/ctModelUniqueProcess.md)),
a higher order latent process
([`ctModelHigherOrder`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md)),
and any models given in `candidates`. These are ranked by the adequacy
check's evidence for them: a slow process for a latent whose updates
stay correlated over time gaps (the `latent` family), a cycle at the
period the `periodic` family points to, and a process of its own for an
item whose errors stay correlated over time (the `residuals` family).

## Usage

``` r
ctModelImprove(
  fit,
  covariates = NULL,
  relax = c("heterogeneity", "covariates", "statedependence", "volatility", "fixed",
    "level", "slow", "unique", "periodic", "order"),
  changes = NULL,
  candidates = NULL,
  alpha = 0.05,
  maxsteps = 10L,
  ntry = 3L,
  criterion = c("LR", "AIC", "BIC"),
  adequacy = TRUE,
  prune = TRUE,
  close = NULL,
  holdout = 0,
  form = "auto",
  periods = NULL,
  nsamples = 100L,
  projection = "auto",
  cores = getOption("mc.cores", 2L),
  checkpoint = NULL,
  verbose = TRUE
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
  `backend = 'julia'`.

- covariates:

  Optional subject level covariates the model may add, as for
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md).

- relax:

  Which relaxations may be made.

- changes:

  Optional rows of
  [`ctModelProposals`](https://cdriveraus.github.io/ctsem/reference/ctModelProposals.md).
  The first round refits with these, in the order given, instead of
  screening, and keeps the first that improves the fit by `criterion`.
  With `maxsteps = 1` this is one step of a search made by hand. Later
  rounds screen as usual.

- candidates:

  Optional list of alternative models (from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md))
  to try when the screen finds nothing but the checks fail.

- alpha:

  Level for the screen (Benjamini-Hochberg within each family), for the
  likelihood ratio decision, and for adequacy.

- maxsteps:

  Maximum number of refit rounds.

- ntry:

  Maximum number of proposals from the score tests refitted in a round,
  one at a time, the strongest first, keeping the first that improves
  the fit. Proposals from the adequacy check are tried in the same way
  until one improves the fit or none is left. Pruning and held-out
  confirmation refit this many at a time.

- criterion:

  `'LR'`, `'AIC'` or `'BIC'`.

- adequacy:

  Use
  [`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
  to decide when to stop.

- prune:

  After the search, try dropping each kept change in turn, and drop it
  if the likelihood does not need it given the others. A change kept
  early can be made unnecessary by one kept later.

- close:

  Optional root mean square residual correlation regarded as a close
  enough fit to stop at (see Selection and confirmation).

- holdout:

  Share of subjects, between 0 and 0.5, held out of the search to
  confirm its changes on.

- form:

  How covariates are screened, as for
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md).
  `'threshold'` lets the search add a covariate as a step.

- periods:

  Periods for the adequacy check's `periodic` family, as for
  [`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md).
  By default a range is scanned.

- nsamples:

  Number of generated datasets for the adequacy checks.

- projection:

  How the adequacy checks account for estimation, as for
  [`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md).

- cores:

  Cores for fitting and testing. With two or more cores, once a refit
  has taken longer than five minutes, the refits for pruning and for
  held-out confirmation run in parallel worker processes, `ntry` at a
  time, each with an equal share of `cores`.
  `options(ctsem.improve.processes = FALSE)` keeps them in this session.

- checkpoint:

  Optional file the search saves its state to after every round and
  after pruning. Called again with the same file, fit and settings, the
  search resumes after the last round saved. A file from a different
  search is refused. With `verbose`, the refit in progress is also saved
  every ten minutes (`<checkpoint>.refit.rds`), with the round, the
  change, the iteration and the raw parameters reached.

- verbose:

  Report progress: each round's steps, each refit's time, and while a
  refit runs, its iteration, objective and gradient every ten minutes.

## Value

An object of class `ctModelImprove`, containing: `fit`, the final fit;
`history`, every refit tried, with the round, the relaxation, its
screening p-value, the change in log likelihood and parameters, the
criterion's verdict and whether it was kept (a row with `family`
`'prune'` says whether the change was dropped, and a `'rank'` row
whether the lower rank was kept); `screens`, the score tests of each
round; `adequacy`, the final checks; `confirmation`, with `holdout`,
each kept change's likelihood ratio test on the held-out subjects; and
`stopped`, why the search ended: `'adequate'` (no adequacy family
detected misfit), `'close fit'` (those that did were shown to be within
`close`), `'not adequate, no candidate improved'` or
`'not adequate, nothing left to try'` (the checks still detect misfit),
`'nothing indicated'` (with `adequacy = FALSE`, the screen found
nothing) or `'maxsteps'`. A refit that failed, or ended below the model
it extends, is in `history` with a note. With `holdout`, `fit` is the
model with the confirmed changes refitted to every subject, and
`selected` is the search's own choice fitted to the subjects it searched
on.

## Details

The relaxations, in terms of the model:

- heterogeneity:

  `indvarying = TRUE` for a parameter and the `T0MEANS` parameters. For
  a fixed `MANIFESTMEANS` cell, a parameter that varies around the fixed
  value as its stated mean.

- covariates:

  A covariate added as a time independent predictor, with effects on the
  parameters its tests pick out. For a predictor the model already has,
  effects on further parameters its tests pick out.

- statedependence:

  A term `d * eta_k * eta_m` in a latent's drift, written into `CINT`
  with `d` a `PARS` parameter, so that a drift coefficient varies with a
  state.

- volatility:

  A latent's innovation sd scaled by `exp(v * eta_m)`, so that the size
  of the innovations depends on a state.

- fixed:

  A cell fixed at zero is freed: a drift cross effect, a diffusion
  covariance or a cross loading.

- level:

  A parameter of the dynamics or error variances that depends on a
  subject's level (a random intercept) is allowed to vary over subjects,
  correlated with the intercepts (see the level tests of
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)).

- slow:

  A slow process driving a latent through its drift, with its effect
  starting at one and its persistence at a tenth of the latent's. If the
  model already has a slow process for another latent, letting that
  process drive this latent too is tried first.

- periodic:

  A cycle of the period the adequacy check points to, driving every
  latent that has innovations. Each subject's phase is set from the
  clock at its first row and can deviate from it. The phase predictors
  are added to the data.

- unique:

  A persistent process of its own, loading on one item. Where the
  `pairs` family finds two items' errors correlated within an occasion,
  a process shared by the two (`ctModelUniqueProcess(shared = TRUE)`) is
  tried too, treating such correlation as a process faster than the
  occasions. Like the other structural changes, this is tried only once
  the score tests indicate nothing.

- errorcov:

  Not by default. A measurement error covariance the score test flags is
  freed, instead of proposing a shared process. This is for designs with
  a reason for correlated errors, such as a method the items share.

- order:

  A latent process made second order.

On a reduced rank level (`poprank`), a heterogeneity result for a
parameter that already varies through the shared dimensions gives it a
variance of its own (`popunique`). Without reduced rank, a significant
`dimension` test makes the parameters in its direction vary together.
Relaxations the score tests cannot tell apart (the same `class`) are not
refitted in the same round.

Where the shared dimensions are in blocks (`popblocks` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)), the
search acts on
[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)'s
block tests. A block that needs a further dimension is given one,
starting along the direction the test found (without blocks, the level's
rank goes up by one). When the dimensions of two blocks covary, the
blocks share them, and whichever block needs fewer new loadings to reach
the other's effects is extended to them. A parameter that a change makes
individually varying joins the block of the effects it should covary
with (the initial means for a heterogeneity result, the intercept for a
level result), unless a block already includes it. A block of full rank
gains a dimension with it.

Each random effect the search adds brings correlations with the others,
and with many random effects and few subjects the data no longer
determine them. When a kept change leaves the population covariance
poorly identified (`fit$identifiability`), the covariance is refitted at
one rank lower (`poprank`), starting from its leading dimensions, and
the lower rank is kept unless the likelihood ratio test needs the
dropped dimension. This repeats until the fit is identified. Later
random effects then load on the shared dimensions, and later
heterogeneity results give own variances. This is for single level
models only.

When the search stops, each change kept is tested again by refitting
without it (`prune`), and dropped if the others make it unnecessary.

A relaxation is kept when its refit improves on the current model by
`criterion`: a likelihood ratio test against a chi-square on the
parameters the change adds, at `alpha` divided by the number of refits
tried at that step (`'LR'`), or a lower AIC or BIC (with subjects as the
sample size). The chi-square reference is only approximate for these
changes. It is conservative for a new random effect, whose variance
starts at zero, and can be wrong in either direction for a new process
or cycle. Screening, rounds and pruning are further choices it does not
count, so the p-values describe the steps the search took and do not
test the final model. A refit that ends below the model it extends has
not reached an optimum, and is recorded as such.

## Error rates

Because the screen and the refit use the same data, the error rate is
close to `alpha` for each decision, and grows with the number of
families screened. In simulations with nothing missing, a one process
model (80 subjects) kept a change in 5 per cent of searches at
`alpha = 0.05` and in none at `alpha = 0.01`, while a model of two
processes with three items each, where state dependence and volatility
are also screened in each process, kept a change in 4 of 20 searches at
`alpha = 0.05`. For larger models `alpha = 0.01` is the safer choice.
`criterion = 'BIC'` was not stricter than `'LR'` at these sizes. Small
effects can be missed or replaced: a drift that varied over subjects was
not found at 80 subjects, and a random diffusion or measurement error
variance was sometimes added instead.

## Selection and confirmation

The result is a model chosen by the data it is fitted to, and its
standard errors and tests do not account for that choice (MacCallum,
Roznowski & Necowitz, 1992). `holdout` sets aside a share of the
subjects. The search runs on the rest, and each change it keeps is then
tested on the held-out subjects alone, by refitting the final model
without it there. `confirmation` says which changes held, and the final
model, refitted to every subject, keeps only those. Every step, and
every refit tried and rejected, is in `history`.

With many observations, the checks can reject a model for departures
that are small. With `close`, the search also stops when every failing
check is a residual or within row correlation family whose size is shown
to be below `close` (see
[`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)).

## Laplace fits

On a Laplace fit, a change that makes the dynamics depend on the state
is refitted with the one-step Laplace approximation
(`optimcontrol$laplace_onestep` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)), and
where the model has random dynamics, it is compared with the current
model refitted the same way. To compare the models with their random
effects integrated more exactly, start the search from a fit with
`intoverpop = 'quadrature'`, so that every refit is refined in the same
way.

## See also

[`vignette('score-tests')`](https://cdriveraus.github.io/ctsem/articles/score-tests.md).
