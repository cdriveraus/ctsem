# Write a folder of summary, diagnostic and interpretation output for a fit

Runs the package's reporting and fit-checking functions over a fitted
model and writes the results into one folder, with a markdown index that
says what each file is and how to read it. Works for either backend.
Components a backend or a model does not support are listed in the index
as not produced, and the run continues.

## Usage

``` r
ctReport(
  fit,
  folder = "ctReport",
  components = "default",
  nsamples = 200,
  cores = getOption("mc.cores", 2L),
  quantiles = c(0.025, 0.5, 0.975),
  times = "auto",
  timeinterval = 1,
  nsubjects = 6,
  lags = 0:3,
  folds = 10,
  profilepoints = 9,
  verbose = TRUE
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)
  (class `ctStanFit` or `ctJuliaFit`).

- folder:

  Directory to write into; created if it does not exist. Files from a
  previous `ctReport` run in the same folder are removed first, so the
  folder always matches its index.

- components:

  Either one of `'quick'` (nothing that needs generated data or further
  model evaluations), `'default'`, `'all'` (adds cross validation and
  LaTeX), or a character vector of component names: `summary`,
  `parmatrices`, `identification`, `profile`, `discretepars`, `network`,
  `predictions`, `residuals`, `adequacy`
  ([`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)),
  `scoretests`
  ([`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
  julia fits), `covcheck`, `postpred`, `tipredeffects`, `crossval`,
  `equations`.

- nsamples:

  Draws used by the components that sample from the fit's uncertainty,
  and by data generation for `adequacy`, `covcheck` and `postpred`.

- cores:

  Cores for the components that can use more than one (`discretepars`,
  `adequacy`, `scoretests`, `covcheck`, `postpred`, `crossval`).

- quantiles:

  Quantiles for the parameter matrices and the discrete parameter plots.

- times:

  Interval grid for `discretepars`. `'auto'` uses zero to the 90th
  percentile of the observed per-subject time span.

- timeinterval:

  Time interval at which the networks, the discrete-time parameter
  matrices and the covariate effects are reported.

- nsubjects:

  Number of subjects plotted by `predictions`, chosen as those with the
  most observations.

- lags:

  Lags for the lagged covariance check.

- folds:

  Folds for `crossval`.

- profilepoints:

  Grid points per parameter in the likelihood profile.

- verbose:

  Report each component as it runs.

## Value

Invisibly, a list with the folder and, per component, whether it was
written, skipped or failed. Called for the files it writes.

## Details

The slow parts are `adequacy`, `covcheck` and `postpred`, which generate
data from the fit (once, shared between them), `scoretests`, which
evaluates the model extended with new terms, and `crossval`, which
refits. `components = 'quick'` leaves them all out.

`profile` evaluates the log probability along each raw parameter,
holding the others at the estimate, over a grid in units of the standard
error implied by that parameter's own curvature. It flags parameters
whose slice peaks away from the estimate (the optimiser did not reach
the maximum), whose shape does not match what a standard error assumes,
and whose reported interval is much wider than their own curvature
implies (they lie on a ridge with other parameters). With the other
parameters held fixed it is cheap, but it is not a profile likelihood
([`ctFitProfile`](https://cdriveraus.github.io/ctsem/reference/ctFitProfile.md)).
`identification` names the ridge directions.

## Which function checks what

- Can the data inform every parameter?:

  [`ctIdentify`](https://cdriveraus.github.io/ctsem/reference/ctIdentify.md)
  before fitting; `fit$identifiability` and
  [`ctFitProfile`](https://cdriveraus.github.io/ctsem/reference/ctFitProfile.md)
  after.

- Did the optimiser, or the sampler, get there?:

  `print(fit)` says whether it converged, `summary(fit)` adds R-hat and
  effective sample sizes for a sampled fit, and
  [`ctTracePlot`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md)
  draws the optimiser's trace or each chain's draws.

- How far did the data move each parameter from its prior?:

  [`ctPlotPosterior`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md).

- Does data generated from the fit look like the data?:

  [`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md),
  and
  [`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
  for the lagged covariances, split by a group or predictor if wanted.

- Are the one-step-ahead errors what the model says they are?:

  [`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
  tests them by family, allowing for the estimation;
  [`ctResiduals`](https://cdriveraus.github.io/ctsem/reference/ctResiduals.md)
  and
  [`ctACFresiduals`](https://cdriveraus.github.io/ctsem/reference/ctACFresiduals.md)
  show them.

- What would relaxing the model change?:

  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
  without refitting;
  [`ctModelProposals`](https://cdriveraus.github.io/ctsem/reference/ctModelProposals.md)
  lists the changes and
  [`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
  refits them.

- Does a richer model fit better, or predict better?:

  [`anova`](https://cdriveraus.github.io/ctsem/reference/anova.ctJuliaFit.md)
  for nested models,
  [`ctLOO`](https://cdriveraus.github.io/ctsem/reference/ctLOO.md) for
  out-of-sample prediction.

- How much do the julia engine's approximations cost?:

  [`ctLaplaceCheck`](https://cdriveraus.github.io/ctsem/reference/ctLaplaceCheck.md)
  for the Laplace approximation of random effects,
  [`ctParticleLik`](https://cdriveraus.github.io/ctsem/reference/ctParticleLik.md)
  for the filter's Gaussian approximation of the states.

- Do the intervals cover as often as they claim?:

  [`ctCoverageCheck`](https://cdriveraus.github.io/ctsem/reference/ctCoverageCheck.md),
  by simulation.

## Examples

``` r
# \donttest{
f <- ctReport(ctstantestfit, folder = file.path(tempdir(), 'ctReport'),
  components = c('summary', 'identification', 'profile'), profilepoints = 5)
list.files(f$folder)
# }
```
