# ctSummaryMatrices

Summarise model-implied parameter matrices from a ctsem fit object.
`ctStanContinuousPars` is maintained as a backward-compatible alias.

## Usage

``` r
ctSummaryMatrices(
  fit,
  calcfunc = quantile,
  calcfuncargs = list(probs = 0.5),
  timeinterval = 1,
  state = NULL,
  ...
)
```

## Arguments

- fit:

  fit object from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)

- calcfunc:

  Function to apply over samples, must return a single value. By default
  the median over all samples is returned using the
  [`quantile`](https://rdrr.io/r/stats/quantile.html) function, but one
  might also be interested in the
  [`mean`](https://rdrr.io/r/base/mean.html) or
  [`sd`](https://rdrr.io/r/stats/sd.html), for instance.

- calcfuncargs:

  A list of additional parameters to pass to calcfunc. For instance,
  with the default of calcfunc = quantile, the probs argument is needed
  to ensure only a single value is returned.

- timeinterval:

  time interval for discrete time parameter matrix computation.

- state:

  The latent state at which state dependent matrix cells are evaluated:
  `'T0MEANS'` (the default), `'mean'`, `'asymptotic'`, or a numeric
  state vector. The default is the population initial state, with time
  dependent predictors at zero. For a model whose dynamics depend on the
  state this can be far from the data, and `'mean'` (the mean smoothed
  latent state) or `'asymptotic'` (the system's fixed point) is usually
  more useful. The initial state and the fixed point are computed for
  each draw. `'mean'` uses the estimate's mean state for every draw, so
  its uncertainty is not included. Values other than the default need
  `backend='julia'` and are an error for a stan fit. No effect for a
  linear model.

- ...:

  arguments passed to the method for the fit's backend.

## Value

A named list of matrices, each summarised over the posterior (or
approximate posterior) draws by `calcfunc`. These are the model's own
matrices (`DRIFT`, `DIFFUSION`, `LAMBDA`, ...), the covariances computed
from them (`DIFFUSIONcov`, `MANIFESTcov`, and `T0cov`, the initial
latent covariance over subjects, which is `T0VAR` plus the population
covariance of any individually varying `T0MEANS`), the asymptotic
`asymCINT` and `asymDIFFUSIONcov`, and `dtDRIFT` at `timeinterval`. Each
covariance has a matching correlation matrix, `DIFFUSIONcor`,
`MANIFESTcor`, `T0cor` and `asymDIFFUSIONcor`, computed for each draw
before `calcfunc` is applied. A variable with zero variance has `NA`
correlations.

On `backend='julia'` the list also holds the Jacobians of the drift and
of the measurement at the evaluation point, `JAx` (latent by latent) and
`Jy` (manifest by latent). These are the local effects of the model
linearised at `state`. They equal `DRIFT` and `LAMBDA` for a linear
model, and differ where a cell is nonlinear in the state: a `DRIFT` cell
`a + q*Y` in the column of process `Y` contributes `a + 2*q*Y` to `JAx`.
`dtDRIFT` is `expm(JAx * timeinterval)` where there is a `JAx`, as in
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
and `expm(DRIFT * timeinterval)` on stan. For a model whose dynamics
depend on the state this is accurate near a fixed point or over a short
interval. Elsewhere, compare with `ctDiscretePars(method = 'simulate')`.

Where parameters vary over subjects, the list also holds their
population covariance on the natural scale, `popcov`, and their
raw-scale correlations, `rawpopcorr`, named by parameter (with several
grouping levels, `popcov.<level>` and `rawpopcorr.<level>`). `popcov`
combines the sds [`summary()`](https://rdrr.io/r/base/summary.html)
reports as `popsd` with those correlations, so it is exact where a
parameter's transform is affine, as for intercepts and initial means.

On an optimised fit, `estimate` holds the same matrices at the fit's
estimate, unrounded. For example `ctSummaryMatrices(fit)$estimate$DRIFT`
is the estimated drift matrix, which can differ from the median of the
draws, and `ctSummaryMatrices(fit, state = s)$estimate$JAx` is the model
linearised at `s` at the estimate. A sampled fit has no `estimate`, and
a `backend='stan'` fit's has no `popcov`.

## Examples

``` r
# \donttest{
#posterior median over all subjects (also reflects mean of unconstrained pars)
ctSummaryMatrices(ctstantestfit)
# }
```
