# Summarise a ctsem fit

Summarise a fit from
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
either backend: population parameters with their uncertainty,
random-effect spreads and correlations, time independent predictor
effects, the parameter matrices, and the log likelihood and AIC.

## Usage

``` r
# S3 method for class 'ctJuliaFit'
summary(
  object,
  timeinterval = 1,
  digits = 3,
  parmatrices = TRUE,
  priorcheck = TRUE,
  residualcov = TRUE,
  ...
)

# S3 method for class 'ctStanFit'
summary(
  object,
  timeinterval = 1,
  digits = 3,
  parmatrices = TRUE,
  priorcheck = TRUE,
  residualcov = TRUE,
  ...
)
```

## Arguments

- object:

  fit object from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), of
  class `ctStanFit` or `ctJuliaFit`.

- timeinterval:

  positive numeric indicating time interval to use for discrete time
  parameter calculations reported in summary.

- digits:

  integer denoting number of digits to report.

- parmatrices:

  if TRUE, also return additional parameter matrices – can be slow to
  compute for large models with many samples.

- priorcheck:

  Whether or not to use `priorchecking` to compare posterior mean and sd
  to prior mean and sd. `backend='stan'` only. `priorcheck = TRUE` is an
  error for a julia fit.

- residualcov:

  Whether or not to show standardised residual covariance. Takes a
  little longer to compute.

- ...:

  Additional arguments to pass to `priorcheckreport`, such as `meanlim`,
  or `sdlim`.

## Value

List containing summary items, with a `print` method for readable
console and knitr output. On an optimised `backend = 'julia'` fit with
uncertainty draws, the parameter tables and `parmatrices` have a `mode`
column beside the draws' mean, sd and quantiles. This is the fit's
estimate on the same scale, the point that likelihood ratio tests and
information criteria refer to.

`popmeans` reports each free parameter. A `T0VAR`, `DIFFUSION` or
`MANIFESTVAR` parameter on the diagonal is a standard deviation. Off the
diagonal it is a parameter that the model's `covmattransform` turns into
correlations. Under the default `'z'` it is an entry of the logarithm of
the correlation matrix. For two variables this is the correlation's
Fisher z ([`tanh()`](https://rdrr.io/r/base/Hyperbolic.html) of it is
the correlation), and for more variables it does not correspond to any
single correlation. The correlations are in `parmatrices`, as `T0cor`,
`DIFFUSIONcor`, `MANIFESTcor` and `asymDIFFUSIONcor` beside the
covariances, each summarised over the draws' own correlations.

Under the default `'z'`, a covariance parameter fixed to 0 gives a
covariance of exactly 0. An estimated parameter that is 0 or near 0 does
not show that the covariance is 0 when there are three or more
variables. To see whether a covariance or correlation is zero, or how
large it is, read `T0cov`, `DIFFUSIONcov`, `MANIFESTcov`,
`asymDIFFUSIONcov` and the correlation matrices.

## Examples

``` r
s <- summary(ctsemtestfit)

# Printing shows every section, some with a note on how to read them.
s

# The population means of the free parameters, and the standard deviations
# over subjects of those that vary, each with its uncertainty.
s$popmeans
s$popsd

# The model matrices. dtDRIFT holds the autoregressive and cross-lagged
# effects over one unit of time (the timeinterval argument).
s$parmatrices[s$parmatrices$matrix %in% c('DRIFT', 'dtDRIFT'), ]

s$loglik
s$aic
```
