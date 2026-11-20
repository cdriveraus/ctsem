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

## Examples

``` r
summary(ctstantestfit)
```
