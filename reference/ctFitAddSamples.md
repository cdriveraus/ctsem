# Append pseudo-posterior samples to an optimized ctStanFit object (deprecated)

Appends normal draws from the fit's covariance to its raw posterior.
Deprecated, use
[`ctOptimUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
with `uncertainty = 'stored'`, which works on both backends.

## Usage

``` r
ctFitAddSamples(fit, nsamples, cores = 2)
```

## Arguments

- fit:

  fit object

- nsamples:

  number of samples desired

- cores:

  number of cores to use

## Value

fit object with extra samples

## Details

The draws come from the fitted covariance, not from a sampler. To sample
the posterior of a fit, use `ctFitUncertainty(fit, 'sample')`.

`ctOptimUncertainty(fit, uncertainty = 'stored', finishsamples = n)`
draws from the same covariance, records the number of draws in
`$uncertainty$settings`, and replaces the existing draws instead of
appending to them, so that draws from different methods are not mixed.

`ctAddSamples` is the older name of this function, and is also
deprecated.

## See also

[`ctOptimUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md),
[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)

## Examples

``` r
if (FALSE) { # \dontrun{
newfit <- ctOptimUncertainty(ctstantestfit, uncertainty = 'stored',
  finishsamples = 30, cores = 1)
} # }
```
