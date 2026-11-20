# Create a data.table to compare data generated from a ctsem fit with the original data.

This function allows for easy comparison of data generated from a fitted
ctsem model with the original data used to fit the model. It provides
options to include residuals in the comparison.

## Usage

``` r
ctPostPredData(fit, residuals = FALSE, nsamples = NA)
```

## Arguments

- fit:

  A fitted ctsem model, from either backend.

- residuals:

  If set to TRUE, includes standardised prior residuals in the
  comparison, as extra rows whose variable names are suffixed
  `' std. res.'`. Needs one Kalman filter pass per generated dataset, so
  it is much slower than the default.

- nsamples:

  Integer or NA. Number of generated datasets to use. If the fit carries
  more than this, draws are subsampled, and if it carries none, this
  many are generated. NA uses the datasets the fit carries, or the
  [`ctGenerateFromFit`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromFit.md)
  default.

## Value

A long data table with one row per generated draw x data row x variable,
carrying the generated `value`, the corresponding `obsValue`, and the
`id`, `Time` and `TimeInterval` of that data row. `LogLik` appears as an
extra 'variable', holding the row-wise log likelihood.

## See also

[`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md),
which turns this into diagnostic plots.

## Examples

``` r
data_comparison <- ctPostPredData(ctstantestfit)
```
