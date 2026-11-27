# Example fit for ctsem

A fit of
[`ctsemtestdat`](https://cdriveraus.github.io/ctsem/reference/ctsemtestdat.md)
with the default backend, julia, used in the examples.
[`summary()`](https://rdrr.io/r/base/summary.html),
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md),
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
[`ctExtract`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md),
[`ctModelLatex`](https://cdriveraus.github.io/ctsem/reference/ctModelLatex.md)
and
[`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
work on it without Julia. Functions that run the Kalman filter, generate
data or fit again need Julia.

## Format

A fit from
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), made
by

    model <- ctModel(LAMBDA = diag(2), manifestNames = c('Y1', 'Y2'),
      latentNames = c('eta1', 'eta2'), TDpredNames = 'TD1')
    ctsemtestfit <- ctFit(ctsemtestdat, model, priors = TRUE)

## See also

[`ctstantestfit`](https://cdriveraus.github.io/ctsem/reference/ctstantestfit.md),
a fit of the same data with `backend = 'stan'`.

## Examples

``` r
summary(ctsemtestfit)
```
