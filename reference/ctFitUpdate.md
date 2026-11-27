# Update a fit to new data, or to the current version of ctsem

Either to include different data, or because you have upgraded ctsem and
the internal data structure has changed. Works on a fit from either
backend.

## Usage

``` r
ctFitUpdate(oldfit, data = NA, recompile = FALSE, refit = FALSE, ...)
```

## Arguments

- oldfit:

  fit object to be updated, from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
  either backend.

- data:

  replacement long format data object. If not supplied, the data the fit
  was made with.

- recompile:

  whether to force a recompile of the Stan model – safer but slower and
  usually unnecessary. Stan fits only.

- refit:

  if TRUE, refits the model using the old estimates as a starting point.
  Only applicable for optimized fits, not sampling. If FALSE, the fit is
  returned with its data and model specification replaced, and its
  estimates and everything computed from them unchanged.

- ...:

  extra arguments to pass to ctFit

## Value

updated fit object, of the same class as `oldfit`.

## Details

The fit is rebuilt by calling
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) again
with the arguments and model it was made with, and with its own data
unless `data` is supplied. Arguments in `...` replace the fit's own.

With `refit=FALSE` and new `data`, the fit keeps its model and
estimates, including its random effects when the new data have only one
subject. To predict for other data without updating the fit, use
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)`(fit, newdata = data)`.

## Examples

``` r
newfit <- ctFitUpdate(ctstantestfit,refit=FALSE)
```
