# Fit again from the other end of a flat parameter's range

A fit can leave a parameter far out where the likelihood hardly changes
with it, such as an auto-effect so fast or so slow that its value no
longer matters, or a variance or random effect standard deviation at
zero. A higher maximum may then have that parameter elsewhere and the
rest of the model arranged around it.
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) names
such parameters in a message. This function fits the model again from
the default starting values with each of them held at the other end of
its range, then lets it go, and keeps the better fit.

## Usage

``` r
ctFitFlatRetry(fit, cores = NULL, ...)
```

## Arguments

- fit:

  A julia fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- cores:

  Number of CPU cores to use. NULL uses the fit's own setting.

- ...:

  Other arguments passed to
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), as
  for
  [`ctFitUpdate`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md).

## Value

The new fit where fitting again raised the log likelihood, otherwise
`fit` unchanged.

## Details

Each parameter costs up to two more optimisations.
`optimcontrol = list(flatretry = TRUE)` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) does
the same during the fit.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- ctFit(datalong, model)
fit <- ctFitFlatRetry(fit)
} # }
```
