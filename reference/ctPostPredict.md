# Posterior predictive checks for a ctsem fit

`ctPostPredict` and `ctStanPostPredict` are aliases for
[`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md),
except that they print the panels by default.

## Usage

``` r
ctPostPredict(fit, ...)
```

## Arguments

- fit:

  A fitted ctsem model.

- ...:

  Passed to
  [`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md).

## Value

A named list of ggplot objects, invisibly.

## Examples

``` r
# \donttest{
ctPostPredict(ctstantestfit, wait = FALSE)
# }
```
