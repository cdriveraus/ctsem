# Chi Square test wrapper for ctsem fit objects.

The p-value of a likelihood-ratio test between two nested model fits,
from either backend.
[`anova`](https://cdriveraus.github.io/ctsem/reference/anova.ctJuliaFit.md)`(fit1, fit2)`
gives the whole test: each fit's log likelihood, the chi-square
statistic, its degrees of freedom and the p-value.

## Usage

``` r
ctChisqTest(fit1, fit2)
```

## Arguments

- fit1:

  One of the fits to be compared.

- fit2:

  Second fit to be compared. The fit with more parameters is tested
  against the other, whichever order they come in.

## Value

The p-value: that of the fit with more parameters in
`anova(fit1, fit2)`. `NA` where the two have the same number of
parameters.

## Examples

``` r
# \donttest{
    df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)

    m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
    m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)

    if (isTRUE(ctJuliaStatus()$available)) {
      f1 <- ctFit(df,m1,cores=1)
      f2 <- ctFit(df,m2,cores=1)

      ctChisqTest(f1,f2)
    }
# }
```
