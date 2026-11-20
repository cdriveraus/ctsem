# Likelihood-ratio tests between nested ctsem fits

[`anova()`](https://rdrr.io/r/stats/anova.html) on two or more optimised
fits of nested models to the same data gives each fit's number of
parameters, AIC, BIC, log likelihood and -2 log likelihood, and tests
each fit against the one above it, with the chi-square statistic (twice
the gain in log likelihood), its degrees of freedom (the extra
parameters) and its p-value. The fits may come from either backend.

## Usage

``` r
# S3 method for class 'ctJuliaFit'
anova(object, ...)

# S3 method for class 'ctStanFit'
anova(object, ...)
```

## Arguments

- object:

  An optimised fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- ...:

  Further optimised fits, of models nested with `object` and fitted to
  the same data.

## Value

A data frame of class `anova`, one row per fit in order of parameters,
with columns `npar`, `AIC`, `BIC`, `logLik`, `-2LL`, `Chisq`, `Df` and
`Pr(>Chisq)`; the last three are `NA` on the first row. Rows are named
by the arguments.

## Details

The fits are ordered by their number of parameters, and each is tested
against the one with the next fewer, so the models must be nested in
that order. [`anova()`](https://rdrr.io/r/stats/anova.html) cannot check
this. A negative statistic gives a warning, and means a model fits worse
than one nested in it: either the models are not nested, or an
optimisation stopped short.

When the smaller model fixes a variance at zero, such as a random
effect's, the chi-square p-value is conservative (Stram & Lee, 1994).
[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
tests a new random effect against the correct reference distribution
without refitting.

The log likelihoods are those
[`logLik()`](https://rdrr.io/r/stats/logLik.html) reports. If the random
effects were integrated more exactly after the fit, by its
random-effects check or a refined fit, that integrated value is used,
and the heading says so. A fit with priors maximised the posterior, so
its log likelihood is taken at the posterior mode, and the heading says
that too. Sampled fits have no single log likelihood. Compare them with
[`ctLOO`](https://cdriveraus.github.io/ctsem/reference/ctLOO.md).

## References

Stram, D. O., & Lee, J. W. (1994). Variance components testing in the
longitudinal mixed effects model. *Biometrics*, 50(4), 1171-1177.

## See also

[`ctChisqTest`](https://cdriveraus.github.io/ctsem/reference/ctChisqTest.md)
for the p-value alone.

## Examples

``` r
# \donttest{
    df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)

    m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
    m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)

    if (isTRUE(ctJuliaStatus()$available)) {
      f1 <- ctFit(df,m1,cores=1)
      f2 <- ctFit(df,m2,cores=1)

      test <- anova(f1, f2)
      test
      test$Chisq[2]
    }
# }
```
