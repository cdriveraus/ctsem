# ctsem

ctsem fits continuous and discrete time dynamic models to longitudinal
data in R: panel data, intensive longitudinal data or a single time
series. Observations can be at any times, and the times can differ
between subjects.

[Website](https://cdriveraus.github.io/ctsem/) · [New in version
4](https://cdriveraus.github.io/ctsem/NEW-IN-4.md) · [Tutorial
paper](https://osf.io/preprints/psyarxiv/4q9ex_v3)

## Features

- Multiple interacting processes, including higher order (such as
  oscillating) and nonlinear dynamics.
- Effects of time-varying inputs, such as interventions.
- Multiple indicators, which can be continuous, binary, ordinal, counts
  or censored.
- Individual differences in any parameter, including at multiple nested
  levels (such as people within couples).
- Covariate effects on any parameter.
- Maximum likelihood, or Bayesian estimation by optimisation or
  sampling.
- Tools for testing a model’s fit and finding how to improve it.
- Latent growth curve, latent change score and cross-lagged panel models
  as special cases.

Version 4 adds an estimation engine written in Julia, now the default.
It is two to three times faster than the Stan backend, uses multiple
cores, and needs no compilation for new models. Random effects are
integrated with a Laplace approximation, which estimates individual
differences in the dynamics more accurately than earlier versions did.

## Installation

``` r

remotes::install_github('cdriveraus/ctsem')
library(ctsem)
ctJuliaInstall()   # once, to set up Julia
```

The first fit after installing takes a few minutes longer while the
engine compiles. Without Julia, `ctFit(..., backend = 'stan')` fits
models as ctsem 3 did.

Compiler problems on Windows

Installing from GitHub compiles the package, which on Windows needs the
Rtools version that matches your R. If compiling fails, check that
`rstan` itself works, and try adding this line to `~/.R/Makevars.win`:

``` R
CXX17FLAGS += -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations
```

The [Stan forums](https://discourse.mc-stan.org/) cover other compiler
problems.

## Example

Milk yield and fat percentage of 53 cows, recorded on most days of a
lactation, from the [doBy](https://cran.r-project.org/package=doBy)
package:

``` r

library(ctsem)
data(milkman_rdm1, package = 'doBy')

milk <- data.frame(
  cow   = milkman_rdm1$cowno,
  day   = milkman_rdm1$dfc,                     # days since calving
  yield = as.numeric(scale(milkman_rdm1$my)),
  fat   = as.numeric(scale(milkman_rdm1$fatpct)))

model <- ctModel(
  type = 'ct',                        # 'dt' for discrete time
  LAMBDA = diag(2),                   # two processes, one indicator each
  manifestNames = c('yield', 'fat'),
  latentNames = c('yield', 'fat'),
  time = 'day',
  id = 'cow')

fit <- ctFit(milk, model)
```

Anything not given to
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
takes a default. The effects between the two processes, for example, are
estimated.

``` r

summary(fit)$popmeans
#>                   mean    sd   2.5%    50%  97.5%   mode
#> T0m_yield       -0.949 0.142 -1.209 -0.950 -0.663 -0.943
#> T0m_fat          0.158 0.185 -0.212  0.155  0.506  0.151
#> drift_yield     -0.041 0.004 -0.049 -0.041 -0.033 -0.041
#> drift_yield_fat -0.020 0.004 -0.027 -0.020 -0.013 -0.020
#> drift_fat_yield -0.006 0.004 -0.013 -0.006  0.002 -0.006
#> drift_fat       -0.015 0.004 -0.024 -0.015 -0.009 -0.015
#> diff_yield       0.136 0.004  0.129  0.136  0.143  0.136
#> diff_fat_yield  -0.591 0.062 -0.708 -0.592 -0.472 -0.592
#> diff_fat         0.116 0.005  0.106  0.116  0.125  0.116
#> mvaryield        0.329 0.003  0.323  0.329  0.335  0.328
#> mvarfat          0.582 0.005  0.573  0.582  0.591  0.582
#> mm_yield        -0.053 0.123 -0.299 -0.052  0.187 -0.054
#> mm_fat           0.141 0.108 -0.063  0.140  0.352  0.142
```

`drift_yield` and `drift_fat` are the auto effects: negative values mean
each process tends back towards its equilibrium. `drift_yield_fat` is
the effect of fat on yield, and `drift_fat_yield` the reverse. Higher
fat is followed by lower yield, while yield has little or no effect on
fat. The `diff` parameters describe the process noise, `mvar` the
measurement error, `mm` the measurement intercepts and `T0m` the initial
states.

``` r

ctPredict(fit, subjects = unique(milk$cow)[1:2], plot = TRUE)
```

![Observed yield and fat percentage for two cows over one lactation,
with the model's expected trajectory and uncertainty band tracking each
series.](reference/figures/README-firstplot-1.png)

Observations for two cows, with the model’s predictions and their
uncertainty.

``` r

ctDiscretePars(fit, plot = TRUE)
```

![Auto and cross regression coefficients plotted against the time
interval. The effect of fat on yield grows steadily more negative; the
reverse effect stays near
zero.](reference/figures/README-firstdiscrete-1.png)

The auto and cross effects as regression coefficients, as they change
with the time interval.

## Covariate effects

This example estimates the effect of diet on the growth of chicks
(`ChickWeight`, in base R):

``` r

chickdata <- ChickWeight
chickdata <- cbind(chickdata, model.matrix(~ Diet - 1, chickdata))  # dummy code diet
colnames(chickdata)[5:8] <- paste0('Diet_', 1:4)
chickdata$weight <- scale(chickdata$weight)
```

``` r

m <- ctModel(
  type = 'ct',
  LAMBDA = diag(1),
  manifestNames = 'weight',
  latentNames = 'Lweight',
  time = 'Time',
  id = 'Chick',
  DRIFT = 'a11, transform=param',                   # allow growth (no sign restriction)
  MANIFESTMEANS = 0,                                # fix the measurement intercept...
  CINT = 'cint, tipreds=c(Diet_2,Diet_3,Diet_4)',   # ...and estimate the process intercept
  TIpredNames = paste0('Diet_', 2:4),               # diet 1 is the baseline
  tipredDefault = FALSE)                            # no other covariate effects

f <- ctFit(chickdata, m, priors = TRUE)
```

In
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
a string is a free parameter and a number a fixed value. Options follow
the parameter’s name, such as `tipreds` for covariate effects or
`indvarying` for random effects;
[`ctPar()`](https://cdriveraus.github.io/ctsem/reference/ctPar.md)
writes these strings.

``` r

summary(f)$tipreds
#>                  mean    sd   2.5%   50% 97.5%  mode     z
#> tip_Diet_2_cint 0.021 0.013 -0.005 0.021 0.046 0.021 1.570
#> tip_Diet_3_cint 0.050 0.013  0.025 0.050 0.074 0.051 3.983
#> tip_Diet_4_cint 0.028 0.013  0.001 0.028 0.052 0.028 2.135
```

Diets 3 and 4 raise the intercept, and so the growth rate, compared with
diet 1.

``` r

ctPredict(f, plot = TRUE, removeObs = TRUE, polygonalpha = 0,
  subjects = as.numeric(chickdata$Chick[!duplicated(ChickWeight$Diet)]))
```

![Predicted weight trajectories for one chick from each of the four
diets, based on the diet covariates alone, separating into four
increasing curves.](reference/figures/README-chickpredict-1.png)

Predicted growth for one chick on each diet, from the diet alone.

## Documentation

- The [website](https://cdriveraus.github.io/ctsem/) has the help pages
  grouped by task, the vignettes and the
  [manual](https://cdriveraus.github.io/ctsem/articles/hierarchicalmanual.pdf).
- A [tutorial paper](https://osf.io/preprints/psyarxiv/4q9ex_v3) works
  through an applied example, with [code and
  materials](https://osf.io/xh6ue/).
- Vignettes: [measurement
  models](https://cdriveraus.github.io/ctsem/articles/measurement-models.html),
  [multilevel random
  effects](https://cdriveraus.github.io/ctsem/articles/laplace-multilevel.html),
  [uncertainty](https://cdriveraus.github.io/ctsem/articles/uncertainty.html),
  [score
  tests](https://cdriveraus.github.io/ctsem/articles/score-tests.html),
  [IRT
  items](https://cdriveraus.github.io/ctsem/articles/irt-items.html),
  [ordinal or
  continuous](https://cdriveraus.github.io/ctsem/articles/ordinal-or-continuous.html),
  [sampling
  weights](https://cdriveraus.github.io/ctsem/articles/sampling-weights.html),
  [empirical
  Bayes](https://cdriveraus.github.io/ctsem/articles/empirical-bayes-fitting.html).
- Blog posts: [quick
  start](https://cdriver.netlify.app/post/ctsem-quick-start/), [dynamic
  systems](https://cdriver.netlify.app/post/introtodynamics/), [discrete
  time models](https://cdriver.netlify.app/post/dtbivariate/),
  [representations of dynamic
  systems](https://cdriver.netlify.app/post/discretetime/), [accelerated
  and multiple group
  designs](https://cdriver.netlify.app/post/accelerated/), [Kalman
  filter and smoother](https://cdriver.netlify.app/post/missingdata/),
  [heterogeneity](https://cdriver.netlify.app/post/heterogeneity/),
  [growth curves](https://cdriver.netlify.app/post/lgc/), [binary
  data](https://cdriver.netlify.app/post/binarydata/).
- Papers on [hierarchical continuous time
  models](https://www.researchgate.net/publication/324093594_Hierarchical_Bayesian_Continuous_Time_Dynamic_Modeling)
  and on [modelling
  interventions](https://www.researchgate.net/publication/328221807_Understanding_the_Time_Course_of_Interventions_with_Continuous_Time_Dynamic_Models).
- [ctsemGUI](https://github.com/cdriveraus/ctsemgui) is a
  point-and-click interface, and
  [ctsemOMX](https://github.com/cdriveraus/ctsemOMX) the original OpenMx
  based ctsem.

## Help

Ask questions in
[discussions](https://github.com/cdriveraus/ctsem/discussions) and
report bugs as [issues](https://github.com/cdriveraus/ctsem/issues). To
cite ctsem, run `citation('ctsem')`.

Older function names

Functions renamed in version 3.11 keep their old names as aliases, so
older code still runs.

| old | new |
|----|----|
| [`ctStanFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) | [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) |
| [`ctStanModel()`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md) | [`ctModelConvertOMX()`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md) |
| [`ctStanGenerate()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromPriors.md) | `ctGenerate(fromPriors = TRUE)` |
| [`ctStanGenerateFromFit()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromFit.md) | [`ctGenerateFromFit()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromFit.md) |
| [`ctStanKalman()`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md) | [`ctKalmanArray()`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md) |
| [`ctStanPlotPost()`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md) | [`ctPlotPosterior()`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md) |
| [`ctStanPostPredict()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredict.md) | [`ctPostPredPlots()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md) |
| [`ctStanSubjectPars()`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md) | [`ctSubjectPars()`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md) |
| [`ctStanTIpredeffects()`](https://cdriveraus.github.io/ctsem/reference/ctTIpredEffects.md) | [`ctTIpredEffects()`](https://cdriveraus.github.io/ctsem/reference/ctTIpredEffects.md) |
| [`ctStanFitUpdate()`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md) | [`ctFitUpdate()`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md) |
| [`ctStanDiscretePars()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md) | [`ctDiscretePars()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md) |
| [`ctStanDiscreteParsPlot()`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md) | [`ctDiscreteParsPlot()`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md) |
| [`ctStanContinuousPars()`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md) | [`ctSummaryMatrices()`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md) |
| [`ctStanParnames()`](https://cdriveraus.github.io/ctsem/reference/ctRawParnames.md) | [`ctRawParnames()`](https://cdriveraus.github.io/ctsem/reference/ctRawParnames.md) |

Renames in version 4 are listed in [New in version
4](https://cdriveraus.github.io/ctsem/NEW-IN-4.md).
