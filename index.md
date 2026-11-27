# ctsem

ctsem fits continuous and discrete time dynamic models to longitudinal
data in R: panel data, intensive longitudinal data or a single time
series. Observations can be at any times, and the times can differ
between subjects.

[Website](https://cdriveraus.github.io/ctsem/) · [Getting
started](https://cdriveraus.github.io/ctsem/articles/ctsem.html) · [New
in version 4](https://cdriveraus.github.io/ctsem/NEW-IN-4.md) ·
[Tutorial paper](https://osf.io/preprints/psyarxiv/4q9ex_v3)

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
package. ctsem takes data in long format: one row per cow and occasion,
with columns for the id, the time and each measured variable. Missing
values are `NA`, and the times can be irregular and differ between
subjects.

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
takes a default. Here the two processes can affect each other, and each
cow has its own average levels and starting values.

``` r

summary(fit)$popmeans
#>                   mean    sd   2.5%    50%  97.5%   mode
#> T0m_yield       -0.940 0.148 -1.227 -0.938 -0.654 -0.943
#> T0m_fat          0.160 0.191 -0.206  0.158  0.527  0.150
#> drift_yield     -0.041 0.004 -0.050 -0.040 -0.033 -0.041
#> drift_yield_fat -0.020 0.004 -0.028 -0.020 -0.013 -0.020
#> drift_fat_yield -0.006 0.004 -0.013 -0.006  0.001 -0.006
#> drift_fat       -0.015 0.004 -0.025 -0.014 -0.008 -0.015
#> diff_yield       0.136 0.004  0.128  0.136  0.143  0.136
#> diff_fat_yield  -0.593 0.059 -0.714 -0.592 -0.481 -0.592
#> diff_fat         0.116 0.005  0.107  0.116  0.126  0.116
#> mvaryield        0.329 0.003  0.323  0.329  0.335  0.328
#> mvarfat          0.582 0.005  0.573  0.582  0.591  0.582
#> mm_yield        -0.051 0.122 -0.286 -0.048  0.184 -0.054
#> mm_fat           0.141 0.105 -0.049  0.138  0.358  0.142
```

The `drift` parameters are rates of change per day, the time unit of the
data. `drift_yield` and `drift_fat` are the auto effects, and
`drift_yield_fat` is the effect of fat on the rate of change of yield,
`drift_fat_yield` the reverse. The `diff` parameters describe the random
shocks to the processes, `mvar` the measurement error, `mm` the average
levels of the measurements and `T0m` the starting values.
`diff_fat_yield` is a parameter for the correlation of the shocks, which
`summary(fit)$parmatrices` reports as `DIFFUSIONcor`.
`summary(fit)$popsd` gives how much the cows differ.

Whether the processes return to equilibrium after a disturbance, and how
quickly, depends on the whole drift matrix. They do when every
eigenvalue of the matrix has a negative real part:

``` r

drift <- ctSummaryMatrices(fit)$estimate$DRIFT
eigen(drift)$values
#> [1] -0.04457955 -0.01083888
log(2) / abs(eigen(drift)$values)    # half-lives in days
#> [1] 15.54855 63.95007
```

A disturbance halves in about 16 days along one combination of the
processes and 64 days along the other. That is slower than either auto
effect alone would give, because the cross effects form a loop: higher
fat lowers yield, and lower yield raises fat.

``` r

ctDiscretePars(fit, times = seq(0, 150, by = 1), plot = TRUE)
```

![Auto and cross regression coefficients plotted against time intervals
of up to 150 days. The autoregressions fall towards zero, and the effect
of fat on yield becomes negative, is strongest at about six weeks, and
then returns towards zero; the effect of yield on fat stays near
zero.](reference/figures/README-firstdiscrete-1.png)

The same effects as regression coefficients, as they change with the
time interval. A higher fat percentage is followed by a lower yield,
most strongly after about six weeks, and yield has little effect on fat.
The [dynamics
vignette](https://cdriveraus.github.io/ctsem/articles/dynamics.html)
explains how to read these.

``` r

ctPredict(fit, subjects = unique(milk$cow)[1:2], plot = TRUE)
```

![Observed yield and fat percentage for two cows over one lactation,
with the model's expected trajectory and uncertainty band tracking each
series.](reference/figures/README-firstplot-1.png)

Observations for two cows, with the model’s predictions and their
uncertainty.

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

The same models can be written in lavaan, Mplus or brms style syntax
with
[`ctModelLavaan()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
[`ctModelMplus()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
and
[`ctModelFormula()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
see the [model syntax
vignette](https://cdriveraus.github.io/ctsem/articles/model-syntax.html).

``` r

summary(f)$popmeans[c('a11', 'cint'), ]
#>       mean    sd  2.5%   50% 97.5%  mode
#> a11  0.032 0.003 0.026 0.032 0.039 0.032
#> cint 0.095 0.008 0.079 0.095 0.111 0.095
summary(f)$tipreds
#>                  mean    sd   2.5%   50% 97.5%  mode     z
#> tip_Diet_2_cint 0.022 0.013 -0.004 0.021 0.049 0.021 1.616
#> tip_Diet_3_cint 0.051 0.013  0.026 0.051 0.077 0.051 3.844
#> tip_Diet_4_cint 0.029 0.013  0.003 0.029 0.054 0.028 2.225
```

The rate at which a chick’s weight changes is `a11` times its current
weight, plus `cint`. `a11` is positive, so the heavier a chick, the
faster it grows. `cint` is the growth rate for diet 1 at the average
weight. It varies between chicks, as a free intercept does by default,
and the diets shift it. Diets 3 and 4 raise the growth rate at any
weight, compared with diet 1, and because heavier chicks grow faster,
the difference compounds over time. The effect of diet 2 is smaller and
its interval includes zero.

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
- Vignettes: [getting
  started](https://cdriveraus.github.io/ctsem/articles/ctsem.html),
  [reading continuous time
  dynamics](https://cdriveraus.github.io/ctsem/articles/dynamics.html),
  [interventions and time dependent
  predictors](https://cdriveraus.github.io/ctsem/articles/interventions.html),
  [model
  syntax](https://cdriveraus.github.io/ctsem/articles/model-syntax.html),
  [measurement
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
  weights](https://cdriveraus.github.io/ctsem/articles/sampling-weights.html).
- Examples on the website: [why continuous
  time](https://cdriveraus.github.io/ctsem/articles/continuous-time.html),
  [a single time
  series](https://cdriveraus.github.io/ctsem/articles/sunspots.html),
  [Kalman filter and
  smoother](https://cdriveraus.github.io/ctsem/articles/kalman-filter.html),
  [latent growth
  curves](https://cdriveraus.github.io/ctsem/articles/growth-curves.html),
  [discrete time
  models](https://cdriveraus.github.io/ctsem/articles/discrete-time.html),
  [accelerated and multiple group
  designs](https://cdriveraus.github.io/ctsem/articles/accelerated-designs.html),
  [binary and continuous
  indicators](https://cdriveraus.github.io/ctsem/articles/binary-data.html),
  [individual differences that come from the
  model](https://cdriveraus.github.io/ctsem/articles/heterogeneity.html).
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
