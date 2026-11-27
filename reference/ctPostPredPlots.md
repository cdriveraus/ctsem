# Posterior predictive checks for a ctsem fit

Compares data generated from the fitted model against the data it was
fitted to, as a set of diagnostic plots. Works for both the stan and
julia backends. `ctPostPredict` and `ctStanPostPredict` are maintained
as aliases.

## Usage

``` r
ctPostPredPlots(
  fit,
  panels = "all",
  variables = NULL,
  nsamples = NA,
  datarows = "all",
  diffsize = 1,
  interval = 0.95,
  resolution = 100,
  jitter = 0,
  residuals = FALSE,
  timebins = 10,
  lags = 0:5,
  lagcor = TRUE,
  notes = TRUE,
  plot = FALSE,
  wait = FALSE
)
```

## Arguments

- fit:

  A fitted ctsem model.

- panels:

  Which panels to produce. `'all'` (the default), `'structure'`,
  `'calibration'`, or a character vector of panel names (see Details).

- variables:

  Character vector of variable names to restrict to, or NULL for all.
  `'LogLik'` is a valid name here.

- nsamples:

  Integer or NA. Number of generated datasets to compare against. NA
  uses the datasets the fit carries.

- datarows:

  Integer vector of data rows to use, or `'all'`.

- diffsize:

  Integer vector \> 0. Lags, in observations, for the rate-of-change
  panels. One pair of panels is produced per lag.

- interval:

  Numeric in (0,1). Width of the central predictive interval used by the
  `PredictedVsObserved` and `IntervalCoverage` panels.

- resolution:

  Positive integer. Grid resolution of the shaded model density in the
  bivariate panels.

- jitter:

  Non-negative numeric. Jitter added to the values in the rate-of-change
  panels, as a proportion of each variable's standard deviation, applied
  to observed and generated alike. Useful when data are rounded or
  ordinal. 0, the default, adds none.

- residuals:

  Logical. Include standardised prior residuals as extra variables.
  Slow, as it needs one Kalman filter pass per generated dataset.

- timebins:

  Positive integer. Number of equal-count bins of time for the
  `MeanTrajectory` panel, which shows the mean within each bin.

- lags:

  Non-negative integer vector. Observation lags for the
  `LaggedCovariance` panel, passed to
  [`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md).

- lagcor:

  Logical. If TRUE (the default) the `LaggedCovariance` panel shows
  correlations, which are usually easier to read; FALSE shows
  covariances.

- notes:

  Logical. If TRUE (the default), each panel carries a caption saying
  what it shows and what a departure from its reference line means. Set
  FALSE for bare plots.

- plot:

  Logical. If TRUE, prints the panels and returns them invisibly. If
  FALSE (the default) returns the list without printing.

- wait:

  Logical. If TRUE and `plot=TRUE`, waits for a keypress between panels.

## Value

A named list of ggplot objects, invisibly if `plot=TRUE`.

## Details

Two families of panel are produced.

**structure** panels check whether the model reproduces the shape of the
data: `Density` (marginal distributions), `MeanTrajectory` (the mean
over subjects against time), `ValueByTime` and `ValueByOccasion` (the
model implied density of each variable against time and against occasion
number, with the observed data as points), `ChangeByValue` and
`ChangeByTime` (rate of change against current value and against time),
and `LaggedCovariance` (the lagged covariance or correlation structure).
`ChangeByValue` is the phase plane. A linear model implies a straight
band with a constant slope there, so curvature in the observed points
relative to the model shading suggests nonlinear dynamics.

`ValueByTime` shows where the model expects individual observations, so
its band is as wide as the data. `MeanTrajectory` shows where it expects
the mean of a bin of observations, which has a much narrower band, so a
bias in the average level is easier to see there.

`LaggedCovariance` shows
[`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
as a panel, with one plot per variable, faceted by the variable it is
paired with, comparing the empirical values to the central 95% of the
generated datasets at each lag. It is the only panel that looks at
dependence across time, and so the one most likely to show a wrong drift
or diffusion when the marginal distributions match. Call
[`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
directly to split it by a group or a predictor.

**calibration** panels check whether the predictive distribution has the
right width and location: `PredictedVsObserved`, `IntervalCoverage`,
`PIT`, and the PIT split by `CalibrationByInterval` and
`CalibrationBySubject`. A problem that depends on the time interval
points to the dynamics, since the measurement model does not depend on
the interval.

The row-wise log likelihood appears in most panels as an extra variable,
`LogLik`, whose distribution should match between model and data like
any other. It is left out of the rate-of-change panels. The observed log
likelihood is evaluated at the point estimate, so it has no parameter
uncertainty, while the generated one does.

`SubjectLogLik` compares each subject's total log likelihood, as a
cumulative distribution across subjects, against the same curve from
each generated dataset. This shows subjects the model fits badly over
many occasions, which the row-wise panels can miss. The reference is the
generated datasets of the same size, so it is valid for any number of
subjects.

The comparison is only as precise as the number of generated datasets
allows. With few draws the picture is coarse, and the PIT takes only
`nsamples + 1` distinct values.

## See also

[`ctPostPredData`](https://cdriveraus.github.io/ctsem/reference/ctPostPredData.md)
for the underlying table,
[`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
for the lagged covariance panel on its own, with splits by group or
predictor, and
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
for the model-implied phase plane on its own.

## Examples

``` r
# \donttest{
plots <- ctPostPredPlots(ctstantestfit)
print(plots$PIT)
# }
```
