# ctDiscretePars

Calculate model implied regressions for a sequence of time intervals (if
ct) or steps (if dt) based on a fit from
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), for
specified subjects. Wrap with print() when used inside for loops!
`ctStanDiscretePars` is maintained as a backward-compatible alias.

## Usage

``` r
ctDiscretePars(
  fit,
  subjects = "popmean",
  times = seq(from = 0, to = 10, by = 0.1),
  nsamples = 200,
  impulseType = "unit",
  standardise = FALSE,
  cov = FALSE,
  plot = FALSE,
  cores = getOption("mc.cores", 2L),
  state = NULL,
  method = "linearise",
  realid = TRUE,
  ...,
  ctstanfitobj,
  observational
)
```

## Arguments

- fit:

  model fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)

- subjects:

  Either 'popmean', to use the population mean parameter, or the
  subjects to use, as ids in the data (see `realid`).

- times:

  Numeric vector of positive values, discrete time parameters will be
  calculated for each. If the fit object is a discrete time model, these
  should be positive integers.

- nsamples:

  Number of samples from the stanfit to use for plotting. Higher values
  will increase smoothness / accuracy, at cost of plotting speed. Values
  greater than the total number of samples will be set to total samples.

- impulseType:

  What a "one unit change in process c" brings with it. Each option
  gives `dtDRIFT(t) %*% C` for a different matrix `C`, whose column c
  says what else changes.

  `'unit'` (the default)

  :   `C = I`. Process c alone changes by one unit (one standard
      deviation when `standardise=TRUE`). This is the *partial*
      regression `E[x(t+u)|x(t)]`, and the only option that depends on
      the dynamics alone.

  `'observed'`

  :   `C = Sigma diag(Sigma)^-1` with `Sigma = asymDIFFUSIONcov`.
      Process c is *observed* one unit above its expectation, and the
      other processes change with it by `Sigma_rc/Sigma_cc`. This is the
      *simple* regression, the one obtained by regressing the observed
      data on x_c alone.

  `'shock'`

  :   `C = Q diag(Q)^-1` with `Q = DIFFUSIONcov`. A system noise
      innovation in c, with the correlated innovations `Q_rc/Q_cc` in
      the other processes. This is an impulse response.

  `'orthogonal'`

  :   From the Cholesky factor of DIFFUSION, so shocks are uncorrelated.
      The result depends on the order of the latent processes.

  These matrices are not symmetric: a one unit change in a process with
  a large variance implies a larger change in one with a small variance
  than the reverse. `'observed'` and `'shock'` are the same only when
  every process decays at the same rate and there are no cross effects.

- standardise:

  Logical, `FALSE` by default. If `TRUE`, report the effect of a one
  standard deviation change in process c in standard deviations of
  process r, with both taken from the stationary within subject
  covariance `asymDIFFUSIONcov`. This makes effects between processes on
  different scales comparable, and is what
  [`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
  reports by default. With `FALSE` the effects are in the processes' own
  units, and `impulseType='unit', standardise=FALSE` gives
  `expm(DRIFT*t)`. Draws for a system with no stationary variance are
  returned as NaN, with a message. The plot title states which scale was
  used.

- cov:

  Logical. If TRUE, covariances are returned instead of regression
  coefficients.

- plot:

  Logical. If TRUE, ggplots output using
  [`ctDiscreteParsPlot`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md)
  instead of returning output.

- cores:

  Number of cpu cores to use for computing subject matrices. If subject
  matrices were saved during fiting, not used.

- state:

  Where to evaluate matrix cells that depend on the latent state or a
  time dependent predictor (see
  [`ctContextDependence`](https://cdriveraus.github.io/ctsem/reference/ctContextDependence.md)).
  `NULL` (the default) or `'T0MEANS'` uses the population initial state,
  with time dependent predictors at zero. For a model whose dynamics
  depend on the state this can be far from the data. `'mean'` uses the
  mean smoothed latent state, `'asymptotic'` the system's fixed point,
  and a numeric vector gives the state directly. The initial state and
  the fixed point are computed for each posterior draw, so their
  uncertainty is included in the intervals. `'mean'` uses the estimate's
  mean state for every draw, so its uncertainty is not included. A
  process that does not change (a latent predictor) is held at its
  initial value when finding the fixed point. Requires `backend='julia'`
  and `subjects='popmean'`. For a linear model the argument has no
  effect.

- method:

  How the regressions are computed. `'linearise'` (the default) computes
  `expm(J*t)`, with `J` the Jacobian of the drift (the model's `JAx`) at
  the evaluation point. For a linear model `J` is the DRIFT matrix and
  the result is exact. For a model whose dynamics depend on the state,
  `J` can differ from DRIFT: a DRIFT cell `a + b*Y` in the column of
  process `Y` contributes `a + 2*b*Y` to the Jacobian. A stan fit has no
  Jacobian, and uses DRIFT evaluated at `T0MEANS`. The linearisation is
  accurate near a fixed point and for intervals that are short relative
  to the dynamics. `'simulate'` integrates the nonlinear system from
  `state`, and from `state` plus a shock of one stationary standard
  deviation in each process in turn (with the accompanying changes
  `impulseType` gives), and takes the difference between the
  trajectories, scaled to one unit (one standard deviation when
  `standardise=TRUE`). It gives the same result as `'linearise'` when no
  cell depends on the state. For a strongly nonlinear system the result
  depends on the size of the shock. `'simulate'` requires
  `backend='julia'` and is slower, so it uses `nsamples` posterior
  draws, with a small default.

- realid:

  If TRUE (the default), subjects are read as the ids in the data, and
  integers that are not ids there are read as ctsem's numbering of the
  subjects, 1 to N, with a message. If FALSE, they are read as that
  numbering.

- ...:

  additional plotting arguments to control
  [`ctDiscreteParsPlot`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md)

- ctstanfitobj:

  Deprecated. Use `fit`.

- observational:

  Deprecated. Use `impulseType`: `FALSE` is `'unit'` and `TRUE` is
  `'observed'`.

## Value

With `plot = FALSE`, an array of class `ctDiscretePars` holding the
discrete time coefficients, by sample, subject, time interval, row and
column. [`summary()`](https://rdrr.io/r/base/summary.html)
([`summary.ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/summary.ctDiscretePars.md))
gives each curve's peak, the time of the peak, the time at which it
falls to half the peak (a half-life), and the area under it, with their
uncertainty. For an optimised fit with `subjects = 'popmean'`, the
coefficients at the fit's estimate are in `attr(x, 'estimate')`, by
subject, time interval, row and column. These are computed from the
matrices `ctSummaryMatrices(fit, state = )$estimate` reports, or for
`'simulate'` simulated at the estimate, and can differ from the mean or
median of the draws. [`summary()`](https://rdrr.io/r/base/summary.html)
reports them in an `estimate` column, and
[`ctDiscreteParsPlot`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md)
draws them as its line. Sampled fits, and subjects other than
`'popmean'`, have no such estimate. With `plot = TRUE`, a ggplot object.

## Details

`impulseType='unit'` with `standardise=FALSE` gives the discrete time
autoregressive and cross-lagged matrix, `expm(DRIFT*t)`.
`impulseType='observed'` with `standardise=TRUE` gives the model implied
*cross-correlation function*, `Cor(x_r(t+u), x_c(t))`, which corresponds
to the empirical autocorrelations
[`ctACF`](https://cdriveraus.github.io/ctsem/reference/ctACF.md)
computes from the data, and equals the latent correlation matrix at
`t=0`.

If plot=TRUE, the function will return a ggplot2 object (and hence needs
to be printed if intended to display within a loop). This can be
modified by the various ggplot2 functions, or displayed using print(x).

## Examples

``` r
data.table::setDTthreads(1) #ignore this line
# The auto and cross effects of the bundled example fit as regression
# coefficients, for time intervals from 0 to 4 units. Choose the range to
# cover the time scale of the processes.
ctDiscretePars(ctsemtestfit, times = seq(0, 4, .1), plot = TRUE)

# Without plot, an array of coefficients. summary() gives each curve's
# peak, the interval at which it peaks, and the interval at which it has
# fallen to half the peak.
dp <- ctDiscretePars(ctsemtestfit, times = seq(0, 4, .1))
summary(dp)

# Cross effects only, standardised by the stationary standard deviations.
# The plot is a ggplot and can be modified.
g <- ctDiscretePars(ctsemtestfit, times = seq(0, 4, .1),
  plot = TRUE, indices = 'CR', standardise = TRUE)
g + ggplot2::labs(title = 'Cross effects')
```
