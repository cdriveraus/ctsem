# Plot an approximate continuous-time ACF object from ctACF

Plot an approximate continuous-time ACF object from ctACF

## Usage

``` r
plotctACF(
  ctacfobj,
  df = "auto",
  quantiles = c(0.025, 0.5, 0.975),
  separateLearnRates = FALSE,
  reducedXlim = 1,
  estimateSpline = TRUE,
  method = c("weighted", "quantile")
)
```

## Arguments

- ctacfobj:

  object

- df:

  basis dimension for the spline. 'auto' uses up to 40 for
  `method='weighted'`, and 10% of the number of time intervals for
  `method='quantile'`.

- quantiles:

  quantiles to plot.

- separateLearnRates:

  if TRUE, estimate the learning rate for the quantile splines for each
  combination of variables. Slower but theoretically more accurate. Only
  used with `method='quantile'`.

- reducedXlim:

  if non-zero, n timesteps are removed from the upper and lower end of
  the x range where the spline estimates are less likely to be
  reasonable.

- estimateSpline:

  if TRUE, a spline is estimated according to `method`, otherwise the
  samples are simply plotted as lines and the other arguments here are
  not used.

- method:

  `'weighted'` (default) fits a spline through the full data estimates,
  weighting each time interval by the inverse of its bootstrap variance,
  so that precisely estimated features such as a sharp peak are followed
  while noisy estimates are smoothed. Needs `nboot > 1` in `ctACF`.
  `'quantile'` fits quantile regression splines to the pooled bootstrap
  samples, which smooths all time intervals alike and can flatten sharp
  peaks.

## Value

a ggplot object

## Details

With a spline, the full data estimate at each time interval is shown as
a point, the spline as a line, and the band spans the outer quantiles.

## Examples

``` r
data.table::setDTthreads(1) #ignore this line
# Example usage:
head(ctstantestdat)
ac=ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',timestep=.5,nboot=5,plot=FALSE)
plotctACF(ac, reducedXlim=0)
```
