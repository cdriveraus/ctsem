# Continuous Time Autocorrelation Function (ctACF)

Computes autocorrelations and cross correlations by time lag, pooled
over subjects, for data observed at irregular times.

## Usage

``` r
ctACF(
  dat,
  varnames = "auto",
  ccfnames = "all",
  idcol = "id",
  timecol = "time",
  plot = TRUE,
  timestep = "auto",
  time.max = "auto",
  nboot = 100,
  scale = FALSE,
  center = FALSE
)
```

## Arguments

- dat:

  The input data in data frame or data table format, in long format.

- varnames:

  Character vector of variable names in the data to compute the ACF for.
  'auto' uses all columns that are not time / id.

- ccfnames:

  Character vector of variable names in the data to compute cross
  correlation for. 'all' uses all variables in varnames, NA uses none.

- idcol:

  The name of the column containing subject IDs (default is 'id').

- timecol:

  The name of the column containing time values (default is 'time').

- plot:

  A logical value indicating whether to create a plot (default is TRUE).

- timestep:

  The width of the time lag bins. 'auto' (default) uses half the median
  time interval between consecutive rows of a subject.

- time.max:

  The largest time lag to compute. 'auto' (default) uses two thirds of
  the mean, over rows, of the time span of the row's subject.

- nboot:

  The number of bootstrap samples, each resampling subjects with
  replacement (default is 100).

- scale:

  if TRUE, divide each variable by its within-subject standard
  deviation.

- center:

  if TRUE, subtract each variable's within-subject mean.

## Value

If 'plot' is TRUE, the function returns a ggplot object of the ACF plot.
If 'plot' is FALSE, it returns a data table with a row for each sample,
variable and time lag: `Sample` (0 for the data, then the bootstrap
samples), `ACF`, `TimeInterval`, `Variable`, `SignificanceLevel` (the
95% bound of a correlation of zero computed from that many independent
pairs), and `Pairs`, the number of pairs the correlation is computed
from.

## Details

Within each subject, every pair of observations is put in the time lag
bin nearest to the time between them, at multiples of `timestep` up to
`time.max`. The value at a lag is the correlation of the pairs in its
bin, pooled over subjects. For an autocorrelation this is a variable at
time t + lag with the same variable at time t. For the cross correlation
`'x_y'` it is x at time t + lag with y at time t, so at a positive lag y
is the earlier. Lags with fewer than three pairs are left out.
Autocorrelations are returned from the first non-zero lag, as the lag 0
autocorrelation is always 1. Cross correlations include lag 0, the
contemporaneous correlation. Bootstrap samples resample whole subjects,
and the plot shows a spline through the estimates with a band from their
bootstrap spread (see
[`plotctACF`](https://cdriveraus.github.io/ctsem/reference/plotctACF.md)).

## Examples

``` r
data.table::setDTthreads(1) #ignore this line
# Example usage:
head(ctstantestdat)
ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',nboot=5)
```
