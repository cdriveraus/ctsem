# Summarise discrete time effect curves

For each curve returned by
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
one per subject and pair of processes, computes summary quantities for
each draw and their distribution over draws.

## Usage

``` r
# S3 method for class 'ctDiscretePars'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- object:

  Output of
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
  with `plot = FALSE`.

- probs:

  Quantiles of each quantity over draws.

- ...:

  Not used.

## Value

A data frame with one row per subject, pair of processes (`row` the
process affected, `col` the one acting) and quantity, giving the mean,
sd and quantiles over draws. If the object carries the fit's estimate
(`attr(object, 'estimate')`, for an optimised fit), the quantity
computed from the estimate's curve is in the `estimate` column. The
quantities are:

- `peak`:

  The value of largest magnitude, with its sign.

- `peaktime`:

  The time interval at which the peak occurs.

- `halftime`:

  The first interval after the peak at which the curve's magnitude has
  fallen to half the peak's, interpolated between intervals. For an
  autoregression, which peaks at zero, this is its half-life, and with
  `impulseType = 'observed', standardise = TRUE` the half-life of the
  autocorrelation. `NA` for a draw whose curve does not fall that far
  within the intervals computed. `beyond` gives the share of such draws,
  and the summary is over the others.

- `area`:

  The area under the curve over the intervals computed. For a unit
  impulse over long enough intervals, this is the cumulative effect.

The intervals are the `times` given to `ctDiscretePars`, so their
spacing limits the precision of `peaktime` and `halftime`, and their
range limits what `halftime` and `area` can capture.

## Examples

``` r
dp <- ctDiscretePars(ctstantestfit, times = seq(0, 5, .1), nsamples = 20)
summary(dp)
```
