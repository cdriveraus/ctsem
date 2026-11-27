# Example data for ctsem

Simulated data from 30 subjects, each measured on 10 occasions, used in
the examples and tests. `ctstantestdat` is the same data under its name
before version 4.

## Format

A matrix with 300 rows and 8 columns: the subject `id`, the measurement
`time`, the measured variables `Y1` and `Y2`, a time dependent predictor
`TD1`, and three time independent predictors `TI1`, `TI2` and `TI3`.

## Examples

``` r
# One row per subject and occasion.
head(ctsemtestdat)
```
