# Add the phase predictors used by ctModelPeriodicComponent

Adds to long format data, for each subject, the cosine and sine of the
time of its first row over `period`, as the columns `<name>_cos0` and
`<name>_sin0` that
[`ctModelPeriodicComponent`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md)
uses to set its cycle.

## Usage

``` r
ctPeriodicPhases(datalong, period, name = "cycle", id = "id", time = "time")
```

## Arguments

- datalong:

  Long format data, one row per observation.

- period:

  The cycle's period, in the data's time units.

- name:

  The cycle's name, as given to
  [`ctModelPeriodicComponent`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md).

- id, time:

  The subject and time columns.

## Value

`datalong` with the two columns added.

## Examples

``` r
d <- data.frame(id = rep(1:2, each = 3), time = c(8, 12, 20, 33, 37, 45))
ctPeriodicPhases(d, period = 24)
```
