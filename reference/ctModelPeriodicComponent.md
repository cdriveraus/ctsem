# Add a cycle, such as a daily rhythm, that drives a model's latent processes

Adds a pair of latents, `<name>_cos` and `<name>_sin`, that rotate with
a fixed `period` and have no noise of their own
(`DRIFT[cos, sin] = -2 pi / period` and
`DRIFT[sin, cos] = 2 pi / period`), so they follow the cosine and sine
of the time of day, or of whatever the period is. Each latent in
`indices` is driven by both, through estimated effects
`DRIFT[k, <name>_cos]` and `DRIFT[k, <name>_sin]`, which set the
amplitude and phase of the cycle's effect. Each driven process responds
with its own dynamics, so a fast process follows the cycle closely and a
slow one smooths and delays it.

## Usage

``` r
ctModelPeriodicComponent(
  ctm,
  period,
  indices = seq_len(ctm$n.latent),
  name = "cycle",
  individual = TRUE
)
```

## Arguments

- ctm:

  A model from `ctModel(type = 'ct')`.

- period:

  The cycle's period, in the data's time units: 24 for a daily rhythm
  with time in hours.

- indices:

  Which latents the cycle drives, by position.

- name:

  Prefix for the added latents, predictors and parameters.

- individual:

  Let each subject's cycle deviate from the clock.

## Value

The extended model.

## Details

The cycle is set to the clock at each subject's first observation, using
two time independent predictors, `<name>_cos0` and `<name>_sin0`, the
cosine and sine of that time over the period, which
[`ctPeriodicPhases`](https://cdriveraus.github.io/ctsem/reference/ctPeriodicPhases.md)
adds to the data. The time variable must therefore count from an origin
shared by all subjects, such as hours since midnight of a common day.
With `individual = TRUE` each subject's cycle may also be shifted from
the clock and differ in strength, through an initial deviation of the
pair with an estimated standard deviation. This allows rhythms that
differ between people, or a time origin that is not shared.

The effects start at zero.

Needs the julia backend.

## See also

[`ctPeriodicPhases`](https://cdriveraus.github.io/ctsem/reference/ctPeriodicPhases.md)
for the predictors it uses,
[`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)'s
`periodic` family, which detects cycles a model is missing, and
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md),
which adds a cycle where that family finds one.

## Examples

``` r
m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'), LAMBDA = diag(2))
m2 <- ctModelPeriodicComponent(m, period = 24)
m2$latentNames
m2$TIpredNames
p <- m2$pars
p[p$matrix %in% 'DRIFT' & (p$row > 2 | p$col > 2) &
  (!is.na(p$param) | p$value != 0), c('matrix', 'row', 'col', 'param', 'value')]
```
