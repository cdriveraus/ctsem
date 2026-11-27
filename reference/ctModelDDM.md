# Drift diffusion model for choice and response time

Builds a ctsem model of choices and response times by the drift
diffusion model, with the drift as a latent process. The model is an
ordinary
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
with one choice and response time indicator (`manifesttype = 5`), so it
is fitted with
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) and can
be changed or extended like any other model, for example with other
indicators of the same latent process.

## Usage

``` r
ctModelDDM(
  name = "rt",
  drift = ~1,
  boundary = ~1,
  ndt = ~1,
  bias = 0.5,
  random = c("drift", "boundary", "ndt"),
  dynamic = FALSE,
  sv = FALSE,
  ...
)
```

## Arguments

- name:

  Name of the data column holding each trial's response time in seconds,
  negative for a response at the lower boundary. With columns `rt` and
  `correct`,
  `dat$rt_signed <- ifelse(dat$correct == 1, dat$rt, -dat$rt)` makes it
  for accuracy coding.

- drift, boundary, ndt:

  One sided formulas of trial level predictors of the drift, the log
  boundary separation and the log nondecision time, such as
  `~ speed + difficulty`. `~ 1` for none. Predictors are numeric columns
  of the data, entered as time dependent predictors. A factor needs
  dummy columns. The drift's predictors shift the drift on their trial
  only.

- bias:

  The relative start point, as a proportion of the boundary separation
  from the lower boundary. A number fixes it, and a one sided formula
  estimates it, with predictors on the logit scale.

- random:

  Which of `'drift'`, `'boundary'`, `'ndt'` and `'bias'` vary over
  subjects.

- dynamic:

  Whether the drift changes over the session. `FALSE` keeps each
  subject's drift constant. `TRUE` makes it a mean reverting
  (Ornstein-Uhlenbeck) process in the data's time column around the
  subject's level, starting from its stationary distribution. The time
  column sets the scale of its rate: trial number, or seconds.

- sv:

  Whether to estimate trial to trial variability of the drift, as a
  standard deviation.

- ...:

  Passed to
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
  such as `id`, `time` or `TIpredNames`.

## Value

A `ctStanModel` object.

## Details

The drift's parameters are named `drift_mean` (its mean, or with
`dynamic = TRUE` `drift_cint`, `drift_auto` and `drift_diffusion` for
the process's intercept, auto effect and diffusion), `boundary`, `ndt`
and `bias`, and each predictor's effect `<parameter>_<predictor>`. The
diffusion coefficient of the decision process is 1.

## See also

[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
whose `WIENER` argument sets the boundary, nondecision time and bias of
a choice and response time indicator directly.

## Examples

``` r
m <- ctModelDDM(name = 'rt_signed', drift = ~ difficulty, boundary = ~ speed,
  dynamic = TRUE)
print(m)
```
