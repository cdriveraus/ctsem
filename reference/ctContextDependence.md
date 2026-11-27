# Context dependence of reported model matrices

A ctsem matrix cell may be an expression of a latent process or a time
dependent predictor, such as `DRIFT[1,1] = '-log1p(exp(param)) * eta2'`.
Such a cell has a value only at a given evaluation context: a latent
state, a set of time dependent predictor values, a time and an interval.

## Details

ctsem's summaries and plots report a number for every cell, so for such
a model they evaluate it at one of two points, depending on the
function:

- **Population matrices** (`summary`, `ctSummaryMatrices`,
  `ctDiscretePars(subjects='popmean')`, `ctTIpredEffects`) are evaluated
  at the population `T0MEANS` state, with time dependent predictors at
  zero.

- **Subject matrices** (`ctSubjectPars`, `ctDiscretePars(subjects=)`)
  are evaluated at each subject's last observed row.

[`ctKalman`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
and
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
evaluate every cell at every step of the filter, so they report no
single point.

References to an *individually varying* parameter are not affected,
since the subject's estimate of such a parameter uses all of that
subject's data. Only references to dynamic processes and to time
dependent predictors make a reported value depend on the evaluation
point.

Functions affected by this attach the cells as
`attr(x, 'contextDependent')`, and give a message naming the point.
[`ctModelIsNonlinear`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md)
reports whether a model has such cells. With `backend='julia'` the
evaluation point can be chosen with the `state` argument of
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md),
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
[`ctExtract`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md)
and
[`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md):
`'T0MEANS'`, `'mean'` (the mean smoothed latent state), `'asymptotic'`
(the fixed point the processes settle at), or a numeric vector with one
value per latent process.
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
and
[`ctStateDependencePlot`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)
show how the dynamics change over the state.

## See also

[`ctModelIsNonlinear`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md),
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md),
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md),
[`ctStateDependencePlot`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
# The auto effect of eta1 changes with the level of eta2, so DRIFT[1, 1]
# has a value only at a given state of eta2.
model <- ctModel(LAMBDA = diag(2),
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
  PARS = c('a|-log1p_exp(param)', 'b'),
  DRIFT = matrix(c('a + b * eta2', 'drift21', 0, 'drift22'), 2, 2),
  CINT = matrix(c('cint1', 'cint2'), 2, 1), MANIFESTMEANS = 0)
model$pars$indvarying <- FALSE
fit <- ctFit(ctsemtestdat, model)
print(ctModelIsNonlinear(fit))

# The population matrices at the default point, which the message names
# as the T0MEANS state, with time dependent predictors at zero. The
# attribute lists the cells that depend on the point.
m <- ctSummaryMatrices(fit)
print(m$estimate$DRIFT)
print(attr(m, 'contextDependent'))

# DRIFT[1, 1] at the fixed point the processes settle at, at the mean of
# the smoothed latent states, and at a chosen state with eta1 at 8 and
# eta2 at 14.
print(ctSummaryMatrices(fit, state = 'asymptotic')$estimate$DRIFT)
print(ctSummaryMatrices(fit, state = 'mean')$estimate$DRIFT)
print(ctSummaryMatrices(fit, state = c(8, 14))$estimate$DRIFT)
}
# }
```
