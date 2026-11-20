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
`attr(x, 'contextDependent')`, and give a message once.
[`ctModelIsNonlinear`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md)
reports whether a model has such cells. With `backend='julia'` the
evaluation point can be chosen with the `state` argument of
[`ctBackendParMatrices`](https://cdriveraus.github.io/ctsem/reference/ctBackendParMatrices.md).

## See also

[`ctModelIsNonlinear`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md),
[`ctBackendParMatrices`](https://cdriveraus.github.io/ctsem/reference/ctBackendParMatrices.md)
