# Write a model in lavaan, Mplus or formula syntax

`ctModelLavaan`, `ctModelMplus` and `ctModelFormula` create the same
model object as
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
from lavaan style syntax, Mplus style syntax, or R formulas in the style
of brms. `ctModelSyntax` writes an existing model in any of the three
styles.

## Usage

``` r
ctModelSyntax(model, style = c("lavaan", "mplus", "formula"))

ctModelLavaan(
  syntax,
  type = "ct",
  id = "id",
  time = "time",
  manifesttype = NULL,
  ...
)

ctModelMplus(
  syntax,
  type = "ct",
  id = "id",
  time = "time",
  manifesttype = NULL,
  ...
)

ctModelFormula(
  ...,
  type = "ct",
  id = "id",
  time = "time",
  manifesttype = NULL,
  family = NULL
)
```

## Arguments

- model:

  A model from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
  or one of the functions here, or a fit.

- style:

  'lavaan', 'mplus' or 'formula'.

- syntax:

  Character string containing the model.

- type:

  'ct' for a continuous time model, 'dt' for discrete time.

- id:

  Name of the subject id column, or several names for nested levels,
  innermost first, as in
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).

- time:

  Name of the time column.

- manifesttype:

  Types of the indicators that are not Gaussian, as a named vector such
  as `c(y3 = 'binary', y4 = 'ordinal')`. Types are 'continuous',
  'binary', 'ordinal', 'count' and 'censored', or the numbers 0 to 4
  that
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
  uses. An unnamed vector with one entry per indicator is also accepted.
  Ordinal indicators need `ncategories`, and censored ones `censormin`
  or `censormax`, each as a named vector.

- ...:

  For `ctModelFormula`, the formulas, followed by any named arguments
  for
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).
  For `ctModelLavaan` and `ctModelMplus`, named arguments for
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).
  Useful ones are `ncategories`, `censormin`, `censormax`,
  `TIpredNames`, `tipredDefault`, `PARS` and `driftScaled`.
  `latentNames`, `manifestNames` and `TDpredNames` set the order of the
  variables. The model matrices themselves are written by the syntax.

- family:

  For `ctModelFormula`, a brms style alternative to `manifesttype`.
  Either one family for every indicator or a named list, given as a
  name, a function or a family object. `gaussian`, `bernoulli`,
  `cumulative` and `poisson` are recognised.

## Value

`ctModelLavaan`, `ctModelMplus` and `ctModelFormula` return a model of
class `ctStanModel`, as
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
does. `ctModelSyntax` returns the code that creates the model, as a
character string that prints as code. Settings that a model holds beside
its matrices, such as `covmattransform`, are written as assignments
after the model is created. A warning says so if the code does not
create the same model.

## Details

**What is estimated by default.** Loadings, effects between processes
and predictor effects are in the model only if they are written. The
rest follows
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md):

- the first indicator of each process has a loading of 1, and other
  loadings that are written are free;

- each process has a free auto effect, free noise variance and free
  initial mean and variance, the initial mean varies between subjects,
  and the noise and initial states of different processes are
  correlated;

- each indicator has a free residual variance and a free intercept, and
  the intercept varies between subjects;

- process intercepts are zero, and residuals of different indicators are
  uncorrelated.

Fix or free any of these with the statements below.

**What the statements mean.** An effect between processes is an entry of
the drift matrix. With `type = 'ct'` it is the effect of one process on
the rate of change of another, and with `type = 'dt'` the effect on its
value at the next occasion. A time dependent predictor acts as an
impulse, shifting the process at each occasion by the effect times the
value of the predictor. A covariance that is not written, or is fixed to
zero, is zero in the fitted model. Processes and indicators need
different names.
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
lets a process take the name of its indicator, but here a statement such
as `y1 ~ 1` would then not say which of the two it is about.

**lavaan style.**


    eta1 =~ y1 + y2 + y3        # eta1 is measured by y1, y2, y3
    eta1 ~ eta1 + a*eta2 + x    # effects on eta1: auto, cross (labelled a), predictor x
    eta1 ~ 1                    # free process intercept
    eta1 ~~ 0*eta2              # no noise covariance
    y1 ~~ y2                    # residual covariance
    y1 ~ 0*1                    # intercept fixed to zero
    y4 | t1 + t2                # y4 is ordinal with three categories
    a ~ age                     # age predicts the parameter a
    a ~~ a                      # a varies between subjects

A number before `*` fixes a parameter, a name labels it, and `NA` frees
it. Two parameters with the same label are equal, as are labels joined
by `==`. A label on the first indicator of a process leaves its loading
at 1, and `NA*y1 + a*y1` frees and labels it (`y1* (a)` in Mplus style).
A variable on the right of `~` that is not a process is a time dependent
predictor. `level: 2` or `level: between` starts a block in which an
indicator's name stands for its intercept and a process's name for its
intercept, so that `y1 ~ age` and `y1 ~~ 0*y1` there are about the
intercept of `y1`. With nested ids, name the level, as in
`level: school`.

**Mplus style.**


    eta1 BY y1 y2 y3;
    eta1 ON eta1 eta2 (a);      ! label a
    eta1 ON x;
    [eta1];                     ! free process intercept
    eta1 WITH eta2@0;
    y1 WITH y2;
    [y1@0];
    [y4$1 y4$2];                ! y4 is ordinal with three categories
    s | eta2 ON eta1;           ! s varies between subjects
    %BETWEEN%
    s ON age;
    y1@0;                       ! the intercept of y1 does not vary

`@` fixes a parameter, `*` frees it, and a name in parentheses labels
the parameters before it on the same line. In a discrete time model an
effect may be written with its lag, `eta1 ON eta1&1`. A continuous time
model has no lag, and one is refused there, because the value of a
lagged effect is on another scale. If the text has sections, only
`MODEL:` is read. With nested ids use `%BETWEEN school%`.

**Formula style.**


    ctModelFormula(
      y1 + y2 + y3 ~ eta1,                 # measurement
      d(eta1) ~ eta1 + eta2 + x,           # dynamics of eta1
      d(eta2) ~ 1 + eta2 + (1 | id),       # with an intercept that varies
      drift_eta1_eta2 ~ age + (1 | id),    # a formula for one parameter
      mm_y1 ~ s(age) + (0 | id))

A formula with `d()` on the left gives the effects on a process. Its
intercept is zero unless `1` is written. Other formulas with a variable
on the left give the processes that the variable measures. A formula
with a parameter name on the left describes that parameter: `1`
estimates it, a number fixes it, predictors on the right are time
independent predictors, `s()` is a smooth effect and `(1 | id)` makes it
vary. `(0 | id)` stops it varying. Random effect terms may also be
written in the other formulas, where `(1 + eta2 | id)` lets the
intercept and the effect of `eta2` vary. Fix or label a term with
`0.5*eta2`, `"a"*eta2` or `NA*eta1`.

**Parameter names.** Every parameter has a name, whether it is free or
fixed. The name is your label if you gave one, and otherwise a default
name made from the matrix and the variables. For processes `eta1` and
`eta2`, indicators `y1` and `y2` and a time dependent predictor `x`, the
default names are:

|                   |                                     |
|-------------------|-------------------------------------|
| `drift_eta1`      | auto effect of eta1                 |
| `drift_eta1_eta2` | effect of eta2 on eta1              |
| `td_eta1_x`       | effect of x on eta1                 |
| `cint_eta1`       | intercept of eta1                   |
| `diff_eta1`       | noise sd of eta1                    |
| `diff_eta2_eta1`  | noise covariance of eta1 and eta2   |
| `T0m_eta1`        | initial mean of eta1                |
| `T0var_eta1`      | initial sd of eta1                  |
| `T0var_eta2_eta1` | initial covariance of eta1 and eta2 |
| `lambda_y2_eta1`  | loading of y2 on eta1               |
| `mm_y1`           | intercept of y1                     |
| `mvary1`          | residual sd of y1                   |
| `mvary2_y1`       | residual covariance of y1 and y2    |

`print(model)` shows the names of the free parameters.

A statement with a parameter's name on the left is about that parameter.
This is how the styles reach the initial state, and random effects and
predictors of any parameter. For the initial mean of eta1, in lavaan
style:


    T0m_eta1 ~ 0*1          # fix it at 0
    T0m_eta1 ~ start*1      # give it the label start
    start ~~ start          # it varies between subjects
    start ~~ 0*start        # it is the same for every subject
    start ~ age             # it depends on age

In Mplus style these are `[T0m_eta1@0];`, `[T0m_eta1] (start);`,
`start;`, `start@0;` and `start ON age;`. In formula style they are
`T0m_eta1 ~ 0`, `T0m_eta1 ~ "start"*1`, `start ~ (1 | id)`,
`start ~ (0 | id)` and `start ~ age`. The default name can be used in
place of the label, and the statements can come in any order. A value
given to a parameter by name is on the scale ctsem uses for that cell,
which for a variance is a standard deviation. Values given with `~~`,
`WITH` or a variance statement are variances. A name that is close to a
parameter name and matches none is taken as a misspelling and refused.

**Other ctsem features.** Wherever a label can be written, a quoted
ctsem cell or a
[`ctPar`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) call
can be written instead, such as
`ctPar("a", transform = "exp(param)", sdscale = 0.5)*eta2` in lavaan and
formula style, or `eta1 ON eta2 (a|exp(param));` in Mplus style. This
gives transforms, prior scales, stated population means and state
dependent expressions. The result is an ordinary model, so `model$pars`,
[`ctModelMatrices`](https://cdriveraus.github.io/ctsem/reference/ctModelMatrices.md),
[`ctModelHigherOrder`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md)
and the other model functions apply to it. Some things belong elsewhere.
Start values are given to
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) as
`inits`, and priors with `priors = TRUE`. A stationary initial state is
`stationary = TRUE` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md). Groups
are compared with a time independent predictor. Random effects are
correlated with each other. `poprank` and `popblocks` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) give
their covariance a factor structure, and `poprank = 0` with
`popunique = TRUE` makes them independent. Higher order lags are added
with
[`ctModelHigherOrder`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md).

## See also

[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
[`ctPar`](https://cdriveraus.github.io/ctsem/reference/ctPar.md),
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)

## Examples

``` r
# The same model four ways
m1 <- ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
  latentNames = c('eta1', 'eta2'), silent = TRUE)

m2 <- ctModelLavaan('
  eta1 =~ y1
  eta2 =~ y2
  eta1 ~ eta1 + eta2
  eta2 ~ eta1 + eta2')

m3 <- ctModelMplus('
  eta1 BY y1;
  eta2 BY y2;
  eta1 eta2 ON eta1 eta2;')

m4 <- ctModelFormula(
  y1 ~ eta1,
  y2 ~ eta2,
  d(eta1) ~ eta1 + eta2,
  d(eta2) ~ eta1 + eta2)

# A cross effect that varies between subjects and depends on age
m5 <- ctModelLavaan('
  eta1 =~ y1
  eta2 =~ y2
  eta1 ~ eta1 + a*eta2
  eta2 ~ eta2
  a ~ age
  a ~~ a')

ctModelSyntax(m5, 'mplus')
ctModelSyntax(m5, 'formula')
```
