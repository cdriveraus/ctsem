# Writing a model in lavaan, Mplus or formula syntax

``` r

library(ctsem)
```

[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
describes a model as a set of matrices. If you know lavaan, Mplus or
brms, you can write the same model in a syntax close to the one you
know:

| Function | Style | Looks like |
|:---|:---|:---|
| [`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md) | matrices | `DRIFT = matrix(c('a11', 'a12', 0, 'a22'), 2, 2)` |
| [`ctModelLavaan()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md) | lavaan | `eta1 ~ eta1 + eta2` |
| [`ctModelMplus()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md) | Mplus | `eta1 ON eta1 eta2;` |
| [`ctModelFormula()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md) | R formulas, as in brms | `d(eta1) ~ eta1 + eta2` |

All four return the same kind of model object. Fitting, summaries, plots
and the functions that extend a model work the same way whichever you
use, and
[`ctModelSyntax()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
translates a model from one style to another.

## One model, four ways

Two processes, each measured by one variable, with effects in both
directions. With
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
this is the default model for a two by two loading matrix.

``` r

m_matrix <- ctModel(LAMBDA = diag(2), silent = TRUE,
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'))

m_lavaan <- ctModelLavaan('
  eta1 =~ Y1                 # eta1 is measured by Y1
  eta2 =~ Y2
  eta1 ~ eta1 + eta2         # effects on eta1
  eta2 ~ eta1 + eta2         # effects on eta2
')

m_mplus <- ctModelMplus('
  eta1 BY Y1;
  eta2 BY Y2;
  eta1 ON eta1 eta2;
  eta2 ON eta1 eta2;
')

m_formula <- ctModelFormula(
  Y1 ~ eta1,                 # Y1 measures eta1
  Y2 ~ eta2,
  d(eta1) ~ eta1 + eta2,     # change in eta1
  d(eta2) ~ eta1 + eta2)

all.equal(m_matrix$pars, m_lavaan$pars)
#> [1] TRUE
all.equal(m_matrix$pars, m_mplus$pars)
#> [1] TRUE
all.equal(m_matrix$pars, m_formula$pars)
#> [1] TRUE
```

A regression between processes is an effect over time. In a continuous
time model (`type = 'ct'`, the default) `eta1 ~ eta2` is the effect of
the current level of `eta2` on the rate of change of `eta1`, an entry of
the drift matrix. With `type = 'dt'` it is the effect of `eta2` at one
occasion on `eta1` at the next. The effect of a process on itself is its
auto effect.

The two kinds of effect are on different scales. An autoregression of
0.8 over one unit of time corresponds to a continuous time auto effect
of log(0.8), about -0.22, and a positive auto effect in continuous time
is a process that grows. Mplus writes a discrete time effect with its
lag, `eta1 ON eta1&1`.
[`ctModelMplus()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
reads that with `type = 'dt'`, and stops if it meets a lag in a
continuous time model.

Processes and indicators need different names.
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
lets a process take the name of its indicator, but in these styles a
statement such as `Y1 ~ 1` would then not say which of the two it is
about. A process measured by `Y1` could be called `eta_Y1`, for example.

Printing the model shows what will be estimated.

``` r

m_lavaan
#> ctsem model (continuous time)
#>   2 latent (eta1, eta2), 2 manifest (Y1, Y2)
#> 
#> 13 free parameters
#> 
#> 4 individually varying: T0m_eta1, T0m_eta2, mm_Y1, mm_Y2
#>   This is a multilevel model. It estimates a 4 x 4 population covariance
#>   over those parameters -- 4 standard deviations and 6 correlations on top of the 13 above.
#>   Set indvarying = FALSE in a cell (see ?ctPar) for any that should not vary.
#> 
#> T0MEANS
#>              
#> eta1 T0m_eta1
#> eta2 T0m_eta2
#> 
#> T0VAR
#>                 eta1       eta2
#> eta1      T0var_eta1          0
#> eta2 T0var_eta2_eta1 T0var_eta2
#> 
#> DRIFT
#>                 eta1            eta2
#> eta1      drift_eta1 drift_eta1_eta2
#> eta2 drift_eta2_eta1      drift_eta2
#> 
#> DIFFUSION
#>                eta1      eta2
#> eta1      diff_eta1         0
#> eta2 diff_eta2_eta1 diff_eta2
#> 
#> CINT
#>       
#> eta1 0
#> eta2 0
#> 
#> LAMBDA
#>    eta1 eta2
#> Y1    1    0
#> Y2    0    1
#> 
#> MANIFESTMEANS
#>         
#> Y1 mm_Y1
#> Y2 mm_Y2
#> 
#> MANIFESTVAR
#>        Y1     Y2
#> Y1 mvarY1      0
#> Y2      0 mvarY2
```

## What is estimated when you say nothing

Loadings, effects between processes and effects of predictors are in the
model only if you write them. Everything else follows the defaults of
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).

| Part of the model | Default |
|:---|:---|
| Loading of the first indicator of a process | fixed to 1 |
| Other loadings that are written | free |
| Auto effect of each process | free |
| Effects between processes | zero unless written |
| Process noise (diffusion) variances and covariances | free |
| Initial means of the processes | free, and different for each subject |
| Initial variances and covariances of the processes | free |
| Intercepts of the indicators | free, and different for each subject |
| Intercepts of the processes | zero |
| Residual variances of the indicators | free |
| Residual covariances of the indicators | zero |

Some of these differ from what lavaan or brms users may expect. The
intercepts of the indicators and the initial means of the processes vary
between subjects, because people usually differ in their typical levels
and in where they start. The noise of different processes is also
correlated. Each of these can be switched off, as shown below.

## Fitting

The model is fitted with
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) like
any other. The data are in long format with one row per subject and
occasion. Here a time dependent predictor `TD1` affects the first
process. A time dependent predictor acts as an impulse. At each occasion
the process shifts by the effect times the value of the predictor, and
the dynamics carry the shift forward. This is also what `d(eta1) ~ TD1`
means in formula style. The example data set is small, so the fit uses
the default priors.

``` r

model <- ctModelLavaan('
  eta1 =~ Y1
  eta2 =~ Y2
  eta1 ~ eta1 + eta2 + TD1
  eta2 ~ eta1 + eta2
')
fit <- ctFit(ctsemtestdat, model, priors = TRUE, cores = 2)
round(summary(fit)$popmeans, 2)
#>                  mean   sd  2.5%   50% 97.5%  mode
#> T0m_eta1         4.30 1.48  1.40  4.32  7.19  4.38
#> T0m_eta2         7.44 3.41  0.61  7.49 14.01  7.65
#> drift_eta1      -0.47 0.27 -1.08 -0.42 -0.13 -0.42
#> drift_eta1_eta2 -0.10 0.14 -0.37 -0.10  0.17 -0.09
#> drift_eta2_eta1  0.39 0.24 -0.04  0.39  0.87  0.39
#> drift_eta2      -0.37 0.21 -0.91 -0.32 -0.10 -0.31
#> diff_eta1        1.98 0.20  1.62  1.96  2.42  1.98
#> diff_eta2_eta1   0.48 0.10  0.27  0.48  0.68  0.48
#> diff_eta2        2.41 0.21  2.03  2.40  2.86  2.41
#> mvarY1           0.61 0.09  0.45  0.60  0.81  0.60
#> mvarY2           0.46 0.14  0.26  0.44  0.78  0.45
#> mm_Y1            2.13 1.45 -0.73  2.12  5.01  2.05
#> mm_Y2            2.33 3.42 -4.26  2.26  9.20  2.12
#> td_eta1_TD1      1.09 0.26  0.58  1.10  1.59  1.09
```

The rows are named by parameter. `drift_eta1_eta2` is the effect of
`eta2` on `eta1`, and `td_eta1_TD1` the effect of the predictor. These
names can also be used in the syntax, as the section on parameter names
shows.

## Fixing, freeing and labelling

A number fixes a parameter, a name labels it, and two parameters with
the same label are equal.

``` r

m <- ctModelLavaan('
  eta1 =~ NA*Y1 + Y2           # NA frees the first loading
  eta2 =~ Y3 + 0.8*Y4          # a loading fixed to 0.8
  eta1 ~ a*eta1 + eta2         # auto effect labelled a
  eta2 ~ a*eta2                # the same label: equal auto effects
  eta1 ~~ 0.25*eta1            # noise variance of eta1 fixed to 0.25
  eta1 ~~ 0*eta2               # uncorrelated noise
  Y1 ~~ Y2                     # a residual covariance
  Y3 ~ 0*1                     # intercept of Y3 fixed to zero
  eta2 ~ 1                     # free intercept for the process eta2
')
```

The same in Mplus style, where `@` fixes, `*` frees and a name in
parentheses labels the parameters before it on the same line:

``` r

m_mplus <- ctModelMplus('
  eta1 BY Y1* Y2;
  eta2 BY Y3 Y4@0.8;
  eta1 ON eta1 (a)
          eta2;
  eta2 ON eta2 (a);
  eta1@0.25;
  eta1 WITH eta2@0;
  Y1 WITH Y2;
  [Y3@0];
  [eta2];
')
all.equal(m$pars, m_mplus$pars)
#> [1] TRUE
```

In formula style the modifier stands before the term, and a label is
quoted. An intercept is written as `1`, removed with `0`, or fixed with
a number. A `d()` formula has no intercept unless `1` is written.

``` r

m_formula <- ctModelFormula(
  Y1 ~ NA*eta1,
  Y2 ~ eta1,
  Y3 ~ 0 + eta2,
  Y4 ~ 0.8*eta2,
  d(eta1) ~ "a"*eta1 + eta2,
  d(eta2) ~ 1 + "a"*eta2,
  diff_eta1 ~ 0.5,             # a parameter by name: the noise sd of eta1
  diff_eta2_eta1 ~ 0,          # no noise covariance
  mvarY2_Y1 ~ 1)               # estimate the residual covariance
all.equal(m$pars, m_formula$pars)
#> [1] TRUE
```

Formula style has no `~~`, so variances and covariances are set through
their parameter names, which are described below. In lavaan and Mplus
style a value given with `~~`, `WITH` or a variance statement is a
variance or covariance, and only zero can be given for a covariance.

A covariance that is not written, or is fixed to zero, is zero in the
fitted model. With `Y1 ~~ Y2` and `Y2 ~~ Y3` in the model, `Y1` and `Y3`
have a residual covariance of zero.

As in lavaan and Mplus, a label on the first indicator of a process
names its loading and leaves it at 1. Write `NA*Y1 + a*Y1` to free and
label it, or `Y1* (a)` in Mplus style.

## Parameter names

Every parameter has a name, whether it is free or fixed. The name is
your label if you gave one, and otherwise a default name made from the
matrix and the variables. For processes `eta1` and `eta2`, indicators
`Y1` and `Y2`, and a time dependent predictor `TD1`, the default names
are:

| Parameter                               | Default name      |
|:----------------------------------------|:------------------|
| Auto effect of `eta1`                   | `drift_eta1`      |
| Effect of `eta2` on `eta1`              | `drift_eta1_eta2` |
| Effect of `TD1` on `eta1`               | `td_eta1_TD1`     |
| Intercept of `eta1`                     | `cint_eta1`       |
| Noise sd of `eta1`                      | `diff_eta1`       |
| Noise covariance of `eta1` and `eta2`   | `diff_eta2_eta1`  |
| Initial mean of `eta1`                  | `T0m_eta1`        |
| Initial sd of `eta1`                    | `T0var_eta1`      |
| Initial covariance of `eta1` and `eta2` | `T0var_eta2_eta1` |
| Loading of `Y2` on `eta1`               | `lambda_Y2_eta1`  |
| Intercept of `Y1`                       | `mm_Y1`           |
| Residual sd of `Y1`                     | `mvarY1`          |
| Residual covariance of `Y1` and `Y2`    | `mvarY2_Y1`       |

In the name of an effect, the process affected comes first. In the name
of a covariance, the variable that comes later in the model comes first.
`print(model)` shows the names of the free parameters, and
`summary(fit)` reports the estimates under the same names.

A statement with a parameter’s name on the left is about that parameter.
This is how the syntax reaches the parts of a ctsem model that lavaan,
Mplus and brms have no words for, such as the initial state of a
process. The table shows the statements for the initial mean of `eta1`.

| To | lavaan style | Mplus style |
|:---|:---|:---|
| Estimate it | `T0m_eta1 ~ 1` | `[T0m_eta1];` |
| Fix it at 0 | `T0m_eta1 ~ 0*1` | `[T0m_eta1@0];` |
| Give it the label `start` | `T0m_eta1 ~ start*1` | `[T0m_eta1] (start);` |
| Let it vary between subjects | `start ~~ start` | `start;` |
| Make it the same for every subject | `start ~~ 0*start` | `start@0;` |
| Let it depend on the covariate `TI1` | `start ~ TI1` | `start ON TI1;` |

The last three statements could also use the default name `T0m_eta1`.
Once a parameter has a label, `print(model)` and `summary(fit)` show the
label, and the statements about it can come in any order. In formula
style, one formula for the parameter can say all of this, as in
`T0m_eta1 ~ "start"*1 + TI1 + (0 | id)`.

``` r

m <- ctModelLavaan('
  eta1 =~ Y1
  T0m_eta1 ~ start*1           # label the initial mean of eta1 "start"
  start ~~ 0*start             # the same for every subject
  start ~ TI1                  # and depending on TI1
  drift_eta1 ~~ drift_eta1     # the auto effect varies between subjects
')

m_mplus <- ctModelMplus('
  eta1 BY Y1;
  [T0m_eta1] (start);
  start@0;
  start ON TI1;
  drift_eta1;
')

m_formula <- ctModelFormula(
  Y1 ~ eta1,
  d(eta1) ~ eta1,
  T0m_eta1 ~ "start"*1 + TI1 + (0 | id),
  drift_eta1 ~ (1 | id))

all.equal(m$pars, m_mplus$pars)
#> [1] TRUE
all.equal(m$pars, m_formula$pars)
#> [1] TRUE
```

A value given to a parameter by name is on the scale ctsem uses for that
cell. For `diff_eta1`, `T0var_eta1` and `mvarY1` this is a standard
deviation, while `~~`, `WITH` and an Mplus variance statement take a
variance. A covariance cell holds a covariance parameter, which
[`ctCov()`](https://cdriveraus.github.io/ctsem/reference/ctCov.md)
computes from a covariance matrix. Fixing it at 0 makes the covariance
zero.

## Differences between subjects

Any parameter can vary between subjects, and any parameter can depend on
a subject level covariate (a time independent predictor). In lavaan and
Mplus style, label the parameter and then write statements about the
label.

``` r

m <- ctModelLavaan('
  eta1 =~ Y1
  eta2 =~ Y2
  eta1 ~ eta1 + cross*eta2
  eta2 ~ eta2
  cross ~~ cross               # the cross effect varies between subjects
  cross ~ TI1                  # and depends on the covariate TI1
')
```

Mplus users will recognise a random slope, followed by a between level
model. Inside `%BETWEEN%` (or `level: 2` in lavaan style) the name of an
indicator stands for its intercept.

``` r

m_mplus <- ctModelMplus('
  eta1 BY Y1;
  eta2 BY Y2;
  cross | eta1 ON eta2;
  eta2 ON eta2;
  %BETWEEN%
  cross ON TI1;
')
all.equal(m$pars, m_mplus$pars)
#> [1] TRUE

# intercepts that depend on TI1, and no random intercept for Y2
m_between <- ctModelMplus('
  eta1 BY Y1;
  eta2 BY Y2;
  %BETWEEN%
  Y1 Y2 ON TI1;
  Y2@0;
')
```

In formula style, random effects are written as in lme4 and brms. A term
such as `(eta2 | id)` lets the effect of `eta2` vary over `id`.
Predictors of a parameter go in a formula of their own, with the
parameter on the left.

``` r

m_formula <- ctModelFormula(
  Y1 ~ eta1,
  Y2 ~ eta2,
  d(eta1) ~ eta1 + "cross"*eta2 + (eta2 | id),
  d(eta2) ~ eta2,
  cross ~ TI1)
all.equal(m$pars, m_formula$pars)
#> [1] TRUE
```

A formula for a parameter can say everything about it. `1` estimates it,
a number fixes it, `(1 | id)` makes it vary, `(0 | id)` stops it
varying, and `s()` gives a smooth effect of a covariate.

``` r

m <- ctModelFormula(
  Y1 ~ eta1,
  Y2 ~ eta2,
  d(eta1) ~ eta1 + eta2,
  d(eta2) ~ eta2,
  drift_eta1_eta2 ~ TI1 + (1 | id),   # varies, and depends on TI1
  mm_Y1 ~ s(TI2, k = 5),              # smooth effect on the intercept of Y1
  mm_Y2 ~ (0 | id))                   # no random intercept for Y2
```

All random effects in a model are correlated with each other. `poprank`
and `popblocks` in
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) give
their covariance a simpler structure.

## Indicators that are not Gaussian

The type of each indicator is an argument. Binary, ordinal, count and
censored indicators are available.

``` r

m <- ctModelLavaan('
  ability =~ item1 + item2 + rating + errors',
  manifesttype = c(item1 = 'binary', item2 = 'binary', rating = 'ordinal',
    errors = 'count'),
  ncategories = c(rating = 5))
m$manifesttype
#> [1] 1 1 2 3
```

Thresholds written in the native way also declare a binary or ordinal
indicator: `rating | t1 + t2 + t3 + t4` in lavaan style, and
`[rating$1 rating$2 rating$3 rating$4];` in Mplus style.
[`ctModelFormula()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
also takes `family`, as in
`family = list(item1 = 'bernoulli', errors = 'poisson')`. See the
vignette on measurement models for what each type assumes.

## More than one level

If subjects are nested in groups, give both id columns, innermost first,
and name the level a parameter varies at.

``` r

m <- ctModelLavaan(id = c('person', 'school'), '
  eta1 =~ Y1
  level: person
  drift_eta1 ~~ drift_eta1       # auto effect varies between persons
  level: school
  Y1 ~~ Y1                       # intercept varies between schools only
')

m_formula <- ctModelFormula(id = c('person', 'school'),
  Y1 ~ eta1 + (1 | school),
  d(eta1) ~ eta1 + (eta1 | person))
all.equal(m$pars, m_formula$pars)
#> [1] TRUE
```

In Mplus style the blocks are `%BETWEEN person%` and `%BETWEEN school%`.

## Other ctsem features

A cell of a ctsem model can also carry a transform, a prior scale for
its random effect, a stated population mean, or an expression of the
latent states. Wherever a label can be written, a
[`ctPar()`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) call
or a quoted ctsem cell can be written instead.

``` r

m <- ctModelLavaan(PARS = 'damp', '
  eta1 =~ Y1
  eta2 =~ Y2
  # a positive cross effect, with a tighter prior on how much it varies
  eta1 ~ eta1 + ctPar(cross, transform = "exp(param)", indvarying = TRUE, sdscale = 0.5)*eta2
  # an auto effect that depends on the level of the process
  eta2 ~ "-exp(damp) * (1 + eta2^2)"*eta2
')
```

In Mplus style the same goes in the parentheses that hold a label, as in
`eta1 ON eta2 (cross|exp(param)|TRUE|0.5);`. The arguments of
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
that are not matrices can be passed to all three functions, for example
`driftScaled`, `PARS` and `tipredDefault`.

The result is an ordinary model, so everything that works on a model
from
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
works here.

``` r

m <- ctModelFormula(Y1 ~ eta1, Y2 ~ eta2,
  d(eta1) ~ eta1 + eta2, d(eta2) ~ eta2)
m2 <- ctModelHigherOrder(m, indices = 1)   # second order dynamics for eta1
m$pars$indvarying[m$pars$param %in% 'drift_eta2'] <- TRUE   # edit a cell directly
```

A model saved with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
can also be opened in the ctsem GUI (the ctsemGUI package), which shows
it as a diagram and as matrices.

## Translating a model

[`ctModelSyntax()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
writes any model in any of the three styles, as code that creates it.
Use it to see how a model you already have looks in another style, or to
learn the names of its parameters.

``` r

m <- ctModelLavaan('
  eta1 =~ Y1 + Y2
  eta2 =~ Y3
  eta1 ~ eta1 + cross*eta2 + TD1
  eta2 ~ eta2
  eta1 ~~ 0*eta2
  cross ~ TI1
  cross ~~ cross
  level: 2
  Y3 ~~ 0*Y3
')
ctModelSyntax(m, 'mplus')
#> ctModelMplus(type = "ct", id = "id", time = "time", TIpredNames = c("TI1"), syntax = '
#>   eta1 BY Y1 Y2;
#>   eta2 BY Y3;
#>   eta1 ON eta1
#>           eta2 (cross)
#>           TD1;
#>   eta2 ON eta2;
#>   eta2 WITH eta1@0;
#>   cross ON TI1;
#>   cross;
#>   mm_Y3@0;
#> ')
ctModelSyntax(m, 'formula')
#> ctModelFormula(
#>   Y1 ~ eta1,
#>   Y2 ~ eta1,
#>   Y3 ~ eta2,
#>   d(eta1) ~ eta1 + "cross"*eta2 + TD1,
#>   d(eta2) ~ eta2,
#>   diff_eta2_eta1 ~ 0,
#>   cross ~ TI1 + (1 | id),
#>   mm_Y3 ~ (0 | id),
#>   type = "ct", id = "id", time = "time", TIpredNames = c("TI1"))
```

It also works on a model built with
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md):

``` r

ctModelSyntax(m_matrix, 'lavaan')
#> ctModelLavaan(type = "ct", id = "id", time = "time", syntax = '
#>   eta1 =~ Y1
#>   eta2 =~ Y2
#>   eta1 ~ eta1 + eta2
#>   eta2 ~ eta1 + eta2
#> ')
```

Settings that a model holds beside its matrices are written as
assignments after the model is created.
[`ctModelSyntax()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
compares the model its code creates with the one it was given, and warns
if they differ.

``` r

m_matrix$covmattransform <- 'cholesky'
ctModelSyntax(m_matrix, 'lavaan')
#> model <- ctModelLavaan(type = "ct", id = "id", time = "time", syntax = '
#>   eta1 =~ Y1
#>   eta2 =~ Y2
#>   eta1 ~ eta1 + eta2
#>   eta2 ~ eta1 + eta2
#> ')
#> model$covmattransform <- "cholesky"
#> model
```

## What is done elsewhere

Some things that are part of the model syntax in lavaan, Mplus or brms
are set in another place in ctsem.

| You may be looking for | In ctsem |
|:---|:---|
| Start values ([`start()`](https://rdrr.io/r/stats/start.html), `*0.5`) | `ctFit(inits = )` |
| Priors (`prior()`, `MODEL PRIORS`) | `ctFit(priors = TRUE)`, and `sdscale` in [`ctPar()`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) |
| Bounds on a parameter (`a > 0`) | a transform, `ctPar("a", transform = "exp(param)")` |
| Defined parameters (`:=`, `MODEL CONSTRAINT`) | computed from the fit, for example with [`ctSummaryMatrices()`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md) |
| Multiple groups (`group:`, `GROUPING`) | a time independent predictor of the parameters that may differ |
| A factor model for the random effects | `poprank` and `popblocks` in [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) |
| Uncorrelated random effects (`||`) | `ctFit(poprank = 0, popunique = TRUE)` |
| Lags beyond one (`y&2`), cycles, trends | [`ctModelHigherOrder()`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md), [`ctModelPeriodicComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md), [`ctModelSlowComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelSlowComponent.md) |
| Stationarity | `ctFit(stationary = TRUE)` |
| Interactions of processes (`XWITH`) | a state dependent cell, see [`?ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md) |
| Estimator and sampling settings (`ANALYSIS:`) | arguments of [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) |

If the syntax contains one of these, the functions stop and say where it
is set. An Mplus input file can be pasted whole, and only its `MODEL:`
section is read.
