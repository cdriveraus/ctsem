# Define a ctsem model

This function is used to specify a continuous time structural equation
model, which can then be fit to data with function
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

## Usage

``` r
ctModel(
  LAMBDA,
  type = "ct",
  n.manifest = "auto",
  n.latent = "auto",
  Tpoints = NULL,
  manifestNames = "auto",
  manifesttype = rep(0, nrow(LAMBDA)),
  ncategories = NULL,
  censormin = NULL,
  censormax = NULL,
  latentNames = "auto",
  id = "id",
  time = "time",
  silent = FALSE,
  T0VAR = "auto",
  T0MEANS = "auto",
  MANIFESTMEANS = "auto",
  MANIFESTVAR = "diag",
  ASYMPTOTES = NA,
  RESPONSE = NA,
  DRIFT = "auto",
  CINT = 0,
  DIFFUSION = "auto",
  n.TDpred = "auto",
  TDpredNames = "auto",
  TDPREDEFFECT = "auto",
  TDPREDMEANS = "auto",
  TDPREDVAR = "auto",
  n.TIpred = "auto",
  TIpredNames = "auto",
  tipredDefault = TRUE,
  PARS = NULL,
  driftScaled = character()
)
```

## Arguments

- LAMBDA:

  n.manifest\*n.latent loading matrix relating latent to manifest
  variables, with latent processes 1:n.latent along the columns, and
  manifest variables 1:n.manifest in the rows.

- type:

  character string. Use 'ct' (continuous time) or 'dt' (discrete time)
  for fitting with
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).
  These return a modern ctsem model object whose `pars` data frame can
  also be edited through the pars-backed `matrices` view. For legacy
  workflows, 'omx' returns a matrix-list model object that can be edited
  directly and then converted via
  [`ctModelConvertOMX`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md).

- n.manifest:

  Number of manifest indicators per individual at each measurement
  occasion / time point. Manifest variables are included as the first
  element of the wide data matrix, with all the 1:n.manifest manifest
  variables at time 1 followed by those of time 2, and so on.

- n.latent:

  Number of latent processes.

- Tpoints:

  Number of time points, or measurement occasions, in the data. This
  will generally be the maximum number of time points for a single
  individual, but may be one extra if sample relative time intervals are
  used, see
  [`ctIntervalise`](https://cdriveraus.github.io/ctsem/reference/ctIntervalise.md).

- manifestNames:

  n.manifest length vector of manifest variable names as they appear in
  the data structure, without any \_Tx time point suffix that may be
  present in wide data. Defaults to Y1, Y2, etc.

- manifesttype:

  n.manifest length vector of manifest variable types: 0 (the default)
  for continuous Gaussian variables, 1 for binary (coded 0/1), 2 for
  ordinal (consecutive integers from 1), 3 for counts (non-negative
  integers, Poisson with a log link, so the latent process gives the log
  rate) and 4 for censored (Gaussian, observed between `censormin` and
  `censormax`). For a count the `MANIFESTVAR` diagonal stays free and is
  a log-scale dispersion, giving a Poisson-lognormal observation with
  variance \\\mu + \mu^2(e^{\sigma^2}-1)\\. Fix it to zero for a plain
  Poisson. A censored value at a limit contributes the probability of
  being at or beyond that limit. Types other than 0 need the julia
  backend, except binary, which stan also fits.

- ncategories:

  n.manifest length integer vector giving the number of categories of
  each ordinal manifest variable, ignored for other types. Required when
  any `manifesttype` is 2, by `ctGenerate` as well as `ctFit`.

- censormin, censormax:

  n.manifest length numeric vectors giving the lower and upper censoring
  limits of each `manifesttype = 4` variable, ignored for other types.
  The defaults are `-Inf` and `Inf`, so a variable censored on one side
  keeps the other limit infinite. The limits are fixed, known values.

- latentNames:

  n.latent length vector of latent variable names (used for naming
  parameters, defaults to eta1, eta2, etc).

- id:

  character string denoting column name containing subject
  identification variables. id data may be of any form, though will be
  coerced internally to an integer sequence rising from 1. Or a
  character vector of column names for nested grouping levels, innermost
  first: `id = c('subject','study')` places each subject in one study.
  The `indvarying` field of a cell (see
  [`ctPar`](https://cdriveraus.github.io/ctsem/reference/ctPar.md)) sets
  the levels its parameter varies at: `TRUE` for every level, `FALSE`
  for none, or the names of levels. For example `'cint||study'` gives an
  intercept that varies between studies, and `'cint||subject,study'` (or
  `TRUE`) one that varies at both levels. `model$pars$indvarying` holds
  these values, one per cell, and can be edited directly. A cell's
  `sdscale` applies at every level it varies at. Intercepts (`T0MEANS`,
  `MANIFESTMEANS`, `CINT`) that are left unspecified vary at every
  level. Levels above the subject require `backend='julia'` with
  `intoverpop='laplace'`, and the levels must be strictly nested.

- time:

  character string denoting column name containing timing data. Timing
  data must be numeric.

- silent:

  Suppress all output to console.

- T0VAR:

  lower triangular n.latent\*n.latent cholesky matrix of latent process
  initial variance / covariance. "auto" freely estimates all parameters.

- T0MEANS:

  n.latent\*1 matrix of latent process means at first time point, T0.
  "auto" freely estimates all parameters.

- MANIFESTMEANS:

  n.manifest\*1 matrix of manifest intercept parameters. "auto" frees
  all parameters.

- MANIFESTVAR:

  lower triangular n.manifest\*n.manifest cholesky matrix of variance /
  covariance between manifests at each measurement occasion (i.e.
  measurement error / residual). "auto" freely estimates variance
  parameters, and fixes covariances between manifests to 0. "free" frees
  all values, including covariances.

- ASYMPTOTES:

  Matrix with one row per binary indicator, in manifest order, and two
  columns giving the lower and upper asymptotes of its response
  probability. The probability is then lower + (upper - lower) times the
  logistic of the linear predictor. `NA` means no asymptote on that
  side, and the default of all `NA` gives the two parameter logistic
  model. A parameter name in the lower column estimates a guessing
  probability (three parameter model), in the upper column a ceiling
  below one, and in both the four parameter model. A number between 0
  and 1 fixes the asymptote, and a name repeated across rows estimates
  one asymptote for those items. A length 2 vector applies to every
  binary indicator, so `c('guess', NA)` estimates one guessing parameter
  for the whole test. Julia backend only.

- RESPONSE:

  Matrix with one row per manifest variable and two columns, setting the
  response scale of censored indicators (manifesttype 4). Rows for other
  indicators keep their defaults. Column 1, `halfwidth`, is half the
  width of the interval a recorded value stands for. A value `y` means
  the latent response lies between the transformed `y - halfwidth` and
  `y + halfwidth`, a value at `censormin` that it lies below
  `censormin + halfwidth`, and a value at `censormax` that it lies above
  `censormax - halfwidth`. The default, 0, treats values as exact. For a
  rating scale use half the step between categories. A step that differs
  between studies can be written as a time independent predictor in the
  cell. Column 2, `power`, sets the transform from the recorded scale to
  the latent response, \\G(r) = L + ((1 + r - L)^p - 1)/p\\, with \\L\\
  the censormin (0 if there is none). The default, 1, is linear. Values
  below 1 make the latent response change less per recorded unit towards
  the top of the scale, and 0 gives the log transform. A parameter name
  estimates the power, and a name shared across items or studies
  estimates one power from all of them. `NA` keeps the default, and a
  length 2 vector fills every censored row, so `c(0.5, 'power')` suits
  integer data with one estimated transform. Julia backend only.

- DRIFT:

  n.latent\*n.latent DRIFT matrix of continuous auto and cross effects,
  relating the processes over time. "auto" freely estimates all
  parameters.

- CINT:

  n.latent \* 1 matrix of latent process intercepts, allowing for non 0
  asymptotic levels of the latent processes. Generally only necessary
  for additional trends and more complex dynamics.

- DIFFUSION:

  lower triangular n.latent\*n.latent cholesky matrix of diffusion
  process variance and covariance (latent error / dynamic innovation).
  "auto" freely estimates all parameters.

- n.TDpred:

  Number of time dependent predictor variables in the dataset.

- TDpredNames:

  n.TDpred length vector of time dependent predictor variable names, as
  they appear in the data structure, without any \_Tx time point suffix
  that may appear in wide data. Default names are TD1, TD2, etc.

- TDPREDEFFECT:

  n.latent\*n.TDpred matrix of effects from time dependent predictors to
  latent processes. Effects from 1:n.TDpred columns TDpredictors go to
  1:n.latent rows of latent processes. "auto" freely estimates all
  parameters.

- TDPREDMEANS:

  Legacy argument for time-dependent predictor means, mainly to support
  data generation workflows. Expected dimensions are (n.TDpred \*
  Tpoints) rows and 1 column. "auto" creates a free-label matrix.

- TDPREDVAR:

  Legacy argument for time-dependent predictor covariance structure,
  mainly to support data generation workflows. Expected dimensions are
  (n.TDpred \* Tpoints) square. "auto" creates a free-label matrix.

- n.TIpred:

  Number of time independent predictors. Each TIpredictor is inserted at
  the right of the data matrix, after the time intervals.

- TIpredNames:

  n.TIpred length vector of time independent predictor variable names,
  as they appear in the data structure. Default names are TI1, TI2, etc.
  Or a formula of such names and smooth terms, such as
  `~ s(age, k = 8) + female`. A smooth term is a penalised spline of a
  subject-level covariate, with its smoothness estimated from the data
  (`backend = 'julia'`). It enters the model as `k - 1` predictor
  columns, `age_lin` (the linear trend) and `age_s1`, ..., built from
  the data when the model is fitted. To give a parameter the smooth
  effect, name the covariate in its cell's `tipreds` field
  (`'cint||||age'`). Under `tipredDefault = TRUE`, cells without a
  `tipreds` field get the smooth effect if their parameter is free and
  varies between subjects.

- tipredDefault:

  Logical. TRUE sets any parameters with unspecified time independent
  predictor effects to have effects estimated, FALSE fixes the effect to
  zero unless individually specified. A cell's effects are unspecified
  only when its `tipreds` field is left out. An empty field
  (`'mypar||||'`) means no predictor effects, which excludes a parameter
  under `tipredDefault = TRUE`. The setting is stored on the model as
  `$tipredDefault`.

- PARS:

  for types 'ct' and 'dt' only. May be of any structure, only needed to
  contain extra parameters for certain non-linear models.

- driftScaled:

  for type 'ct' only: which of "DIFFUSION" and "CINT" to parameterise by
  what they add over one time unit, given the process's own drift
  diagonal `a`. A diffusion sd `s` becomes `s * sqrt(2 a / expm1(2 a))`,
  so that `s^2` is the variance the process gains from its own noise
  over one time unit. An intercept `u` becomes `u * a / expm1(a)`. The
  model and its likelihood are unchanged. The rescaled parameters depend
  less on the drift diagonal, which helps uncertainty estimation and
  sampling when fast mean reversion and larger noise are hard to tell
  apart. Only the diagonal of DIFFUSION is rescaled, so its correlations
  keep their usual meaning. Each rescaled parameter, with its name,
  transform, random effects and predictor effects, moves to a new cell
  of `PARS`, and its original cell becomes an expression of that cell
  and the drift diagonal. Not available with
  `covmattransform = 'cholesky'`.

## Value

For `type = 'ct'` or `'dt'`, a model of class `ctStanModel`, with one
row per matrix cell in `$pars` and the same cells as matrices in
`$matrices`. Either can be edited. For `type = 'omx'`, a list of
matrices of class `ctsemInit`, for
[`ctGenerate`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)
or
[`ctModelConvertOMX`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md).

## Details

A matrix cell can be one of four things:

- a **number**, which fixes the cell;

- a **parameter name**, which frees it. The same name in two cells is
  one parameter, which is how equality constraints are written;

- a name **with fields**, giving the transform, the levels of `id` at
  which the parameter varies, the prior scale of that variation, and the
  time independent predictors that affect it;

- an **expression**, computed from the latent states, the time dependent
  predictors and the `PARS` cells.

Fields can be given in a fixed order separated by `|`
(`'mm||TRUE|0.5'`), or by name (`'mm, indvarying=TRUE, sdscale=0.5'`).
With nested ids, `indvarying` may name levels, as in `'mm||study'` or
`'mm, indvarying=c(subject, study)'`. A cell containing `|` is read in
the ordered form. Otherwise a cell with a field of the form `name=value`
is read in the named form, and an unknown field name is an error.
[`ctPar`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) writes
these strings from named arguments.

Every name in an expression must be a latent state, a time dependent
predictor, or a parameter declared in `PARS`. For example
`DRIFT = matrix('-exp(dr)*eta1')` needs `PARS = c('dr')`. To transform a
single parameter, the transform field is shorter:
`DRIFT = matrix('dr|-exp(param)')`. With the julia backend a cell may
also use a time independent predictor by name, which gives the subject's
value of it, as in `MANIFESTVAR = 'rsd + mvar'` with
`PARS = c('mvar|1e-10 + 5 * log1p_exp(2 * param)')`. Missing values of
such a predictor are imputed or sampled like other missing predictor
values (`optimcontrol$tipredMissing`). Population summaries evaluate the
cell at a predictor value of zero. A cell cannot refer to a manifest
variable or to the time interval. Predictor *effects* on a parameter are
set with the `tipreds` field.

A latent state or time dependent predictor name in a cell refers to that
state or predictor, and any number of matrices may refer to the same
one. A name declared in `PARS` likewise refers to the `PARS` cell
wherever else it appears, and that cell holds the parameter and its
transform.

Any other name repeated across cells is one parameter, with one value
and one random effect, so it needs one transform. Cells that would
transform it differently are refused, including cells whose *default*
transforms differ. The `DRIFT` diagonal is bounded negative and the
off-diagonals are not, so `DRIFT = matrix('a', 2, 2)` is refused. To use
one parameter in cells that need different transforms, declare it in
`PARS` (`PARS = c('a')`) and write the transform in each cell, as in
`DRIFT[1,1] = '-exp(a)'` and `DIFFUSION[1,1] = 'exp(a)'`. If one
transform suits every cell, give it in `PARS`
(`PARS = c('a|-exp(param)')`) and write only the name in each cell. Two
`PARS` cells with the same name are one parameter and must have the same
transform.

A predictor effect belongs to the parameter, so cells holding one
parameter must have the same `tipreds`. Naming an effect
(`'mypar||||age=shared'`) gives every parameter with that effect name,
for that predictor, one shared coefficient. Because the coefficient acts
on the raw scale, parameters with different transforms cannot share it.
The same effect name under two predictors gives two coefficients. A
number instead of a name (`'mypar||||age=0.5'`) fixes the effect,
shifting each subject's raw parameter value by 0.5 times their `age`
before the transform. Fixed effects are fitted by `backend = 'julia'`
and used by `ctGenerate`. The stan backend refuses them.

A cell that is exactly a latent state or time dependent predictor name
is read as a reference to it, with a warning, since a parameter of that
name would give a different model. Write `'state[1]'` or
`'tdpreds[rowi, 1]'` to avoid the warning. In `PARS` such a cell is an
error, because a bare name there declares a parameter. For a `PARS` cell
that depends on a state, write the state in an expression, such as
`'log1p(exp(eta1))'`.

## See also

[`ctPar`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) for the
fields a matrix cell can carry,
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) for
fitting the result.

## Examples

``` r
 ### Frequentist example:
 ### impulse and level change time dependent predictor 
 ### example from Driver, Oud, Voelkle (2015)
 data('ctExample2')
 tdpredmodel <- ctModel(type='omx', n.manifest = 2, n.latent = 3, n.TDpred = 1, 
 Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
 TDpredNames = 'MoneyInt', 
 latentNames = c('LeisureTime', 'Happiness', 'MoneyIntLatent'),
 LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3))

 tdpredmodel$DIFFUSION[, 3] <- 0
 tdpredmodel$DIFFUSION[3, ] <- 0
 tdpredmodel$T0VAR[3, ] <- 0
 tdpredmodel$T0VAR[, 3] <- 0
 tdpredmodel$CINT[3] <- 0
 tdpredmodel$T0MEANS[3] <- 0
 tdpredmodel$DRIFT[3, ] <- 0
 
 
###Bayesian example:
model<-ctModel(type='ct',
n.latent=2, latentNames=c('eta1','eta2'),
n.manifest=2, manifestNames=c('Y1','Y2'),
n.TDpred=1, TDpredNames='TD1', 
n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
LAMBDA=diag(2))
```
