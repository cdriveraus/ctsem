# ctGenerate

This function generates data according to the specified ctsem model
object.

## Usage

``` r
ctGenerate(
  ctmodelobj,
  n = 100,
  burnin = 0,
  dtmean = 1,
  logdtsd = 0,
  dtmat = NA,
  Tpoints = NULL,
  wide = FALSE,
  backend = c("auto", "r", "julia", "stan"),
  intoverstates = "auto",
  fromPriors = FALSE,
  nsamples = 200,
  datastruct = NA,
  cores = getOption("mc.cores", 2L),
  values = NULL,
  n.subjects
)
```

## Arguments

- ctmodelobj:

  ctsem model object from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).

- n:

  Number of subjects to output. For a model with grouping levels above
  the subject (`ctModel(id = c('subject', 'study'))`), one count per
  level in the same order, innermost first. For example `c(100, 10)`
  gives 100 subjects split evenly over 10 studies. A level not given has
  one group.

- burnin:

  Number of initial time points to discard (to simulate stationary data)

- dtmean:

  Positive numeric. Median time interval (delta T) to use. Intervals are
  drawn as `exp(rnorm(n, log(dtmean), logdtsd))`, so the mean interval
  is `dtmean * exp(logdtsd^2/2)`, which is longer than `dtmean` when
  `logdtsd` is not zero.

- logdtsd:

  Numeric. Standard deviation of the log time interval. Zero gives an
  equal-interval design.

- dtmat:

  Either NA, or numeric matrix of n rows and Tpoints-1 columns,
  containing positive numeric values for all time intervals between
  measurements. If not NA, dtmean and logdtsd are ignored.

- Tpoints:

  Optional number of time points to generate. If supplied, this
  overrides any `Tpoints` stored in `ctmodelobj`. If not supplied,
  `ctGenerate` uses `ctmodelobj$Tpoints` when available.

- wide:

  Logical. Output in wide format?

- backend:

  Which generator to use: `'r'`, `'julia'`, or `'auto'` (the default).
  `'julia'` generates from the model itself. Each subject's parameters
  are drawn from the population at every level, the latent trajectory is
  simulated from the process, and each observation is drawn from its
  measurement model. Nonlinear dynamics, random effects on any
  parameter, and binary, ordinal and count indicators are all supported.
  It needs a `ctModel(type='ct')` or `'dt'` model. `'r'` simulates a
  linear Gaussian system at the values the model states, with individual
  differences only through `TRAITVAR` and `MANIFESTTRAITVAR`, and fills
  free parameters with fixed defaults. It does not support nonlinear
  models or non-Gaussian indicators, and reports which random effects or
  TI predictor effects of the model it leaves out. `'auto'` uses `'r'`
  unless the model or the call needs `'julia'`: a nonlinear model, a
  non-Gaussian indicator, initial conditions set to `'stationary'`,
  stated population means (`ctPar(mean=)`), more than one count in `n`,
  `values`, or `datastruct`.

- intoverstates:

  For `backend='julia'`: `'auto'` (the default, meaning `FALSE`), `TRUE`
  or `FALSE`. Mainly for checking the Kalman filter's approximation.
  `FALSE` simulates the latent trajectory and then each observation
  given the state at its row. `TRUE` draws each row from the filter's
  one-step-ahead prediction for the same subjects, the density that a
  fit integrating over the states maximises. The two agree in
  distribution for linear dynamics with Gaussian indicators.
  `backend='r'` always simulates the states and ignores the argument.

- fromPriors:

  Draw the parameters from the model's priors, and return `nsamples`
  datasets instead of one. This gives the prior predictive distribution
  and needs no fit or data. Without a `datastruct`, a balanced design of
  `n` subjects at `Tpoints` occasions `dtmean` apart is used.

- nsamples:

  With `fromPriors=TRUE`, how many datasets to draw.

- datastruct:

  A long format data frame giving the design to generate for, in place
  of the balanced design from `n`, `Tpoints` and `dtmean`. It gives the
  subject and any group ids, the times, TD and TI predictor values, and
  which observations exist. A manifest value that is NA stays missing,
  and any other value, or a manifest column left out, is generated. A TD
  predictor the model names but `datastruct` lacks is set to zero, and a
  TI predictor it lacks is drawn standard normal per subject. The same
  data frame is returned with the generated values filled in. Needs
  `backend='julia'`, or `fromPriors=TRUE`, where a fit passed as
  `ctmodelobj` supplies its own design.

- cores:

  With `fromPriors=TRUE`, cpu cores to use.

- values:

  For `backend='julia'`, the population values to generate from, named
  and scaled as [`summary()`](https://rdrr.io/r/base/summary.html)
  reports them. A named numeric vector gives free parameters' values on
  their natural scale, e.g. `c(mm = 2, drift = -0.3)`. A list gives
  [`summary()`](https://rdrr.io/r/base/summary.html)'s sections:
  `popmeans` (as above), `popsd` (random effect sds on the natural
  scale, named by parameter), `rawpopcorr` (random effect correlations,
  named `'a__b'`) and `tipreds` (linearised TI predictor effects, named
  `'tip_<predictor>_<parameter>'`), as in
  `list(popmeans = c(mm = 2), popsd = c(mm = 0.5), tipreds = c(tip_TI1_mm = 1))`.
  With grouping levels, `popsd.<level>` and `rawpopcorr.<level>` give a
  level's values. Values not given are set to their prior's centre: raw
  zero for a mean, `sdscale` for a raw sd, and zero for a correlation or
  a free effect. A value fixed in the model cannot also be given here.

- n.subjects:

  Deprecated; use `n`.

## Value

From the R generator, a long format matrix with one row per subject and
time point, or wide format with `wide = TRUE`. From the julia generator,
the same long format, carrying the attributes described under Generating
values. With `fromPriors = TRUE`, a list: `Y`, an `nsamples` by data
rows by manifest variables array of generated data, and `llrow`, the
matching array of row log likelihoods.

## Details

Covariance related matrices are treated as Cholesky factors.
TRAITTDPREDCOV and TIPREDCOV matrices are not accounted for, at present.
The first 1:n.TDpred rows and columns of TDPREDVAR are used for
generating tdpreds at each time point.

## Generating values

With `backend='julia'`, each value the model leaves free is set as
follows, and all of them are listed in a message:

- a free parameter is at its prior's centre, raw zero, unless `values`
  names it;

- the population sd of an individually varying parameter is its
  `sdscale` on the raw scale, at every level it varies at, unless
  `RAWPOPVAR` (`RAWPOPVAR_<level>`) gives a number. Random effects are
  uncorrelated unless `RAWPOPVAR` gives its off-diagonal values, which
  [`ctCov`](https://cdriveraus.github.io/ctsem/reference/ctCov.md)
  writes from a covariance matrix;

- a time independent predictor is drawn standard normal per subject. An
  effect fixed in the model (`'mm||||TI1=0.5'`) shifts the raw parameter
  by that much per unit of the predictor, and an effect left free is
  zero.

`fromPriors=TRUE` draws all of these from the priors.

The julia generator also returns, as attributes, what it generated from:
`attr(, 'subjectpars')`, each subject's parameters on their natural
scale with its id and groups; `attr(, 'grouppars')`, for a model with
grouping levels, each group's parameters at each level, as for a subject
at the group's average; and `attr(, 'states')`, the latent states at
every row by id and time (not with `intoverstates=TRUE`, which draws no
states).

## Examples

``` r
#a 2 process model, each process measured by a noisy indicator, with stable
#individual differences in process levels
generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
 MANIFESTVAR=diag(.1,2),
 LAMBDA=diag(1,2),
 DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
 DIFFUSION=matrix(c(1,.2,0,4),2),
 CINT=matrix(c('cint1, indvarying=TRUE, sdscale=0.05',
   'cint2, indvarying=TRUE, sdscale=0.08'),nrow=2),
 MANIFESTMEANS=matrix(0,ncol=1,nrow=2),
 T0MEANS=matrix(0,ncol=1,nrow=2),
 T0VAR=diag(1,2))

#the R generator: every subject at the values the model states, with free
#cells filled by fixed defaults and no individual differences
data <- ctGenerate(generatingModel,n=15,burnin=10)

# \donttest{
#the julia engine: each subject's CINT drawn from the population, with the
#means stated on the natural scale and the sds given by sdscale on the raw
#scale (CINT is 10 * raw, so 0.5 and 0.8)
if (isTRUE(ctJuliaStatus()$available)) {
  data <- ctGenerate(generatingModel,n=15,burnin=10,
    backend='julia',values=c(cint1=1,cint2=0))
}
# }

#the prior predictive: what the model says data could look like before it
#has seen any. No fit, and no data.
# \donttest{
priorpred <- ctGenerate(generatingModel, fromPriors = TRUE, Tpoints = 6,
  n = 10, nsamples = 20, cores = 2)
str(priorpred$Y)
# backend = 'auto' generates with julia where it is available. Julia keeps
# running between calls, holding what it has compiled; this ends it.
if (isTRUE(ctJuliaStatus()$available)) ctJuliaKill("session")
# }
```
