# Get Kalman filter estimates from a ctStanFit object

Get Kalman filter estimates from a ctStanFit object. `ctStanKalman` is
an alias for this function. `ctKalman`, despite the similar name, is an
alias for
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md),
which returns a long data frame.

## Usage

``` r
ctKalmanArray(
  fit,
  nsamples = NA,
  pointest = TRUE,
  collapsefunc = NA,
  cores = getOption("mc.cores", 2L),
  subjects = "all",
  timestep = "asdata",
  maxtime = "asdata",
  standardisederrors = FALSE,
  subjectpars = TRUE,
  indvarstates = FALSE,
  removeObs = F,
  realid = TRUE,
  newdata = NULL,
  ...
)
```

## Arguments

- fit:

  fit object from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), from
  any backend.

- nsamples:

  either NA (to extract all) or a positive integer from 1 to maximum
  samples in the fit.

- pointest:

  If TRUE, uses the posterior mode as the single sample.

- collapsefunc:

  function to apply over samples, such as `mean`

- cores:

  Integer number of cpu cores to use. Only needed if savescores was set
  to FALSE when fitting. `backend='stan'` only, and an error for a julia
  fit.

- subjects:

  'all' (the default), or the subjects to compute for, as ids in the
  data (see `realid`). The output numbers them in the order given,
  sorted.

- timestep:

  Either a positive numeric value, 'asdata' to use the times in the
  dataset, or 'auto' to select a timestep automatically (resulting in
  some interpolation but not excessive computation).

- maxtime:

  only relevant if timestep is not 'asdata'. Positive numeric denoting
  max time for computations.

- standardisederrors:

  If TRUE, computes standardised errors for prior, upd, smooth
  conditions.

- subjectpars:

  if TRUE, state estimates are not returned, instead, predictions of
  each subjects parameters are returned, for parameters that had random
  effects specified.

- indvarstates:

  if TRUE, do not remove indvarying states from output

- removeObs:

  Logical or integer. If TRUE, observations (but not covariates) are set
  to NA, so only expectations based on parameters and covariates are
  returned. If a positive integer N, every N observations are retained
  while others are set NA for computing model expectations – useful for
  observing prediction performance forward further in time than one
  observation.

- realid:

  If TRUE (the default), subjects are read as the ids in the data, and
  integers that are not ids there are read as ctsem's numbering of the
  subjects, 1 to N, with a message. If FALSE, they are read as that
  numbering.

- newdata:

  Optional long format data, in the same layout as the data the model
  was fitted to, to compute for instead of the fit's own data. It can
  hold new subjects, the fit's subjects with other covariate values, or
  a record cut short. The fit's estimates (or, with `pointest=FALSE`,
  its draws) are used without refitting. See
  [`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md).

- ...:

  additional arguments to collpsefunc.

## Value

list containing Kalman filter elements, each element in array of
iterations, data row, variables. llrow is the log likelihood for each
row of data.

For a binary or ordinal manifest variable the returned manifest
quantities are the \*linear predictor\* \\\nu = \Lambda \eta + \tau\\,
on the logit scale, not a probability and not a category. Turn one into
a probability with `plogis(nu)` for binary, or into cumulative category
probabilities with `plogis(tau_k - nu)` for ordinal, taking the
thresholds from `summary(fit)$popmeans`. The manifest covariance for
such a variable is likewise the covariance of the linear predictor.

## Examples

``` r
k=ctKalmanArray(ctstantestfit,subjectpars=TRUE,collapsefunc=mean)
```
