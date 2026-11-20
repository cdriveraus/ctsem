# Kalman filter and smoother estimates from a julia backend fit

Prior, filtered and smoothed estimates of the latent states and the
observations, for every row of data, from the same forward pass the
engine uses for the likelihood.
[`ctKalmanArray`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md)
(alias `ctStanKalman`) dispatches here for a `ctJuliaFit`, so call that
unless you want the julia route specifically.
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
(alias `ctKalman`) returns the same information as a long data frame.

## Usage

``` r
ctBackendKalman(
  fit,
  subjects = "all",
  timestep = "asdata",
  maxtime = "asdata",
  removeObs = FALSE,
  pointest = TRUE,
  nsamples = NA,
  collapsefunc = NA,
  standardisederrors = FALSE,
  subjectpars = FALSE,
  indvarstates = FALSE,
  randomEffects = NULL,
  ...
)
```

## Arguments

- fit:

  A `ctJuliaFit`.

- subjects:

  `'all'`, a vector of subject ids, or integer positions into the fitted
  subjects.

- timestep:

  `'asdata'` to use the observed times, `'auto'` to choose an
  interpolation step, or a positive number.

- maxtime:

  Only used when interpolating: the largest time to compute for.

- removeObs:

  `TRUE` withholds every observation so that only expectations given
  parameters and covariates are returned; a positive integer N keeps
  every Nth.

- pointest:

  Use the point estimate as the single sample.

- nsamples:

  Number of posterior draws to use when `pointest=FALSE`.

- collapsefunc:

  Function applied over draws, e.g. `mean`.

- standardisederrors:

  Also return residuals standardised by the prior observation
  covariance.

- subjectpars:

  Also return each subject's own model matrices.

- indvarstates:

  Keep the augmented individual-difference states in the latent output
  rather than trimming to the real processes. On a single level of
  Laplace random effects, where the effects are not states, they are
  appended as if they were: predicted and updated from the subject's
  rows before and up to each row, smoothed from its whole record, on the
  raw scale as deviations from the population value. Their covariances
  with the latent states are not computed and are `NA`.

- randomEffects:

  For an \`intoverpop='laplace'\` fit, which levels of random effect the
  trajectories are built from. Naming one of the model's id columns
  includes that level and every level outside it, so the innermost id
  gives each subject its own trajectory and an outer id gives the group
  mean, shared by every subject in the group. \`'population'\` includes
  none, giving the trajectory implied by the population parameters
  alone. Defaults to the innermost id. Ignored for other fits.

  With one level of random effects, left in, the predicted and updated
  passes take each subject's effects at their posterior given its rows
  up to each row (\`ctsem_laplace_sequential\`), and the smoothed pass
  at their mode given its whole record. Otherwise every pass takes the
  whole-record modes. A message says which once per call.

  The whole-record modes are the ones the fit arrived at, and they stay
  fixed however this call changes the rows being filtered. Selecting
  subjects, interpolating a time grid or withholding observations with
  `removeObs` therefore leaves them alone: with every observation
  withheld the smoothed pass still uses that subject's own random
  effects.

- ...:

  Passed to `collapsefunc`.

## Value

A list of arrays, the same shape
[`ctKalmanArray`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md)
returns for Stan fits, with `subjectMatrices` attached when
`subjectpars=TRUE`.

## Details

Internal. `ctBackend*` means internal throughout the package, so this is
not exported; users reach it through
[`ctKalmanArray`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md)
and
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md),
which dispatch here for a `ctJuliaFit`. The `randomEffects` notes below
describe behaviour those two inherit.

Individual-difference parameters come out of the smoothed initial state:
ctsem carries an individually varying parameter as an augmented latent
state with no drift and no diffusion, so its smoothed t0 estimate is the
subject's value for it.

## Examples

``` r
# \donttest{
# ctBackendKalman(fit, timestep = .1)
# }
```
