# Extract samples from a ctsem fit

Returns the posterior (or approximate posterior) draws of a fit's
parameters, as arrays named for the model matrices.

## Usage

``` r
ctExtract(
  object,
  subjectMatrices = FALSE,
  cores = getOption("mc.cores", 2L),
  nsamples = "all",
  subjects = "all",
  state = NULL,
  realid = TRUE,
  ...
)
```

## Arguments

- object:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
  either backend: sampled, or optimized with its approximate posterior
  draws.

- subjectMatrices:

  Calculate subject specific system matrices? On a `backend='julia'` fit
  with `intoverpop='laplace'`, or a sampled fit, each draw's subject
  matrices combine the population draw with a draw of the subject's own
  effects from their posterior given its data. Their spread over draws
  then shows the subject's full uncertainty, as
  [`ctSubjectPars`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md)
  reports by default.

- cores:

  Only used if subjectMatrices = TRUE . For faster computation use more
  cores.

- nsamples:

  either 'all' or an integer denoting number of random samples to
  extract.

- subjects:

  either 'all', or the subjects to extract, as ids in the data (see
  `realid`).

- state:

  The latent state at which state dependent cells of the `pop_*`
  matrices are evaluated: `'T0MEANS'` (the default), `'mean'`,
  `'asymptotic'`, or a numeric state vector. Values other than the
  default need `backend='julia'` and are an error for a stan fit. No
  effect for a linear model. Only the `pop_*` arrays are affected. The
  subject matrices from `subjectMatrices=TRUE` come from the Kalman
  filter and are at each subject's last observed row.

- realid:

  If TRUE (the default), subjects are read as the ids in the data, and
  integers that are not ids there are read as ctsem's numbering of the
  subjects, 1 to N, with a message. If FALSE, they are read as that
  numbering.

- ...:

  not used.

## Value

A list of posterior sample arrays. The point at which the `pop_*`
matrices were evaluated is recorded in the `'evaluatedAt'` attribute,
without a message.
[`ctModelIsNonlinear`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md)
shows whether the point matters for the model.

Where parameters vary over subjects, the list also holds, per draw,
`popsd` (draws by parameters), `rawpopcov` and `rawpopcorr` (the
covariance and correlations on the raw scale, draws by parameters by
parameters) and `popcov` (the covariance on the natural scale). On a
`backend='julia'` fit they are named by parameter, and a model with
several grouping levels gives each level's as `rawpopcov.<level>` and so
on. There, `popcov` combines `popsd` with the raw correlations, which is
exact where a parameter's transform is affine.

## Examples

``` r
# The draws of every population matrix, an array per matrix with the draw
# first.
e <- ctExtract(ctsemtestfit)
dim(e$pop_DRIFT)

# Where the matrices were evaluated, which matters only for a model whose
# matrices depend on the state (see ?ctContextDependence).
attr(e, 'evaluatedAt')
```
