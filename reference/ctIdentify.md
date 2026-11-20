# Check which parameters a dataset can inform, before fitting

Reports the directions of the parameter space the data carry no
information about, and the parameters involved in them, without fitting
the model.

## Usage

``` r
ctIdentify(
  datalong,
  model,
  inits = NULL,
  nstart = 3L,
  spread = 0.5,
  priors = FALSE,
  cores = getOption("mc.cores", 2L),
  verbose = 0L,
  rtol = 1e-13,
  seed = 1L,
  ...
)
```

## Arguments

- datalong:

  Long format data, as for
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- model:

  Model from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).

- inits:

  Optional evaluation point. Supplying one uses that single point,
  without the comparison across points described in Details.

- nstart:

  Number of points to evaluate. The first is the origin of the
  unconstrained space, and the rest are dispersed draws around it.

- spread:

  Standard deviation of the additional evaluation points on the raw
  scale.

- priors:

  Whether to include the prior. `FALSE` by default, because a prior
  informs every direction, so a model the data do not identify looks
  fine with priors on.

- cores:

  As for
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- verbose:

  As for
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- rtol:

  A direction counts as uninformed when the information along it is
  below this fraction of the best informed direction's. The default
  separates zero from small, so that only structurally unidentified
  directions are reported. Raise it to report weakly identified
  directions too, including some a fit would settle.

- seed:

  Seed for the dispersed evaluation points, so that identical calls give
  identical output. The session's random stream is restored afterwards.
  `NULL` draws from the session's stream instead.

- ...:

  Further arguments to
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) that
  shape the model, such as `intoverpop` or `poprank`, with ctFit's
  defaults.

## Value

An object of class `ctIdentify`, with `$parameters` (names in the flat
directions at every evaluation point), `$partial` and `$structural`
(those split into the partially and completely unidentified),
`$sometimes` (in a flat direction at some points only), `$rotating`
(whether the flat directions change between points), `$nweak`,
`$condition`, `$starts` (the detail for each point), and `$effects`, the
information about the random effects at a supplied `inits` (`NULL`
otherwise). Printing it summarises what was found.

## Details

This finds structural non-identification, such as parameters that only
appear together, a latent process nothing measures, or a scale that is
never fixed. These are properties of the model that hold wherever it is
evaluated. A direction can be flat at one point and well informed at the
estimate, so several points are used, and only parameters in the flat
directions at every point are reported as uninformed.

The model is prepared by
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), with
any of its arguments given in `...`, so what is checked is the model a
fit with the same arguments would estimate.

A random effect on a variance parameter (DIFFUSION, MANIFESTVAR) under
`intoverpop = 'augmented'` with a full rank population covariance
(`poprank = NA`) is partially identified. The population covariances it
implies with the other varying parameters are determined, but their
split into a scale and correlations is not. The printed summary
separates partial from complete non-identification, because the advice
differs: fixing a value discards information in the first case, and is
the fix in the second.

Requires the julia backend, which supplies the per-subject scores. The
information has rank at most the number of subjects, so a model with
more free parameters than subjects shows flat directions for that reason
alone, and the printed summary says when that applies.

How well the data determine each random effect's population standard
deviation, and how much of each subject's effect its own data determine,
is checked only at a supplied `inits`, such as a fit's `estimate$raw`.
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) makes
the same check at its estimate, stored as `fit$identifiability$effects`.
It depends on the population standard deviations, so at the default
evaluation points it would describe the starting values instead of the
data.

## See also

[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) runs
the equivalent check on a finished fit, evaluated at the estimate and
stored as `fit$identifiability`.

## Examples

``` r
# \donttest{
# ctIdentify() always uses the julia backend internally, whatever the
# eventual fit will use, so the example is inert where Julia is absent --
# including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  gen <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 8))
  datalong <- ctGenerate(gen, n = 40, Tpoints = 8, backend = 'r')

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0)))

  result <- ctIdentify(datalong, model, cores = 1)
  print(result)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
