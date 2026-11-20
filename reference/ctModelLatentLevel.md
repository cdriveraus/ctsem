# Put subjects' levels in the latent processes instead of the indicators

In the usual model a subject's level is in random manifest intercepts
(`MANIFESTMEANS`), and the latent state is the subject's deviation from
it. For linear dynamics this is the same model as one where the latent
process holds the level itself, through a random latent intercept
(`CINT`). The two differ once the dynamics depend on the state, since a
term that uses the state then uses the deviation in the first model and
the level in the second. Persistence that differs near the bottom of a
scale, or innovations that grow with a subject's usual level, depend on
the level, and can only be expressed in the second model.
`ctScoreTest(tests = 'level')` tests for such dependence between
subjects.

## Usage

``` r
ctModelLatentLevel(ctm, latents = NULL)
```

## Arguments

- ctm:

  A model from `ctModel(type = 'ct')`.

- latents:

  Which latents, by name or position; by default every one.

## Value

The model, with the latents whose level moved in attribute `moved`.

## Details

The function changes each latent in `latents` that has a marker (an
indicator with loading fixed at 1 on that latent and 0 on all others)
whose `MANIFESTMEANS` intercept varies over subjects. `CINT` becomes a
free parameter, named `cint_` and the latent's name, with the
intercept's random effect and predictor effects. The marker's intercept
is fixed at 0, and the latent's initial mean is set to its stationary
value (`T0MEANS` `'stationary'`), so each subject starts at their own
level. The intercepts of other indicators are left as they are. A latent
without such a marker, or whose `CINT` is already free, is left
unchanged, with a message.

## See also

[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
for the tests of level dependence.

## Examples

``` r
m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'), LAMBDA = diag(2), CINT = matrix(0, 2),
  MANIFESTMEANS = matrix(c('mm1', 'mm2'), 2))
lm <- ctModelLatentLevel(m)
attr(lm, 'moved')
p <- lm$pars
p[p$matrix %in% c('CINT', 'MANIFESTMEANS', 'T0MEANS'), c('matrix', 'row', 'param',
  'value', 'indvarying')]
```
