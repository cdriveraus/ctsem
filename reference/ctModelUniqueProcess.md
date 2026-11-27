# Give indicators processes of their own

Each manifest variable in `items` gets a latent, named `u_` and the
variable's name, that loads only on that variable with a fixed loading
of one, and has its own persistence and innovations. The variable then
has, besides what it shares with others through the model's latents and
its measurement error, a specific part that persists over time. This
suits an item whose wording picks up something besides the construct,
consistently over time.

## Usage

``` r
ctModelUniqueProcess(ctm, items, shared = FALSE, name = NULL)
```

## Arguments

- ctm:

  A model from `ctModel(type = 'ct')`.

- items:

  Which manifest variables, by name or position.

- shared:

  If TRUE, one process for all of `items`, named `u_` and their names
  joined by `_`, instead of one each.

- name:

  Names for the new latents, replacing the defaults: one for a shared
  process, otherwise one per item. Useful for adding a second process on
  the same items at another time scale.

## Value

The extended model.

## Details

Without such a process, the specific part can only go into the
measurement error, which is independent over time, or into the shared
latents, whose dynamics it then distorts. In the first case the item's
prediction errors stay correlated over nearby times, which the
`residuals` family of
[`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
shows. The specific process and the measurement error are distinguished
by persistence over time.

The process has no intercept (the item's own intercept carries its
level), ctsem's default drift transform, so it is mean reverting, and
innovations independent of the other latents. It starts from its
stationary distribution, with mean zero and the variance its persistence
and innovations imply (the `T0VAR` diagonal `'stationary'`, see
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)).

With `shared = TRUE` the items share one process, with a loading of one
on the first item and free loadings on the others, representing what
moves them together beyond the model's other latents. Measurement errors
that are correlated within an occasion are often better modelled this
way, as a process faster than the occasions, unless the design gives a
reason for correlated error, such as a method shared by the items.

## See also

[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md),
which adds one for the items the adequacy check points to.

## Examples

``` r
m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2', 'Y3'),
  latentNames = 'eta', LAMBDA = matrix(c(1, 'l2', 'l3'), 3))
m2 <- ctModelUniqueProcess(m, 'Y2')
m2$latentNames
p <- m2$pars
p[p$col %in% 2 & p$matrix %in% c('LAMBDA', 'DRIFT', 'DIFFUSION'),
  c('matrix', 'row', 'col', 'param', 'value')]
ctModelUniqueProcess(m, c('Y2', 'Y3'), shared = TRUE)$latentNames
```
