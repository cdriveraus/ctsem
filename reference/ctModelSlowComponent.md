# Add slow processes that drive a model's latent processes

Each latent in `indices` gets a slow latent, named `slow_` and the
latent's name, that drives it through the drift. An estimated effect
`DRIFT[k, slow_k]` shifts the level the latent returns to with the slow
process, so a process can drift over weeks while fluctuating over hours.
The latent's equilibrium is `DRIFT[k, slow_k] / |DRIFT[k, k]|` times the
slow process, so a person with faster dynamics follows the slow process
less closely.

## Usage

``` r
ctModelSlowComponent(ctm, indices = seq_len(ctm$n.latent))
```

## Arguments

- ctm:

  A model from `ctModel(type = 'ct')`.

- indices:

  Which latents get a slow process, by position.

## Value

The extended model.

## Details

The slow process is always slower than the latent it drives. Its drift
is the latent's drift times `inv_logit(r)`, with `r` a `PARS` parameter
named `drift_` and the slow process's name. If both rates were free,
they could swap, since with only the latent measured the two processes
are interchangeable. The latent's drift parameter is declared in `PARS`
for this, keeping its name, transform, individual variation and
predictor effects, so where the latent's rate varies over people the
slow process's rate varies with it. Where the latent's drift is fixed,
the slow process's drift is a fraction of that value. Where the latent's
drift is an expression, the slow process's persistence is left free.

A slow process has innovations of unit variance (the effect carries its
scale), correlated with those of other slow processes, including ones
added earlier, and uncorrelated with the others. It has no loadings, no
intercept, and no individual variation or predictor effects. Its initial
mean and variance are its stationary ones (`T0MEANS` and the `T0VAR`
diagonal `'stationary'`, see
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)), zero
and `1 / (2 |a|)` for persistence `a`. Its initial covariance with the
latent it drives is free, and with other states zero.

The effect needs a starting value away from zero, because at zero
neither the slow process nor its effect has a gradient. When refitting
from a model without the slow process, start the effect at one (either
sign) and the fraction small.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
starts it at a tenth.

## See also

[`ctModelHigherOrder`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md)
for the other structural extension the search tries.

## Examples

``` r
m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'), LAMBDA = diag(2))
m2 <- ctModelSlowComponent(m, 1:2)
m2$latentNames
p <- m2$pars
p[p$matrix %in% c('DRIFT', 'DIFFUSION') & (p$row > 2 | p$col > 2) &
  (!is.na(p$param) | p$value != 0), c('matrix', 'row', 'col', 'param', 'value')]
```
