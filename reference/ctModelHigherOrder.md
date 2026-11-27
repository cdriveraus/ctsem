# Raise the order of a ctsem model's latent processes.

Each latent in `indices` gets a derivative latent (named `d` and the
latent's name), which drives the level with a fixed coefficient of one.
Each becomes a damped oscillator, a second order process whose
autocorrelation can fall below zero and rise again.

## Usage

``` r
ctModelHigherOrder(
  ctm,
  indices,
  diffusion = TRUE,
  crosseffects = FALSE,
  cint = FALSE,
  explosive = FALSE
)
```

## Arguments

- ctm:

  ctModel, `type = 'ct'` or `type = 'omx'`.

- indices:

  Vector of integers, which latents to raise the order of.

- diffusion:

  Move the process noise from the level to its derivative?

- crosseffects:

  Move the effects of the other processes on this one to its derivative
  (divided by \\\tau\\, for a `type = 'ct'` model)? For a `type = 'omx'`
  model, its effects on them move as well.

- cint:

  Move the continuous intercept to the derivative (divided by \\\tau\\,
  for a `type = 'ct'` model)?

- explosive:

  For a `type = 'omx'` model, allow explosive (non equilibrium
  returning) processes. A `type = 'ct'` model's level keeps its own
  drift's transform, which decides that.

## Value

extended ctModel

## Details

A model from `ctModel(type = 'ct')` is rewritten on its parameter table,
with a new inertia time \\\tau \> 0\\ for each latent (`inertia_` and
the latent's name, in `PARS`): \$\$\tau \ddot{x} + \dot{x} = a x + q \\
dW,\$\$ where \\a\\ is the level's own drift and \\q\\ its diffusion,
both keeping their names and settings (transform, individual variation,
predictor effects). The derivative's row of `DRIFT` is \\a/\tau\\ and
\\-1/\tau\\, and its diffusion is \\q/\tau\\. As \\\tau\\ goes to zero
the process becomes the original first order one, so the new model
contains the old one with every old parameter at its value. For \\\tau
\ge 1/(4\|a\|)\\ the process oscillates. The derivative starts at its
equilibrium given the process: mean zero, uncorrelated with other
states, and variance \\q^2/(2\tau)\\, so the inertia time is the only
new parameter. With the noise left on the level (`diffusion = FALSE`),
or no noise, the derivative's initial mean, variance and covariances are
free.

A model in the older matrix-list form from `ctModel(type = 'omx')` is
extended differently, with the derivative pulled back by the level
through an estimated coefficient (negative unless `explosive`), and the
level's own autoregression as the derivative's damping. This gives the
same family of processes, but the first order model is only reached as
the damping, the coefficient and the diffusion all go to infinity.

## Examples

``` r
om <- ctModel(LAMBDA=diag(1,2),DRIFT=0, 
  MANIFESTMEANS=0,type='omx',Tpoints=4)
  
om <- ctModelHigherOrder(om,1:2)
print(om$DRIFT)

print(om$pars)
```
