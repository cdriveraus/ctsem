# Phase portrait of a fitted continuous time model

Draws the model's vector field over a plane of the latent state space,
showing at each point the direction and speed in which the deterministic
part of the system moves. For a model whose matrices depend on the
latent state, this shows behaviour that a single DRIFT matrix cannot
describe.

## Usage

``` r
ctPhasePortrait(
  fit,
  latents = 1:2,
  extent = "data",
  gridsize = 15,
  state = "asymptotic",
  tipreds = NULL,
  nullclines = TRUE,
  fixedpoint = TRUE,
  trajectories = 0,
  plot = TRUE
)
```

## Arguments

- fit:

  A ctJuliaFit, or a ctStanFit of a linear model.

- latents:

  Two latent processes, by name or index, forming the plane.

- extent:

  `'data'` (the default) spans the smoothed latent states, `'sd'` spans
  two stationary standard deviations either side of T0MEANS, or supply a
  list of two numeric ranges.

- gridsize:

  Arrows per side. Computing time grows with its square.

- state:

  Where the \*other\* latent processes are held, as in
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md):
  `'asymptotic'`, `'mean'`, `'T0MEANS'`, or a numeric state.

- tipreds:

  Time independent predictor values, or NULL for the population.

- nullclines:

  Draw the curves where each process stops changing?

- fixedpoint:

  Mark the asymptotic state when it is inside the region?

- trajectories:

  Number of subjects' smoothed trajectories to overlay, or 0 for none.

- plot:

  If FALSE, return the field, nullclines and trajectories as data
  instead of a plot.

## Value

A ggplot, or when `plot=FALSE` a list with `field`, the `latents` on the
plane, `stateLabel` naming where the off-plane processes were held, the
`held` process names, and `nullclines`, `fixedpoint` and `trajectories`
when asked for. The plot names `stateLabel` in its subtitle whenever
anything is held off-plane.

## Details

Three things are drawn over the field, each of which can be turned off:

- **Nullclines**, the curves where one process stops changing. Where the
  nullclines of two processes cross, the whole system is at rest. That
  is a fixed point, and a nonlinear system may have several.

- **The fixed point** found by `state='asymptotic'`, marked when it lies
  inside the plotted region. Only one fixed point is found, starting
  from T0MEANS. The nullcline crossings show whether there are others.

- **Smoothed trajectories** for a few subjects, so the field can be
  compared with where the data went.

With more than two latent processes the others are held at `state`, so
the plot shows one slice of the state space, which changes with `state`.

A linear model has an affine field, straight nullclines and exactly one
fixed point. It can be drawn from a stan or a julia fit, but shows
nothing beyond the DRIFT matrix. A nonlinear stan fit is refused.

## See also

[`ctContextDependence`](https://cdriveraus.github.io/ctsem/reference/ctContextDependence.md),
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)

## Examples

``` r
# The arrows show the direction and speed in which the two processes of the
# bundled fit move from each point of the plane. The lines are the
# nullclines, where one process stops changing, and they cross at the fixed
# point. The smoothed trajectories of three subjects are drawn over it.
ctPhasePortrait(ctstantestfit, gridsize = 9, trajectories = 3)

# The same as data: one row of the field per arrow, with its position
# (x, y), direction (dx, dy) and speed, and the fixed point.
pp <- ctPhasePortrait(ctstantestfit, gridsize = 7, plot = FALSE)
head(pp$field)
pp$fixedpoint
```
