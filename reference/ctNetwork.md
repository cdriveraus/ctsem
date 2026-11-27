# Temporal and contemporaneous networks of a continuous time model

Computes the directed temporal network and the undirected
contemporaneous network used in network psychometrics, from a ctsem
model specification or from a fit on either backend.

In a continuous time model both networks depend on the time interval,
and their shape as well as their scale changes with `dt`.

## Usage

``` r
ctNetwork(
  x,
  dt = 1,
  state = NULL,
  standardise = TRUE,
  impulseType = "unit",
  threshold = 0,
  networks = c("temporal", "contemporaneous"),
  plot = FALSE,
  quiet = FALSE,
  ...,
  observational
)
```

## Arguments

- x:

  A `ctStanModel` from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
  or a fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) on
  either backend.

- dt:

  Time interval for the temporal and contemporaneous networks, a whole
  number of steps for a discrete time model. A single value. To compare
  intervals, call
  [`ctNetworkPlot`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
  with a vector.

- state:

  Where state dependent cells are evaluated, when `x` is a fit:
  `'T0MEANS'` (the default, as in
  [`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md)),
  `'mean'` for the mean smoothed latent state, `'asymptotic'` for the
  system's fixed point, or a numeric state vector. Values other than the
  default need `backend='julia'` and are an error for a stan fit. No
  effect for a linear model or an unfitted specification.

- standardise:

  If TRUE (the default), temporal edges are in standard deviation units
  of the processes, using `asymDIFFUSIONcov`, so edges between processes
  on different scales are comparable. Set to FALSE, with a message, when
  the system has no stationary variance.

- impulseType:

  What a one unit change in a process brings with it; see
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md).
  `'unit'` (the default) changes that process alone, giving the partial
  regression.

- threshold:

  Edges with `abs(weight)` at or below this are dropped from `$edges`
  and from any plot. The matrices are not thresholded.

- networks:

  Which networks appear in `$edges` and in the plot. Any of
  `'temporal'`, `'contemporaneous'`, `'asymptotic'`, `'measurement'`.
  All four matrices are always computed.

- plot:

  If TRUE, returns the
  [`ctNetworkPlot`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
  of the result instead of the result.

- quiet:

  Suppress the messages about how the numbers were computed.

- ...:

  Passed to
  [`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md)
  when `x` is a fit, and to
  [`ctNetworkPlot`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
  when `plot=TRUE`.

- observational:

  Deprecated. Use `impulseType`.

## Value

A list of class `ctNetwork`:

- `temporal`, `contemporaneous`, `asymptotic`, `measurement`:

  Weight matrices with the model's names, in the form `qgraph::qgraph()`
  takes. `asymptotic` may be `NULL`.

- `edges`:

  Data frame with columns `network`, `from`, `to`, `weight`, `directed`
  and `self`, thresholded and restricted to `networks`.

- `nodes`:

  Data frame of `name` and `type` (`'latent'` or `'manifest'`).

- `DRIFT`, `DIFFUSIONcov`, `asymDIFFUSIONcov`, `innovation`:

  The matrices the networks were computed from. On a julia fit `DRIFT`
  is the drift's Jacobian at the evaluation point, which equals DRIFT
  when no cell depends on the state.

With attributes `dt`, `standardise`, `impulseType`, `threshold`,
`networks`, `continuoustime`, `source`, `stateLabel` (the evaluation
point, for a fit) and `stateDependent` (the state dependent cells, or
`NULL`).

## Details

Three networks over the latent processes, plus the measurement model:

- `temporal`:

  `expm(DRIFT * dt)`, the discrete time autoregressive and cross-lagged
  matrix at interval `dt`. Directed, with the autoregressions as
  self-loops on the diagonal. Cell `[i,j]` is the effect of process `j`
  now on process `i` one interval later, so an edge runs from the column
  to the row. This is
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
  at a single interval, with the same `impulseType` and `standardise`
  arguments.

- `contemporaneous`:

  Partial correlations of the innovation covariance accumulated over
  `dt`,
  `Q(dt) = int_0^dt expm(DRIFT s) DIFFUSIONcov t(expm(DRIFT s)) ds`.
  This is the continuous time counterpart of the graphicalVAR
  contemporaneous network. As `dt` shrinks it tends to the partial
  correlations of `DIFFUSIONcov`, and as `dt` grows to those of the
  asymptotic covariance.

- `asymptotic`:

  Partial correlations of `asymDIFFUSIONcov`, the long run
  within-subject covariance. This is the network a cross-sectional study
  of a stationary system would find. `NULL` if the system is not
  stationary.

- `measurement`:

  `LAMBDA`, as a bipartite manifest by latent graph, directed from
  latent to manifest.

From a fit, the matrices are the posterior medians that
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md)
reports. From an unfitted specification, a fixed cell takes its value,
and a free cell takes the value its transform gives at a raw parameter
of zero. The number of DRIFT, DIFFUSION and LAMBDA cells set this way is
reported.

## State dependent and nonlinear models

When a DRIFT, DIFFUSION or LAMBDA cell is an expression of a latent
process, the network differs across the state space. The networks
reported for a fit are those of the system linearised at one point, and
elsewhere the edge weights differ, as can which edges are present. The
temporal and contemporaneous networks then use the drift's Jacobian at
that point (`JAx`), and the measurement network uses the measurement
Jacobian (`Jy`). These differ from DRIFT and LAMBDA where a cell is
nonlinear in the state: a cell `a + b*Y` in the column of process `Y`
has slope `a + 2*b*Y`. A stan fit has no Jacobians, so its edges are the
cells' values at `T0MEANS`, with a message saying so.

When such cells are present, `ctNetwork` gives a message naming the
evaluation point, and records the cells in `attr(x, 'stateDependent')`.
[`ctNetworkPlot`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
puts the point in the plot's subtitle. Use `state` to change the point,
and
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
or
[`ctStateDependencePlot`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)
to see how the system varies over the state space. For a linear model
there is no message.

From an unfitted specification such cells are reported and left out.

## See also

[`ctNetworkPlot`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
to draw it,
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
for the temporal network over a range of intervals,
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
and
[`ctStateDependencePlot`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)
for how a state dependent system varies away from the evaluation point.

## Examples

``` r
# Two processes, b driven by a, drawn at one interval.
m <- ctModel(type='ct', n.latent=2, n.manifest=2,
  manifestNames=c('y1','y2'), latentNames=c('a','b'),
  LAMBDA=diag(2),
  DRIFT=matrix(c(-1, 0, .5, -2), 2, 2, byrow=TRUE),
  DIFFUSION=matrix(c(1, 0, .3, 1), 2, 2, byrow=TRUE))

net <- ctNetwork(m, dt=1, standardise=FALSE)
round(net$temporal, 4)          # exp(-1) and exp(-2) on the diagonal
round(net$contemporaneous, 4)
net$edges
```
