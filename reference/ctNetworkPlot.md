# Draw the temporal and contemporaneous networks of a continuous time model

Plots the networks
[`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
computes, in the usual style of network psychometrics: signed edges with
width by magnitude, arrowheads and self-loops on the directed temporal
network, and no arrowheads on the undirected contemporaneous network.

Given more than one `dt`, the networks are drawn side by side at each
interval, showing how the network a study finds depends on how often it
measures.

## Usage

``` r
ctNetworkPlot(
  x,
  dt = 1,
  state = NULL,
  networks = c("temporal", "contemporaneous"),
  threshold = 0,
  engine = c("ggplot", "qgraph"),
  maxwidth = 3,
  arrowsize = 2.6,
  nodesize = 10,
  labelsize = 3.2,
  poscolour = "#2166AC",
  negcolour = "#B2182B",
  title = "auto",
  ...
)
```

## Arguments

- x:

  A `ctStanModel`, a fit, or a `ctNetwork` from
  [`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md).

- dt:

  One or more time intervals. Ignored, with a message, when `x` is
  already a `ctNetwork`.

- state:

  Where state dependent cells are evaluated; see
  [`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md).
  Ignored, with a message, when `x` is already a `ctNetwork`. For a
  state dependent system the point appears in the plot's subtitle. A
  linear system gets no subtitle.

- networks:

  Which networks to draw; see
  [`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md).
  Each network and each `dt` gets its own panel.

- threshold:

  Edges with `abs(weight)` at or below this are not drawn.

- engine:

  `'ggplot'` (the default) or `'qgraph'`.

- maxwidth:

  Line width of the largest edge drawn.

- arrowsize:

  Arrowhead length in mm.

- nodesize, labelsize:

  Node and label size.

- poscolour, negcolour:

  Colours for positive and negative edges. The defaults are a
  colour-blind safe blue and red. For the common green and red, use
  `'darkgreen'` and `'red'`.

- title:

  Panel title, or NULL for none. `'auto'` names the source and the
  interval.

- ...:

  Passed to
  [`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md),
  or to `qgraph::qgraph()` when `engine='qgraph'`.

## Value

A ggplot2 object, which can be modified further or printed. With
`engine='qgraph'`, the qgraph object, invisibly, after drawing it.

## Details

Drawn with ggplot2. `engine='qgraph'` passes the weight matrix to
`qgraph::qgraph()` instead, if that package is installed. The matrices
from
[`ctNetwork`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
can also be passed to qgraph directly, as in
`qgraph::qgraph(ctNetwork(fit)$temporal, directed=TRUE)`.

## Examples

``` r
m <- ctModel(type='ct', n.latent=3, n.manifest=3,
  manifestNames=c('y1','y2','y3'), latentNames=c('a','b','c'),
  LAMBDA=diag(3),
  DRIFT=matrix(c(-.4,0,0, .3,-.5,0, 0,.4,-.6), 3, 3, byrow=TRUE),
  DIFFUSION=matrix(c(1,0,0, .3,1,0, 0,-.4,1), 3, 3, byrow=TRUE))

# The indirect path a -> b -> c is invisible at a short interval and present
# at a long one -- the same model, two networks.
ctNetworkPlot(m, dt=c(.2, 2), networks='temporal')
```
