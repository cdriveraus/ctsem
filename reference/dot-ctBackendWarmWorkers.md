# Start worker sessions and get them ready to sample

Launches \`workers\` background R sessions and sets each one evaluating
one gradient for \`object\`, which is what forces Julia to build the
objective and compile the engine for this model's shape. Returns
immediately: the sessions work while the caller carries on.

## Usage

``` r
.ctBackendWarmWorkers(
  object,
  workers,
  values = NULL,
  threads = 1L,
  layer = FALSE
)
```

## Arguments

- object:

  A prepared \`ctJuliaModel\` or a \`ctJuliaFit\`.

- workers:

  How many sessions to start.

- values:

  Parameter vector to evaluate at. Any finite point compiles the same
  code; a fit's own estimate is the natural choice when there is one.

- threads:

  Julia threads each worker's session starts with, one for all or one
  per worker (\`.ctBackendWorkerThreads()\`). Fixed when the session
  starts, so a pool warmed at other widths is replaced rather than
  reused.

- layer:

  Whether the work the workers will do evaluates the layer
  (\`.ctJuliaObjective(layer = TRUE)\`, SAEM's kernel) rather than the
  fit's objective, so each builds and compiles the objective its chains
  use.

## Value

A list of future handles, or \`NULL\` when warming is unavailable.

## Details

The sessions persist. \`future::multisession\` keeps its pool across
calls, so the warmed sessions are the same ones a later chain runs in –
which is the point, and the reason this sets a plan rather than making a
one-shot cluster.
