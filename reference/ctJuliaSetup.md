# Configure the Julia engine used by ctsem

Prepares the Julia engine (ContinuousTimeSEM.jl) that ships with ctsem.
The engine's Julia package dependencies have to be installed, which
needs the network on a machine where they are not already present. As
with
[`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md),
ctsem asks before downloading in an interactive session, skips the step
with a message when the dependencies are already installed, and declines
in a non-interactive session unless `agree = TRUE` or the environment
variable `CTSEM_JULIA_AGREE=yes` is set.

## Usage

``` r
ctJuliaSetup(
  project = NULL,
  julia_bin = NULL,
  threads = NULL,
  force = FALSE,
  agree = NULL
)
```

## Arguments

- project:

  Optional local ContinuousTimeSEM.jl checkout to use instead of the
  copy that ships with ctsem. Rarely needed, mainly for development.

- julia_bin:

  Optional Julia binary directory.

- threads:

  Number of Julia threads. Julia fixes this when the session starts, so
  it takes effect only if no Julia session is running yet, or with
  `force = TRUE`, which restarts the session. `ctFit`'s `cores` cannot
  exceed it. `NULL`, the default, starts a session with
  `getOption("mc.cores", 2)` threads unless `JULIA_NUM_THREADS` is set,
  or with the `cores` of the call that starts the session if that is
  more. A number given here is also the most threads used by calls that
  do not set `cores`. A fit that asks for more cores than the session
  has threads says so (see
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)).

- force:

  Reconfigure an existing Julia session.

- agree:

  `TRUE` to consent to installing the engine's Julia package
  dependencies without being asked, `FALSE` to refuse. `NULL`, the
  default, asks in an interactive session and refuses otherwise (see
  [`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)).

## Value

A Julia-engine status list, invisibly.

## Details

The first call downloads and precompiles the dependencies (roughly 120
MB) and compiles the engine, which takes several minutes. Later sessions
reuse both, and the first fit of a new model in a session then takes
only seconds longer than later fits. Models with a categorical
indicator, or with more than four random effects per subject under the
Laplace approximation, compile a little more the first time.
`ctFit(forcerecompile = TRUE)` compiles a heavily state-dependent model
for faster gradients, at the cost of a minute or more of compiling.

Pressing Escape (or Ctrl-C) during a julia fit returns to the R prompt
at once, and keeps the Julia session and everything it has compiled.
Julia stops the interrupted work at its next iteration. Work without
iterations, such as compiling for a new model, finishes in the
background, and if it is still running when the next julia call starts,
that call says so and waits. Escape during that wait stops Julia
completely.
[`ctJuliaProcesses`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)
lists the Julia processes ctsem has running, including any left by R
sessions that have ended, and
[`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)
stops them.

If JuliaConnectoR or Julia itself is missing, this offers to install it,
as
[`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)
does. Use `ctJuliaInstall` to set up the backend deliberately, or from a
script.

## See also

[`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md),
[`ctJuliaStatus`](https://cdriveraus.github.io/ctsem/reference/ctJuliaStatus.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ctJuliaSetup()                         # prepare the engine, asking first
ctJuliaSetup(threads = 8, force = TRUE) # restart the session wider
} # }
```
