# Report Julia backend availability

Reports what `backend='julia'` would find, without installing anything
or starting a Julia session. It is safe to call on a machine without
Julia or JuliaConnectoR.
[`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)
installs anything it reports as missing.

## Usage

``` r
ctJuliaStatus(project = NULL, julia_bin = NULL)
```

## Arguments

- project:

  Optional local ContinuousTimeSEM.jl checkout to use instead of the
  copy that ships with ctsem. Rarely needed, mainly for development.

- julia_bin:

  Optional Julia binary directory.

## Value

A list describing the Julia engine:

- `available`:

  TRUE if a julia fit could run now without asking for anything: Julia
  starts and the engine's Julia packages are installed, or consent to
  install them was given with `CTSEM_JULIA_AGREE=yes`.

- `connectoR`:

  Whether JuliaConnectoR is installed at a suitable version.

- `julia_bin`, `julia`:

  The Julia directory in use, and the Julia version there.

- `project`:

  The engine's project directory.

- `engine`:

  The version of the engine shipped with this ctsem.

- `threads`:

  The number of threads of a running Julia session, or `NA` if none is
  running.

## Examples

``` r
status <- ctJuliaStatus()
status$available
```
