# Set up the Julia backend

Installs everything `ctFit(backend='julia')` needs, in one call and
without restarting R: the JuliaConnectoR package, Julia, and the Julia
packages that the engine shipped with ctsem depends on. Steps that are
already done are skipped, so on a machine that already has Julia this
only prepares the engine.

## Usage

``` r
ctJuliaInstall(
  threads = NULL,
  version = NULL,
  agree = NULL,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- threads:

  Number of Julia threads. Julia fixes this when it starts, so it
  applies to the session this call starts (see
  [`ctJuliaSetup`](https://cdriveraus.github.io/ctsem/reference/ctJuliaSetup.md)).

- version:

  Julia version to install. Defaults to the version the engine shipped
  with ctsem was tested with. The download is checksum-verified only for
  that default.

- agree:

  `TRUE` to consent to the installs without being asked, `FALSE` to
  refuse them. `NULL`, the default, asks in an interactive session and
  refuses otherwise.

- force:

  Reinstall Julia even if one is already available, and restart any
  running Julia session.

- quiet:

  Suppress download and installation progress output.

## Value

A Julia-engine status list, invisibly; see
[`ctJuliaStatus`](https://cdriveraus.github.io/ctsem/reference/ctJuliaStatus.md).

## Details

Each step that installs something asks first, naming the download and
where it goes. In a non-interactive session nothing is installed unless
`agree = TRUE` or the environment variable `CTSEM_JULIA_AGREE=yes` is
set.

Julia is downloaded only if none can be found. ctsem looks in
`JULIA_BINDIR`, then at a Julia it installed itself, then on the `PATH`,
then among juliaup's installations. A Julia that ctsem installs goes
under `tools::R_user_dir("ctsem", "data")`, and later sessions find it
without any environment variable. The `PATH` and shell profile are left
unchanged. To use a different Julia afterwards, set `JULIA_BINDIR`.

## See also

[`ctJuliaSetup`](https://cdriveraus.github.io/ctsem/reference/ctJuliaSetup.md),
[`ctJuliaStatus`](https://cdriveraus.github.io/ctsem/reference/ctJuliaStatus.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ctJuliaInstall()            # ask, then install whatever is missing
ctJuliaInstall(threads = 4) # and start the Julia session with four threads
} # }
```
