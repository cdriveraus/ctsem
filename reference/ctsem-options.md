# Options ctsem reads

Settings that change how ctsem behaves for a whole session. Set them
with [`options()`](https://rdrr.io/r/base/options.html), usually once at
the top of a script or when a front end starts.

## Value

Nothing. This page documents options.

## Progress and output

- `ctsem.progress.overwrite`:

  Whether progress reporting rewrites a single line (`TRUE`) or prints
  separate lines now and then (`FALSE`). If unset, ctsem decides from
  where the output goes. In a Shiny session, under
  [`sink()`](https://rdrr.io/r/base/sink.html) or in a knitr chunk the
  line is not rewritten. Set it if that choice is wrong for your front
  end.

- `ctsem.cluster.outfile`:

  Where parallel workers send their output. If unset, they are silent.
  `""` sends their output to the console, which is useful when debugging
  a worker.

## Per-call alternatives

Progress and reporting can also be set for each fit, and these settings
take precedence: `verbose` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) (0
quiet, 1 reports progress, 2 keeps every progress line),
`optimcontrol$callback` for a function called while the fit runs, and
`optimcontrol$progress` to turn progress on or off for one fit.

## Julia

The Julia backend is mostly configured through environment variables:
`JULIA_BINDIR` names the Julia binary directory, and `JULIA_NUM_THREADS`
the number of threads of the next session that starts.
`ctJuliaSetup(project=)` points the session at a different engine. See
[`ctJuliaSetup`](https://cdriveraus.github.io/ctsem/reference/ctJuliaSetup.md).

- `ctsem.julia.restart`:

  Whether a fit asking for more `cores` than the running Julia session
  has threads may restart the session to get them. `FALSE` by default,
  because Julia fixes its number of threads when it starts, and a
  restart loses the code it has compiled. With `FALSE`, such a fit runs
  with the threads available and says so.
