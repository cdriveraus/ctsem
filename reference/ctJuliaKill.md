# Stop Julia processes started for ctsem

Stops the Julia processes
[`ctJuliaProcesses`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)
lists. By default only those left by R sessions that have ended, which
nothing can be using.

## Usage

``` r
ctJuliaKill(
  which = c("orphaned", "session", "workers", "mine", "all"),
  pid = NULL,
  ask = interactive()
)
```

## Arguments

- which:

  Which to stop. `"orphaned"` (the default) stops processes whose R
  session has ended. `"session"` stops this R session's engine, which
  the next julia call replaces with a new one that has to compile again.
  `"workers"` stops this session's sampling workers, the background R
  sessions that sampling with `processes = TRUE` keeps for the next
  call, with their Julia processes. `"mine"` stops all three, and
  `"all"` every process listed, including the engines of other running R
  sessions, which then fail.

- pid:

  Process ids to stop instead of choosing by `which`. Only processes
  [`ctJuliaProcesses()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)
  lists can be stopped this way.

- ask:

  Whether to ask before stopping another R session's engine, as
  `which = "all"` or `pid` can.

## Value

The rows of
[`ctJuliaProcesses()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)
that were stopped, invisibly.

## See also

[`ctJuliaProcesses`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)

## Examples

``` r
# \donttest{
ctJuliaKill()             # engines left by R sessions that have ended
ctJuliaKill("session")    # this session's engine, e.g. stuck compiling
# }
```
