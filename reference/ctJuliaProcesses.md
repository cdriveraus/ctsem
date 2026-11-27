# Julia processes started for ctsem

Lists the Julia processes on this machine that ctsem started: this R
session's engine, the engines of its sampling workers, those of other R
sessions, and any left running by an R session that has ended. A Julia
process can keep running after it is no longer needed, for example after
Escape interrupts a step that cannot be interrupted, or after R crashes.
Look here when the machine is busier than it should be.
[`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)
stops them.

## Usage

``` r
ctJuliaProcesses(interval = 0.5)
```

## Arguments

- interval:

  Seconds over which to measure CPU use, reported in the `busy` column.
  `0` skips the measurement.

## Value

A data frame with one row per process, and columns:

- `pid`:

  The process id.

- `role`:

  `"this session"`, `"worker of this session"`, `"another R session"`,
  `"orphaned"` (its R session has ended),
  `"JuliaConnectoR, unrecorded"`, or `"child of <pid>"`.

- `busy`, `cpu_percent`:

  CPU use over `interval`.

- `cpu_seconds`, `memory_mb`:

  CPU time and memory used.

- `julia_threads`:

  The number of Julia threads ctsem started it with.

- `os_threads`:

  The number of operating system threads.

- `r_pid`:

  The R process it serves.

- `started`:

  When it started.

- `note`:

  Says when this session's engine is still finishing an interrupted
  call.

## Details

ctsem records each Julia session it starts. Unrecorded JuliaConnectoR
sessions, such as one started by other code or by an older ctsem, are
listed too, as are Julia processes they started (usually to precompile
the engine). CPU, memory, thread and parent details need the ps package.
Without it, only the processes and their owners are listed.

## See also

[`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md),
[`ctJuliaStatus`](https://cdriveraus.github.io/ctsem/reference/ctJuliaStatus.md)

## Examples

``` r
# \donttest{
ctJuliaProcesses()
# }
```
