# One line describing every chain

The whole report on one line, so that it can be overwritten in place: a
carriage return goes to the start of the last visual row, so a block of
\`chains\` lines cannot be updated without cursor movement that not
every console honours. One line can, and the same line is what the
single-process path prints.

## Usage

``` r
.ctBackendProcessLine(infos, chains, elapsed, width = getOption("width", 80L))
```

## Arguments

- infos:

  One entry per chain from \[.ctBackendReadProgressFile()\], with
  \`NULL\` for a chain that has not written yet.

- chains:

  Number of chains, so that the line can say when it is describing fewer
  of them than are running.

- elapsed:

  Seconds each chain has spent in its current phase, for the rate the
  estimate extrapolates from.

- width:

  Console width to truncate to. A line that wraps cannot be overwritten
  in place – the earlier rows are left behind as debris.

## Value

One string, or \`NULL\` when no chain has reported yet.

## Details

What is on it, and why in this order. The phase and the per-chain
iteration counts answer "is this going to finish"; the counts share
their total and their phase label whenever the chains agree on them,
which is nearly always and is what keeps the line inside a terminal
width. Then the time remaining, from the \*slowest\* chain, because that
is when the run ends. Then the log posterior per chain, which is what
says whether the draws are going anywhere sensible – and it is per chain
rather than summarised because a single chain stuck in a bad region is
the failure this reveals, and an average hides exactly that. Divergences
last, and only once there are some.
