# Fold the deprecated sampling arguments into `sampleControl`

[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) took
`iter`, `chains` and `control` as arguments of its own, which put three
of a sampler's settings in the signature and the rest in a list – and
made `control` mean rstan's control list on one backend and the julia
sampler's settings on the other. They are all entries of `sampleControl`
now.

## Usage

``` r
.ctSampleControlResolve(
  sampleControl = list(),
  given = character(0),
  iter = 1000L,
  chains = 2L,
  control = list()
)
```

## Arguments

- sampleControl:

  The new list.

- given:

  Names the caller actually used, from `names(match.call())` in the
  calling function – a default cannot be told from a value that happens
  to equal it any other way.

- iter, chains, control:

  The deprecated arguments, as they arrived.

## Value

A list of `iter`, `chains` and `control` to carry on with.

## Details

Resolved here rather than threaded through: everything downstream, on
both backends, still receives `iter`, `chains` and `control` as it
always did, so the rename cannot change what a fit does. In particular
`control` still reaches `stan` unchanged on the stan path.

Precedence is `sampleControl` first, because it is the argument being
kept. A deprecated argument that was \*explicitly supplied\* is used
only for what `sampleControl` does not say, and is reported either way –
a silent deprecation teaches nobody, and a silently ignored one is
worse.
