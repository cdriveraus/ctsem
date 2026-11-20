# A comparison between a fit's own answer and a more accurate one

The common shape behind \[ctLaplaceCheck()\], \[ctParticleLik()\] and
the curvature certification on \`fit\$uncertainty\$certification\`.

## Usage

``` r
.ctFitGap(
  reference,
  gap,
  tolerance = NA_real_,
  gap_se = NA_real_,
  nsubjects = NA_integer_,
  remedy = NULL,
  detail = NULL
)
```

## Arguments

- reference:

  short name for what the fit was compared against: \`"curvature"\`,
  \`"quadrature"\` or \`"particle"\`.

- gap:

  log likelihood difference, better minus reported. Positive means the
  better computation finds more likelihood than the fit reports.

  One asymmetry between the three, which is a property of the references
  and not of this shape: the curvature gap is \`1/2 g' H^-1 g\` and so
  is non-negative by construction, while the quadrature and particle
  gaps are signed. A negative particle gap says the assumed-density
  filter reports \*more\* likelihood than the particle filter finds,
  which is a finding rather than an error – hence \`resolved\` tests
  \`abs(gap)\`, not \`gap\`.

- tolerance:

  the bar, in log likelihood units, or \`NA\` where the route states
  none.

- gap_se:

  standard error of \`gap\`, where the reference is stochastic.

- nsubjects:

  used to report the gap per subject, which is what says whether a total
  is large.

- remedy:

  one sentence on what to do if the gap matters.

- detail:

  optional per-parameter or per-row table the route already produces.

## Value

an object of class \`ctFitGap\`.
