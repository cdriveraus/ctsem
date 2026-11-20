# Summarise empirical Bayes subject-wise ctsem fits

Summarise empirical Bayes subject-wise ctsem fits

## Usage

``` r
# S3 method for class 'ctEmpiricalBayesFit'
summary(
  object,
  use = c("rawest", "rawposterior"),
  probs = c(0.025, 0.5, 0.975),
  sdscale = c("unit", "rawsd"),
  minsd = 1e-06,
  digits = 4,
  ...
)
```

## Arguments

- object:

  Object returned by
  [`ctEmpiricalBayesFit`](https://cdriveraus.github.io/ctsem/reference/ctEmpiricalBayesFit.md).

- use:

  `'rawest'` to summarise final-pass subject point estimates, or
  `'rawposterior'` to pool final-pass subject raw posterior samples.

- probs:

  Quantiles to report for transformed parameters.

- sdscale:

  How to set `model$pars$sdscale` when reconstructing the adjusted
  single-subject empirical Bayes model from the final empirical raw
  distribution. `'unit'` keeps any later random-effect SDs on the
  EB-standardised raw scale. `'rawsd'` uses the final empirical raw SDs
  directly. Retained for compatibility;
  [`summary()`](https://rdrr.io/r/base/summary.html) no longer returns a
  reconstructed model.

- minsd:

  Lower bound used for empirical raw SDs before model adjustment.

- digits:

  Number of digits for printed summary tables.

- ...:

  Unused.

## Value

Compact list containing fit settings, transformed-parameter summaries,
outlier diagnostics, transformed-parameter covariance / correlation
matrices, and `subjectpars`: each subject's final-pass point estimates
on the parameters' own scale, one row per subject with its id in
`subject`, unrounded, for comparisons across subjects. The subject fit
lists, raw estimate matrices, pass maps, and adjusted model remain on
the original `ctEmpiricalBayesFit` object.
