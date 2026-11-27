# ctsem

ctsem fits continuous and discrete time dynamic models to longitudinal
data: panel data, intensive longitudinal data or a single time series,
with observations at any times. The processes can interact, be of higher
order or nonlinear, and be measured by continuous, binary, ordinal,
count or censored variables. Any parameter can differ between subjects,
at one or more levels, and depend on covariates. Estimation is by
maximum likelihood, or Bayesian by optimisation or sampling, with an
engine written in Julia (the default) or with Stan.

## Details

A model is specified with
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md), or
in lavaan, Mplus or formula syntax with
[`ctModelLavaan`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
[`ctModelMplus`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
and
[`ctModelFormula`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
and fitted with
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).
[`summary()`](https://rdrr.io/r/base/summary.html),
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
and
[`ctSubjectPars`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md)
describe a fit, and
[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
[`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
and
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
check and improve it.
[`ctJuliaInstall`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)
sets up the Julia engine.

[`vignette('ctsem')`](https://cdriveraus.github.io/ctsem/articles/ctsem.md)
works through a first analysis, and the website,
<https://cdriveraus.github.io/ctsem/>, has the vignettes and the help
pages grouped by task. For the original OpenMx based ctsem, see the
ctsemOMX package. `citation('ctsem')` gives the references to cite.

## References

Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous time
structural equation modeling with R package ctsem. Journal of
Statistical Software, 77(5), 1-35.
[doi:10.18637/jss.v077.i05](https://doi.org/10.18637/jss.v077.i05)

Driver, C. C., & Voelkle, M. C. (2018). Hierarchical Bayesian continuous
time dynamic modeling. Psychological Methods, 23(4), 774-799.
[doi:10.1037/met0000168](https://doi.org/10.1037/met0000168)

## See also

Useful links:

- <https://cdriveraus.github.io/ctsem/>

- <https://github.com/cdriveraus/ctsem>

- Report bugs at <https://github.com/cdriveraus/ctsem/issues>

## Author

**Maintainer**: Charles Driver <charles.driver2@uzh.ch> \[copyright
holder\]

Authors:

- Charles Driver <charles.driver2@uzh.ch> \[copyright holder\]

- Manuel Voelkle \[copyright holder\]

- Han Oud \[copyright holder\]

Other contributors:

- Trustees of Columbia University \[copyright holder\]
