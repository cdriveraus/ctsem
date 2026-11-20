# ctsem

ctsem is an R package for continuous time structural equation modelling
of panel (N \> 1) and time series (N = 1) data, using either a
frequentist or Bayesian approach, or middle ground forms like maximum a
posteriori.

The general workflow begins by specifying a model using the
[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
function, then fitting it to data using
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) (alias
[`ctStanFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)).
For original OpenMx / SEM functionality from the first ctsem versions,
use the ctsemOMX package. For most purposes, the forms in ctsem are more
robust and flexible. For examples, see
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md). For
citation info, please run `citation('ctsem')` .

## References

https://www.jstatsoft.org/article/view/v077i05

Driver, C. C., & Voelkle, M. C. (2018). Hierarchical Bayesian continuous
time dynamic modeling. Psychological Methods. Advance online
publication.http://dx.doi.org/10.1037/met0000168

Stan Development Team (2018). RStan: the R interface to Stan. R package
version 2.17.3. http://mc-stan.org

\#' @keywords internal

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
