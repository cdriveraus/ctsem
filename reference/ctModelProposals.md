# Changes the score tests and adequacy checks propose

Lists the changes
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
would consider for a fit: the relaxations its score tests indicate and,
given an adequacy check, the structural changes the check's families
point to (a process of an item's own or shared by two items, a slow
process, a cycle, or second order dynamics). Rows can be passed to
`ctModelImprove(changes = )`, which refits with them and judges each by
the search's criterion, so that a search can be made one step at a time,
by hand.

## Usage

``` r
ctModelProposals(
  fit,
  score = NULL,
  adequacy = NULL,
  covariates = NULL,
  relax = c("heterogeneity", "covariates", "statedependence", "volatility", "fixed",
    "level", "slow", "unique", "periodic", "order"),
  alpha = 0.05,
  cores = getOption("mc.cores", 2L)
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
  `backend = 'julia'`.

- score:

  A
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
  of `fit`. Run when neither it nor `adequacy` is given.

- adequacy:

  Optionally a
  [`ctFitAdequacy`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
  of `fit`, for the structural changes.

- covariates, relax, alpha, cores:

  As for
  [`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md).

## Value

A data frame of class `ctModelProposals`, one row per change, strongest
first within the relaxations and within the structural changes: `change`
(as a search's `history` names it), `family`, `p` (the evidence ranked
by), `class`, and `proposal`, what
[`ctModelImprove()`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
needs to make the change.

## Details

Every indicated change is listed, including relaxations the tests cannot
tell apart (the same `class`), of which a search round refits only the
strongest. The structural changes are listed beside the relaxations,
while a search round tries them only once the screen finds nothing.

## See also

[`vignette('score-tests')`](https://cdriveraus.github.io/ctsem/articles/score-tests.md).
