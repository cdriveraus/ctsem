# `.ctMatsetupFreeRows()` replaced eleven inline spellings of "which matsetup row
# is a free population parameter", across five files, which did not agree with
# each other. These tests pin each replaced predicate to the expression it
# replaced, over an exhaustive grid of the five columns the predicates read --
# rather than over whatever combinations one example model happens to produce,
# which is how the variants got to disagree unnoticed in the first place.
#
# Cheap and always runs: no fit, no julia, no data.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

# Every combination of the columns the predicates discriminate on. `when`
# carries the real enum including the 100 wildcard and the -999 never-used
# sentinel; see R/ctMatsetup.R for what each value means.
.ms_grid <- function() {
  g <- expand.grid(
    when = c(-999, -1, 0, 1, 2, 3, 4, 100),
    param = c(0L, 1L, 7L),
    copyrow = c(0L, 1L, 3L),
    tipred = c(0L, 1L),
    indvarying = c(0L, 2L),
    KEEP.OUT.ATTRS = FALSE)
  g$row <- 1L
  g$col <- 1L
  g$matrix <- 5L
  g$parname <- "p"
  g
}

test_that("the grid actually exercises every case the predicates split on", {
  ms <- .ms_grid()
  expect_equal(nrow(ms), 8 * 3 * 3 * 2 * 2)
  expect_true(any(ms$when == 100))
  expect_true(any(ms$when == -1))
  expect_true(any(ms$copyrow > 0))
  expect_true(any(ms$tipred > 0 & ms$indvarying == 0))
  expect_true(any(ms$indvarying > 0 & ms$tipred == 0))
})

test_that("free rows are the population rows carrying a parameter", {
  ms <- .ms_grid()
  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms),
    ms$when %in% c(0, -1) & ms$param > 0)
})

test_that("free = FALSE drops the parameter requirement, for the numbering range", {
  # R/ctData.R's `nparams` takes max() over the numbering, where a fixed row
  # still carries a valid number.
  ms <- .ms_grid()
  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms, free = FALSE),
    ms$when %in% c(0, -1))
})

test_that("defining = TRUE excludes rows that copy another row's parameter", {
  ms <- .ms_grid()
  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms, defining = TRUE),
    ms$when %in% c(0, -1) & ms$param > 0 & ms$copyrow < 1)
})

test_that("varying = TRUE keeps parameters that differ over subjects", {
  ms <- .ms_grid()
  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms, defining = TRUE, varying = TRUE),
    ms$when %in% c(0, -1) & ms$param > 0 & ms$copyrow < 1 &
      (ms$tipred > 0 | ms$indvarying > 0))
})

test_that("the PARS wildcard is reached only by asking for it by name", {
  # R/ctData.R's whenvecp is the one caller that wants `when == 100`, because
  # PARS are needed at every when. Every other caller must not get those rows:
  # a PARS carrier row answering a parameter lookup is what cost a drift
  # diagonal its negative definiteness once already -- see the post-mortem in
  # R/ctModelWriter.R above the whenvecp block.
  ms <- .ms_grid()

  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms, when = c(0, 100), defining = TRUE, varying = TRUE),
    ms$when %in% c(0, 100) & ms$copyrow < 1 &
      (ms$tipred > 0 | ms$indvarying > 0) & ms$param > 0)

  # The default must not admit a single when == 100 row.
  expect_false(any(ctsem:::.ctMatsetupFreeRows(ms)[ms$when == 100]))
  expect_false(any(ctsem:::.ctMatsetupFreeRows(ms, defining = TRUE)[ms$when == 100]))
})

test_that("when = 0 alone excludes the -1 population rows", {
  # R/ctData.R's laplaceprior block, which is `when == 0` and not c(0,-1).
  ms <- .ms_grid()
  expect_equal(
    ctsem:::.ctMatsetupFreeRows(ms, when = 0, defining = TRUE),
    ms$param > 0 & ms$when == 0 & ms$copyrow < 1)
  expect_false(any(ctsem:::.ctMatsetupFreeRows(ms, when = 0)[ms$when == -1]))
})

test_that("the never-used sentinel is never selected", {
  ms <- .ms_grid()
  for (args in list(list(), list(defining = TRUE), list(varying = TRUE),
                    list(free = FALSE), list(when = c(0, 100)))) {
    keep <- do.call(ctsem:::.ctMatsetupFreeRows, c(list(ms), args))
    expect_false(any(keep[ms$when == -999]))
  }
})

test_that("ctMatsetupFreePars returns one ordered row per parameter", {
  ms <- .ms_grid()
  out <- ctsem:::ctMatsetupFreePars(ms)
  expect_false(any(duplicated(out$param)))
  expect_equal(out$param, sort(out$param))
  expect_true(all(out$param > 0))
  expect_true(all(out$when %in% c(0, -1)))
})

test_that("no inline spelling of the predicate survives in R/", {
  # The point of the helper is that there is one reading of this question. A new
  # inline copy is how the seven disagreeing spellings accumulated, so this is
  # the check that stops it happening again.
  rdir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"),
    mustWork = FALSE)
  skip_if_not(dir.exists(rdir), "package source not available")
  src <- unlist(lapply(list.files(rdir, pattern = "\\.[Rr]$", full.names = TRUE),
    readLines, warn = FALSE))
  src <- src[!grepl("^\\s*#", src)]
  expect_equal(grep("when\\s*%in%\\s*c\\(0\\s*,\\s*-1\\)", src, value = TRUE),
    character(0))
  expect_equal(grep("when\\s*%in%\\s*c\\(0\\s*,\\s*100\\)", src, value = TRUE),
    character(0))
})
