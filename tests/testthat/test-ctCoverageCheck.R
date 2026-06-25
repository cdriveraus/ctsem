
  library(ctsem)
  library(testthat)


test_that("ctModelCoverage_check generation cores default to fit_cores", {
  expect_identical(
    formals(ctModelCoverage_check)$generate_cores,
    quote(fit_cores)
  )
})
