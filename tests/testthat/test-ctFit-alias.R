suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("ctStanFit remains an alias for ctFit", {
  expect_identical(ctStanFit, ctFit)
  expect_identical(formals(ctStanFit), formals(ctFit))
})

test_that("ctStan-prefixed helper names remain aliases", {
  expect_identical(ctStanModel, ctModelConvertOMX)
  expect_identical(ctStanContinuousPars, ctContinuousPars)
  expect_identical(ctStanDiscretePars, ctDiscretePars)
  expect_identical(ctStanDiscreteParsPlot, ctDiscreteParsPlot)
  expect_identical(ctStanFitUpdate, ctFitUpdate)
  expect_identical(ctStanGenerate, ctGenerateFromPriors)
  expect_identical(ctStanGenerateFromFit, ctGenerateFromFit)
  expect_identical(ctStanKalman, ctKalmanArray)
  expect_identical(ctStanParnames, ctRawParnames)
  expect_identical(ctStanPlotPost, ctPlotPosterior)
  expect_identical(ctStanPostPredict, ctPostPredict)
  expect_identical(ctStanSubjectPars, ctSubjectPars)
  expect_identical(ctStanTIpredeffects, ctTIpredEffects)
})
