suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("ctStanFit remains an alias for ctFit", {
  expect_identical(ctStanFit, ctFit)
  expect_identical(formals(ctStanFit), formals(ctFit))
})

test_that("ctFit accepts model argument and deprecates ctstanmodel", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    LAMBDA=matrix(1),
    silent=TRUE)
  dat <- data.frame(id=1, time=1:2, Y1=c(0, NA))

  expect_type(ctFit(datalong=dat, model=model, fit=FALSE), 'list')
  expect_warning(
    expect_type(ctFit(datalong=dat, ctstanmodel=model, fit=FALSE), 'list'),
    'ctstanmodel argument is deprecated')
  expect_error(
    ctFit(datalong=dat, model=model, ctstanmodel=model, fit=FALSE),
    'Use only one')
})

test_that("ctFit stan backend rejects the retired optimcontrol$is flag with a clear error", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    LAMBDA=matrix(1),
    silent=TRUE)
  dat <- data.frame(id=1, time=1:2, Y1=c(0, NA))

  expect_error(
    ctFit(datalong=dat, model=model, fit=FALSE, backend='stan',
      optimcontrol=list(is=TRUE)),
    regexp="optimcontrol\\$uncertainty")
})

test_that("ctStan-prefixed helper names remain aliases", {
  expect_identical(ctStanModel, ctModelConvertOMX)
  expect_identical(ctStanContinuousPars, ctSummaryMatrices)
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

test_that("the defunct model-update stubs error and name their replacement", {
  # Both are exported, so both are reachable from user code that has not been
  # updated. What matters is that they stop rather than return something, and
  # that the message says where to go.
  expect_error(ctFitUpdateModel(), 'ctFitUpdate')
  expect_error(ctStanUpdModel(), 'ctFitUpdate')
  expect_error(ctFitUpdateModel(), 'defunct')
  expect_error(ctStanUpdModel(), 'defunct')
})
