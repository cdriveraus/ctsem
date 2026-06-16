suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("ctStanModel exposes a pars-backed matrix view", {
  model <- ctModel(type='ct',
    n.latent=2, latentNames=c('eta1','eta2'),
    n.manifest=2, manifestNames=c('Y1','Y2'),
    LAMBDA=diag(2),
    silent=TRUE)
  
  expect_s3_class(model, 'ctStanModel')
  expect_true('matrices' %in% names(model))
  expect_type(model[['matrices']], 'character')
  expect_true(is.list(ctModelMatrices(model)))
  expect_equal(ctModelMatrices(model)$DRIFT, model$matrices$DRIFT)
  expect_equal(model$matr$DRIFT, model$matrices$DRIFT)
  
  model$matrices$DRIFT[1,2] <- 'cross'
  driftrow <- model$pars[
    model$pars$matrix == 'DRIFT' & model$pars$row == 1 & model$pars$col == 2,]
  expect_equal(driftrow$param, 'cross')
  expect_true(is.na(driftrow$value))
  
  model$matrices$DRIFT[2,1] <- -0.2
  driftrow <- model$pars[
    model$pars$matrix == 'DRIFT' & model$pars$row == 2 & model$pars$col == 1,]
  expect_true(is.na(driftrow$param))
  expect_equal(driftrow$value, -0.2)
  expect_true(is.na(driftrow$transform))
  expect_false(driftrow$indvarying)
})

test_that("ctStanModel matrix view can free fixed matrix elements", {
  model <- ctModel(type='ct',
    n.latent=2, latentNames=c('eta1','eta2'),
    n.manifest=2, manifestNames=c('Y1','Y2'),
    LAMBDA=diag(2),
    silent=TRUE)
  
  model$matrices$LAMBDA[1,1] <- 'loading'
  lambdarow <- model$pars[
    model$pars$matrix == 'LAMBDA' & model$pars$row == 1 & model$pars$col == 1,]
  
  expect_equal(lambdarow$param, 'loading')
  expect_true(is.na(lambdarow$value))
  expect_false(is.na(lambdarow$transform))
  expect_equal(lambdarow$sdscale, 1)
})

test_that("ctStanModel detached matrix edits update pars when assigned back", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    LAMBDA=matrix(1),
    silent=TRUE)
  
  mats <- model$matrices
  mats$CINT[1,1] <- 'cint'
  model$matrices <- mats
  
  cintrow <- model$pars[
    model$pars$matrix == 'CINT' & model$pars$row == 1 & model$pars$col == 1,]
  expect_equal(cintrow$param, 'cint')
  expect_true(cintrow$indvarying)
  expect_equal(cintrow$sdscale, 1)
  
  mats <- ctModelMatrices(model)
  mats$CINT[1,1] <- 0
  ctModelMatrices(model) <- mats
  
  cintrow <- model$pars[
    model$pars$matrix == 'CINT' & model$pars$row == 1 & model$pars$col == 1,]
  expect_true(is.na(cintrow$param))
  expect_equal(cintrow$value, 0)
  expect_false(cintrow$indvarying)
})
