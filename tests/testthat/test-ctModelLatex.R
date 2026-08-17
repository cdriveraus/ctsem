suppressWarnings(suppressPackageStartupMessages(library(ctsem)))
library(testthat)

test_that("ctModelLatex subject distribution excludes calculated pars", {
  m <- suppressWarnings(suppressMessages(ctModel(
    type='ct',
    manifestNames=c('response', 'log_rt'),
    latentNames=c('ddm_drift', 'ddm_boundary_raw', 'ddm_ndt_raw'),
    TDpredNames=c('condition_num', 'stimulus_direction'),
    TDPREDEFFECT=matrix(0, nrow=3, ncol=2),
    LAMBDA=matrix(0, nrow=2, ncol=3),
    manifesttype=c(1, 0),
    MANIFESTMEANS=matrix(c(0, 0), nrow=2),
    MANIFESTVAR=matrix(c(0, 0, 0, 'rt_residual_sd'), nrow=2, byrow=TRUE),
    DRIFT=matrix(0, nrow=3, ncol=3),
    CINT=matrix(0, nrow=3, ncol=1),
    DIFFUSION=matrix(0, nrow=3, ncol=3),
    T0MEANS=matrix(c(
      'drift_t0||TRUE',
      'boundary_raw_t0||TRUE',
      'ndt_raw_t0 | 0.1 + param | TRUE'
    ), nrow=3),
    T0VAR=diag(1e-5, 3),
    PARS=c('gamma_v'),
    silent=TRUE)))
  
  m$pars$param[m$pars$matrix == 'PARS'] <- 'exp(gamma_v)'
  m$pars$indvarying[m$pars$matrix == 'PARS'] <- TRUE
  
  latex <- ctModelLatex(m, tex=FALSE, compile=FALSE, open=FALSE)
  
  expect_match(latex, 'drift\\\\_t0')
  expect_match(latex, 'boundary\\\\_raw\\\\_t0')
  expect_match(latex, 'ndt\\\\_raw\\\\_t0')
  expect_no_match(latex, 'exp\\\\(gamma_v\\\\).*_i')
})

test_that("ctModelLatex handles scalar numeric T0VAR display covariance", {
  m <- suppressMessages(ctModel(type='ct', manifestNames='Y1', LAMBDA=diag(1)))
  expect_error(ctModelLatex(m, tex=FALSE, compile=FALSE, open=FALSE), NA)
})
