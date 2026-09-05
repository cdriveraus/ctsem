# A control-list name that only one backend can honour must be refused by name
# on the other, not accepted and ignored. Every assertion below is a name that
# used to be silent: the julia path read a handful of optimcontrol names and
# never looked at the rest, backendcontrol reached nothing on the stan path, and
# `cores` on the two accessors was consumed before the julia branch returned.
#
# None of these needs Julia or a fit -- the checks run before any data
# preparation -- so this file is deliberately cheap and always runs.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

.cs_model <- function() ctModel(type = 'ct',
  n.latent = 1, latentNames = 'eta1',
  n.manifest = 1, manifestNames = 'Y1',
  LAMBDA = matrix(1), silent = TRUE)

.cs_data <- function() data.frame(id = 1, time = 1:2, Y1 = c(0, NA))

.cs_fit <- function(...) ctFit(datalong = .cs_data(), model = .cs_model(),
  fit = FALSE, ...)

test_that("stan-only optimcontrol names are refused by name on backend='julia'", {
  # The message must name the argument and say what to use instead; a refusal
  # that only says no leaves the caller where the silence did.
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(tol = 1e-4)),
    "optimcontrol\\$tol is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(tol = 1e-4)),
    "backendcontrol\\$g_tol")

  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(initsd = .1)),
    "optimcontrol\\$initsd is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(initsd = .1)),
    "inits")

  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stallretries = 5)),
    "optimcontrol\\$carefulfit")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stalltol = 1)),
    "optimcontrol\\$carefulfit")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(subsamplesize = .5)),
    "optimcontrol\\$subsamplesize is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(parsteps = 1L)),
    "optimcontrol\\$parsteps is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(nsubsets = 4)),
    "optimcontrol\\$nsubsets is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(lproughnesstarget = .3)),
    "optimcontrol\\$lproughnesstarget is only available for backend='stan'")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stochasticTolAdjust = 10)),
    "optimcontrol\\$stochasticTolAdjust is only available for backend='stan'")
})

test_that("optimcontrol$stochastic is refused on julia only when it asks for something", {
  # FALSE is what the julia optimiser does anyway, and several existing calls
  # pass it, so refusing it would break honest code for no gain. TRUE asks for a
  # stochastic phase that does not exist.
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stochastic = TRUE)),
    "optimcontrol\\$stochastic=TRUE is only available for backend='stan'")
  expect_type(.cs_fit(backend = 'julia', optimcontrol = list(stochastic = FALSE)),
    'list')
})

test_that("julia-only optimcontrol names are refused by name on backend='stan'", {
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(gradient = 'adjoint')),
    "optimcontrol\\$gradient is only available for backend='julia'")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(datastart = FALSE)),
    "optimcontrol\\$datastart is only available for backend='julia'")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(datastart = FALSE)),
    "inits")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(callback = function(...) NULL)),
    "optimcontrol\\$callback is only available for backend='julia'")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(saveEffects = TRUE)),
    "optimcontrol\\$saveEffects is only available for backend='julia'")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(saveEffects = TRUE)),
    "ctSubjectPars")
})

test_that("names both backends honour are accepted on both", {
  for(be in c('stan','julia')){
    expect_type(.cs_fit(backend = be, optimcontrol = list(estonly = TRUE,
      carefulfit = FALSE, finishsamples = 20, uncertainty = 'hessian',
      uncertaintyDraws = 'auto', uncertaintyControl = list())), 'list')
  }
})

test_that("an unrecognised optimcontrol name is refused on both backends", {
  # stanoptimis() has no ... , so a misspelling was already an error on stan --
  # but only with optimize=TRUE, since optimcontrol never reaches it when
  # sampling. On julia nothing looked at all.
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(tolerance = 1e-4)),
    "Unrecognised optimcontrol name")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(tolerance = 1e-4)),
    "Unrecognised optimcontrol name")
})

test_that("backendcontrol is refused by name on backend='stan'", {
  expect_error(.cs_fit(backend = 'stan', backendcontrol = list(g_tol = 1e-6)),
    "backendcontrol is only used with backend='julia'")
  expect_error(.cs_fit(backend = 'stan', backendcontrol = list(g_tol = 1e-6)),
    "optimcontrol")
  # The default empty list must stay silent -- every stan call passes it.
  expect_type(.cs_fit(backend = 'stan', backendcontrol = list()), 'list')
})

test_that("cores is refused by name on ctKalmanArray and ctSubjectPars for a julia fit", {
  # A stub carrying only the class: both refusals fire before anything reads
  # the fit, which is the point -- the argument is wrong whatever the fit holds.
  stub <- structure(list(), class = 'ctJuliaFit')
  expect_error(ctKalmanArray(stub, cores = 4),
    "ctKalmanArray\\(cores=\\) is only available for backend='stan'")
  expect_error(ctKalmanArray(stub, cores = 4), "ctFit\\(cores=\\)")
  expect_error(ctSubjectPars(stub, cores = 4),
    "ctSubjectPars\\(cores=\\) is only available for backend='stan'")
  expect_error(ctSubjectPars(stub, cores = 4), "ctFit\\(cores=\\)")
  # Not supplying it must get past the check: the stub then fails further in,
  # with a message that is not about cores.
  for(call in list(function() ctKalmanArray(stub), function() ctSubjectPars(stub))){
    msg <- tryCatch({call(); ''}, error = function(e) conditionMessage(e))
    expect_false(grepl('cores=', msg, fixed = TRUE))
  }
})
