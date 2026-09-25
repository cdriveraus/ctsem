# `optimcontrol` is one vocabulary: a name in it means the same thing on both
# backends, and the shared ones are honoured by both. What is left is a
# capability one backend has and the other does not, and there the *value*
# decides -- a value asking for a missing capability is refused by name, a value
# describing what the backend already does is accepted.
#
# All but one of these run without Julia or a fit -- the checks run before any
# data preparation -- so this file is deliberately cheap and always runs. The
# exception fits a ten-subject model, because what it checks is that a value
# reaches the loop it caps, and is guarded by skip_without_julia().

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

.cs_model <- function() ctModel(type = 'ct',
  n.latent = 1, latentNames = 'eta1',
  n.manifest = 1, manifestNames = 'Y1',
  LAMBDA = matrix(1), silent = TRUE)

.cs_data <- function() data.frame(id = 1, time = 1:2, Y1 = c(0, NA))

.cs_fit <- function(...) ctFit(datalong = .cs_data(), model = .cs_model(),
  fit = FALSE, ...)

test_that("the shared optimiser vocabulary is accepted on both backends", {
  # These used to be split: `tol` and `initsd` were refused on julia, and
  # `g_tol`, `x_tol`, `maxiter` and `lbfgs_memory` existed only as
  # `backendcontrol` names that reached nothing on stan.
  shared <- list(tol = 1e-6, g_tol = 1e-7, x_tol = 1e-9, maxiter = 500,
    lbfgs_memory = 20, initsd = .05, estonly = TRUE, carefulfit = FALSE,
    finishsamples = 20, uncertainty = 'hessian', uncertaintyDraws = 'auto',
    uncertaintyControl = list())
  for(be in c('stan','julia')){
    expect_type(.cs_fit(backend = be, optimcontrol = shared), 'list')
  }
})

test_that("stanoptimis() takes the shared stopping rules, with defaults that change nothing", {
  f <- formals(stanoptimis)
  for(nm in c('g_tol','x_tol','maxiter','lbfgs_memory')) expect_true(nm %in% names(f))
  # NULL means "leave the optimizer as it was", so an existing stan call is
  # unaffected by their existence.
  for(nm in c('g_tol','x_tol','maxiter','lbfgs_memory')) expect_null(f[[nm]])
  expect_equal(f$tol, 1e-8)
})

test_that("a stan-only capability is refused on julia when the value asks for it", {
  # The julia optimiser is L-BFGS over the full data: no stochastic gradient
  # phase, and no gradient bar for calling a fit stalled -- it judges a stall
  # by its progress.
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stochastic = TRUE)),
    "optimcontrol\\$stochastic asks for stochastic gradient descent")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(nsubsets = 4)),
    "optimcontrol\\$nsubsets")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(subsamplesize = .5)),
    "optimcontrol\\$carefulfit")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(parsteps = 1L)),
    "optimcontrol\\$parsteps")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stalltol = 1)),
    "optimcontrol\\$stallwindow")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(lproughnesstarget = .3)),
    "optimcontrol\\$lproughnesstarget")
  expect_error(.cs_fit(backend = 'julia', optimcontrol = list(stochasticTolAdjust = 10)),
    "optimcontrol\\$stochasticTolAdjust")
})

test_that("a value describing what the backend already does is accepted", {
  # stochastic=FALSE is the deterministic optimiser julia already runs;
  # nsubsets=1 and subsamplesize=1 are no split and no subset; parsteps=c() asks
  # for no stepwise pass.
  expect_type(.cs_fit(backend = 'julia', optimcontrol = list(stochastic = FALSE,
    nsubsets = 1, subsamplesize = 1, parsteps = c())), 'list')
  # Stan's autodiff is reverse mode, so 'adjoint' is a true description of it;
  # 'forward' names a second gradient only the julia engine has.
  expect_type(.cs_fit(backend = 'stan', optimcontrol = list(gradient = 'adjoint',
    datastart = FALSE, callback = NULL, saveEffects = FALSE, progress = FALSE,
    tipredMissingIncludeOutcome = TRUE)), 'list')
})

test_that("a julia-only capability is refused on stan when the value asks for it", {
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(gradient = 'forward')),
    "optimcontrol\\$gradient selects the julia engine's gradient")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(datastart = TRUE)),
    "optimcontrol\\$initsd")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(callback = function(...) NULL)),
    "optimcontrol\\$callback")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(saveEffects = TRUE)),
    "ctSubjectPars")
  expect_error(.cs_fit(backend = 'stan', optimcontrol = list(progress = TRUE)),
    "optimcontrol\\$progress")
  expect_error(.cs_fit(backend = 'stan',
    optimcontrol = list(tipredMissingIncludeOutcome = FALSE)),
    "optimcontrol\\$tipredMissingIncludeOutcome")
})

test_that("stallretries is one name, and both backends take it", {
  # It was in the capability split twice, stan-only and julia-only, and a list
  # lookup returns the first entry: julia refused every value but 0 with stan's
  # message, and stan accepted it and then dropped it with the julia-only names
  # on the way to stanoptimis(), which restarted twice whatever was asked.
  split <- ctsem:::.ctOptimcontrolSplit()
  expect_equal(anyDuplicated(names(split)), 0L)
  expect_length(intersect(ctsem:::.ctOptimcontrolShared, names(split)), 0L)
  for(be in c('stan','julia')){
    expect_type(.cs_fit(backend = be, optimcontrol = list(stallretries = 1L)), 'list')
  }
  # What ctFit() hands stanoptimis() on the stan path, bound there by name.
  tostan <- ctsem:::.ctOptimcontrolForStan(list(stallretries = 1L, tol = 1e-6))
  expect_equal(tostan$stallretries, 1L)
  expect_true('stallretries' %in% names(formals(stanoptimis)))
})

test_that("stallretries is the cap on the julia escape loop", {
  skip_without_julia()
  # An escape that always has somewhere to go -- back to where the stage
  # stopped, which a resumed stage cannot do worse than -- so the loop runs
  # until the cap stops it, and the escapes offered are the cap it read.
  offered <- 0L
  testthat::local_mocked_bindings(.ctBackendStallEscape = function(result, ...) {
    offered <<- offered + 1L
    as.numeric(result$minimizer)
  }, .package = 'ctsem')
  set.seed(1)
  dat <- do.call(rbind, lapply(1:10, function(i) {
    eta <- as.numeric(stats::filter(stats::rnorm(6), 0.6, method = 'recursive'))
    data.frame(id = i, time = 0:5, Y1 = eta + stats::rnorm(6, 0, .5))
  }))
  escapes <- function(n) {
    offered <<- 0L
    fit <- suppressWarnings(suppressMessages(ctFit(dat, .cs_model(),
      backend = 'julia', cores = 1, verbose = 0,
      optimcontrol = list(stallretries = n, estonly = TRUE, carefulfit = FALSE))))
    c(offered = offered, recorded = fit$optim$stall_escapes)
  }
  expect_equal(escapes(1L), c(offered = 1L, recorded = 1L))
  # Above the default of 2, so it is the value that stops the loop.
  expect_equal(escapes(3L), c(offered = 3L, recorded = 3L))
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

test_that("backendcontrol is gone, and says where its settings went", {
  # It reaches `...`, so without this it would be accepted in silence -- which
  # is the fault the merge was for.
  expect_error(.cs_fit(backend = 'julia', backendcontrol = list(g_tol = 1e-6)),
    "backendcontrol has been merged into optimcontrol")
  expect_error(.cs_fit(backend = 'stan', backendcontrol = list(g_tol = 1e-6)),
    "ctJuliaSetup")
  expect_false('backendcontrol' %in% names(formals(ctFit)))
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
