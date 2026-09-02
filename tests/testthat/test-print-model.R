# Printing a model.
#
# The thing worth testing is not the layout but the claim: a model whose
# `indvarying` defaults leave four parameters varying is a multilevel model
# estimating a 4x4 population covariance, and printing it should say so. That
# was true before and invisible -- present as a logical column in a dumped data
# frame, absent as a fact.

test_that("it states the model's size and what will be estimated", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), Tpoints = 5))
  out <- capture.output(print(m))
  text <- paste(out, collapse = "\n")
  expect_match(text, "continuous time")
  expect_match(text, "2 latent \\(eta1, eta2\\)")
  expect_match(text, "free parameters")
  # The implication, not just the flag.
  expect_match(text, "multilevel")
  expect_match(text, "4 x 4 population covariance")
  expect_match(text, "6 correlations")
  # Matrices are named on both axes, so a cell does not have to be located by
  # counting rows.
  expect_match(text, "DRIFT")
  expect_true(any(grepl("^eta1 ", out)))
})

test_that("a fixed effects model says so rather than staying silent", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0),
    Tpoints = 5))
  text <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(text, "fixed effects model")
  expect_false(grepl("multilevel", text))
})

test_that("predictors and non-Gaussian indicators are named", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = "eta1",
    LAMBDA = matrix(1, 2, 1), n.TIpred = 1, TIpredNames = "age", Tpoints = 5))
  m$manifesttype <- c(1L, 0L)
  text <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(text, "time independent predictor: age")
  # Binary indicators change the measurement model entirely, so a model
  # carrying one should not look like any other model when printed.
  expect_match(text, "binary")
  expect_match(text, "Y1")
})

test_that("printing returns the model invisibly and does not alter it", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    Tpoints = 5))
  invisible(capture.output(result <- print(m)))
  expect_identical(result$pars, m$pars)
})

test_that("matrices can be suppressed for a model too large to read", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    Tpoints = 5))
  text <- paste(capture.output(print(m, matrices = FALSE)), collapse = "\n")
  expect_false(grepl("DRIFT", text))
  expect_match(text, "ctModelMatrices")
})

# A fit gained a second class ('ctFit') after 3.11.1, which turned every
# `class(x) %in% ...` condition into a length-two logical and so into an error
# under R >= 4.2. The bundled ctstantestfit was saved before that change and
# carries one class, so no test using it can catch this; the classes are set
# explicitly here for that reason. `R CMD check --run-donttest` caught it as an
# error inside plot(fit).
test_that("post-fit functions accept a fit carrying both of its classes", {
  skip_on_cran()
  data("ctstantestfit", package = "ctsem")
  fit <- ctstantestfit
  class(fit) <- c("ctStanFit", "ctFit")
  expect_length(class(fit), 2L)

  # Each of these read the class with a membership test and errored with
  # "the condition has length > 1".
  expect_error(suppressWarnings(suppressMessages(ctPlotPosterior(fit))),
    regexp = NA)
  expect_error(suppressWarnings(suppressMessages(
    ctModelLatex(fit, folder = tempdir(), open = FALSE, compile = FALSE))),
    regexp = NA)

  # And the guard still refuses something that is neither.
  expect_error(ctPlotPosterior(list(a = 1)), regexp = "ctStanFit")
})
