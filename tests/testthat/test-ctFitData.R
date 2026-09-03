# ctFitData(): a backend-neutral accessor for the data a fit was built from.
#
# Motivated by an asymmetry found while assessing what a downstream package
# can reach off a fit object: `fit$data` means two different things depending
# on backend -- the prepared (wide) standata list for a stan fit, the original
# long data.frame for a julia fit. Same field name, unrelated shape. This
# tests that ctFitData() returns the same shape (the long data.frame) either
# way, using the shared internal `.ctFitLongData()`.

test_that("ctFitData() returns the original long data for a stan fit", {
  data(ctstantestfit)
  d <- ctFitData(ctstantestfit)

  expect_s3_class(d, "data.frame")
  expect_true(all(c("id", "time", "Y1", "Y2") %in% names(d)))

  # Compared against the reconstruction ctsem's own internals use, rather
  # than a hard-coded row count that would drift if the shipped example fit
  # ever changes.
  expect_equal(d, ctsem:::.ctFitLongData(ctstantestfit))
})

test_that("ctFitData() rejects non-fit objects clearly", {
  expect_error(ctFitData(list()), "not a ctsem fit")
})

test_that("ctFitData() agrees in shape with a julia fit's own $data", {
  skip_on_cran()
  skip_without_julia()

  set.seed(31)
  dat <- do.call(rbind, lapply(1:6, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * 0.5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))

  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 1, optimcontrol = list(estonly = TRUE))))

  d <- ctFitData(fit)
  expect_s3_class(d, "data.frame")
  expect_equal(sort(names(d)), sort(names(fit$data)))
  expect_equal(nrow(d), nrow(dat))
})
