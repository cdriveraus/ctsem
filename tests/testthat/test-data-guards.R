# ctFit checks the columns a model names before reading any of them: the
# missing ones are listed together with their roles, and text that does not
# read as a number is refused rather than turned into NA, which dropped those
# observations silently.

.guard_model <- function() {
  suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2), TIpredNames = "TI1"))
}

test_that("missing columns are listed together, with their roles", {
  dat <- ctstantestdat[, setdiff(colnames(ctstantestdat), c("Y2", "TI1"))]
  expect_error(ctFit(dat, .guard_model(), fit = FALSE),
    "Columns not found in the data: Y2 (manifest); TI1 (time independent predictor).", fixed = TRUE)
})

test_that("values that are not numbers are refused, not turned into NA", {
  dat <- as.data.frame(ctstantestdat)
  dat$Y1 <- as.character(dat$Y1)
  dat$Y1[c(3, 8)] <- c("n/a", "missing")
  expect_error(ctFit(dat, .guard_model(), fit = FALSE),
    "Column Y1 holds values that are not numbers: 'n/a', 'missing'.", fixed = TRUE)

  # Numbers held as text, and logical columns, still read.
  dat <- as.data.frame(ctstantestdat)
  dat$Y1 <- as.character(dat$Y1)
  dat$TI1 <- dat$TI1 > 0
  expect_no_error(suppressWarnings(suppressMessages(ctFit(dat, .guard_model(), fit = FALSE))))
})
