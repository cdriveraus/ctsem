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
  expect_error(ctFit(dat, .guard_model(), fit = FALSE, backend = 'stan'),
    "Columns not found in the data: Y2 (manifest); TI1 (time independent predictor).", fixed = TRUE)
})

test_that("values that are not numbers are refused, not turned into NA", {
  dat <- as.data.frame(ctstantestdat)
  dat$Y1 <- as.character(dat$Y1)
  dat$Y1[c(3, 8)] <- c("n/a", "missing")
  expect_error(ctFit(dat, .guard_model(), fit = FALSE, backend = 'stan'),
    "Column Y1 holds values that are not numbers: 'n/a', 'missing'.", fixed = TRUE)

  # Numbers held as text, and logical columns, still read.
  dat <- as.data.frame(ctstantestdat)
  dat$Y1 <- as.character(dat$Y1)
  dat$TI1 <- dat$TI1 > 0
  expect_no_error(suppressWarnings(suppressMessages(ctFit(dat, .guard_model(), fit = FALSE, backend = 'stan'))))
})

# sameInitialTimes = TRUE adds the rows to the data ctFit() works on, so both
# backends read them. Added to stan's prepared data alone, they stopped every
# julia fit ("Prepared TD predictor rows do not match the fitted data").
test_that("sameInitialTimes adds an empty first row where a subject has none", {
  dat <- as.data.frame(ctstantestdat)
  dat$TD1[is.na(dat$TD1)] <- 0
  dat <- dat[!(dat$id %in% c(2, 5) & dat$time == 0), ]
  m <- suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2), TDpredNames = "TD1",
    TIpredNames = "TI2"))
  added <- ctsem:::.ctSameInitialTimes(dat, m)
  expect_equal(nrow(added), nrow(dat) + 2)
  first <- added[added$id %in% c(2, 5) & added$time == 0, ]
  expect_true(all(is.na(first[, c("Y1", "Y2")])))
  expect_equal(first$TD1, c(0, 0))
  expect_equal(first$TI2, dat$TI2[match(c(2, 5), dat$id)])
  # A replayed fit's data have the rows already.
  expect_identical(ctsem:::.ctSameInitialTimes(added, m), added)
  # No warning of missing TD predictors: the added rows have no input.
  s <- expect_no_warning(suppressMessages(ctFit(dat, m, fit = FALSE,
    backend = 'stan', sameInitialTimes = TRUE)))
  expect_equal(s$standata$ndatapoints, nrow(added))
})

test_that("a julia fit with sameInitialTimes is the fit of the data with the rows added", {
  skip_without_julia()
  dat <- image_data(nsub = 12L, nobs = 6L, seed = 1L)
  dat <- dat[!(dat$id %in% c(3, 4, 9) & dat$time < 2), ]
  byhand <- rbind(dat, data.frame(id = c(3, 4, 9), time = 0, y1 = NA, y2 = NA))
  byhand <- byhand[order(byhand$id, byhand$time), ]
  set.seed(1)
  fa <- suppressMessages(ctFit(dat, image_model(), cores = 1, verbose = 0,
    sameInitialTimes = TRUE))
  set.seed(1)
  fb <- suppressMessages(ctFit(byhand, image_model(), cores = 1, verbose = 0))
  expect_equal(nrow(ctsem:::.ctFitLongData(fa)), nrow(byhand))
  expect_equal(ctsem:::.ctFitOptimisedLoglik(fa), ctsem:::.ctFitOptimisedLoglik(fb),
    tolerance = 1e-8)
  expect_equal(fa$estimate$raw, fb$estimate$raw, tolerance = 1e-6)
})
