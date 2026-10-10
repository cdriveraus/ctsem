# kalmanvec selects from .ctKalmanSeries, which is declared rather than read
# off the output so that a caller can offer it before predicting. So the two
# have to agree: every series the output holds is in the set, and nothing in
# the set is missing from an output that asked for everything.

test_that("the declared series are the ones the prediction output holds", {
  k <- ctPredict(ctstantestfit, subjects = 1, standardisederrors = TRUE)
  series <- setdiff(unique(as.character(k$Element)), "llrow")
  series <- series[!grepl("cov$", series)]
  expect_setequal(series, .ctKalmanSeries)

  expect_s3_class(plot(k, kalmanvec = c("y", "ysmooth"), plot = FALSE), "ggplot")
  expect_error(plot(k, kalmanvec = c("y", "ysmoothed"), plot = FALSE), "kalmanvec takes")
})

# A row whose first variable is missing is standardised by the factor of its
# observed block. The observed block of the full factor -- what this used to
# take -- divided the second variable's error by its sd conditional on the
# missing first one.
test_that("standardised errors on a partly missing row use the observed block", {
  k <- suppressMessages(ctKalmanArray(ctstantestfit, standardisederrors = TRUE))
  rows <- which(is.na(k$y[, 1]) & rowSums(!is.na(k$y[, -1, drop = FALSE])) > 0)
  skip_if(!length(rows), "the test fit has no row with only a later variable observed")
  for (r in rows) {
    obs <- which(!is.na(k$y[r, ]))
    S <- matrix(k$ypriorcov[1, r, obs, obs], length(obs))
    expect_equal(as.numeric(k$errstdprior[1, r, obs]),
      as.numeric(solve(t(chol(S)), k$errprior[1, r, obs])), tolerance = 1e-8)
  }
})
