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
