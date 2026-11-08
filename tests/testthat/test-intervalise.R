library(ctsem)
library(testthat)

context("intervalise")

test_that("intervals", {
  data('longexample')
  
  #Then convert to wide format
  wideexample <- ctLongToWide(datalong = longexample, id = "id",
                              time = "time", manifestNames = c("Y1", "Y2", "Y3"),
                              TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
  
  #Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
  wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3,
                        n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"),
                        TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )
  
  dt <- matrix(t(diff(t(wideexample[,paste0('T',0:2)]))), nrow=2,ncol=2,
            dimnames = list(c("1", "2"), c("dT1", "dT2")))
  dt[1,2] = .001
  expect_equal(wide[,c('dT1','dT2')], dt)
})

test_that("columns named as ctsem names them are read by name, whatever their order", {
  # Three people, three occasions, two manifests, two time dependent
  # predictors and one time independent one, every value distinct.
  Tp <- 3
  ordered <- c(ctWideNames(n.manifest = 2, Tpoints = Tp, n.TDpred = 2)[1:(Tp * 4)],
    paste0("T", 0:(Tp - 1)), "TI1")
  x <- matrix(NA_real_, 3, length(ordered), dimnames = list(1:3, ordered))
  for (p in 1:3) {
    x[p, 1:(Tp * 4)] <- 100 * p + seq_len(Tp * 4)
    x[p, paste0("T", 0:(Tp - 1))] <- c(0, 1, 2.5) + p
    x[p, "TI1"] <- -p
  }
  reference <- suppressMessages(ctIntervalise(x, Tpoints = Tp, n.manifest = 2,
    n.TDpred = 2, n.TIpred = 1))
  # Each predictor's occasions together, as a listing might put them.
  shuffled <- x[, c(ordered[1:6], "TD1_T0", "TD1_T1", "TD1_T2", "TD2_T0",
    "TD2_T1", "TD2_T2", "TI1", paste0("T", 0:(Tp - 1)))]
  expect_message(out <- ctIntervalise(shuffled, Tpoints = Tp, n.manifest = 2,
    n.TDpred = 2, n.TIpred = 1), "wide order by name")
  expect_identical(out, reference)
  # A data frame too, with its names checked as given.
  frame <- suppressMessages(ctIntervalise(as.data.frame(shuffled),
    Tpoints = Tp, n.manifest = 2, n.TDpred = 2, n.TIpred = 1))
  expect_equal(unname(as.matrix(frame)), unname(reference))
  expect_identical(colnames(frame), colnames(reference))
})

test_that("an occasion with no time loses that occasion's values and no other's", {
  Tp <- 10
  nm <- c(paste0(rep(c("Y1", "Y2"), Tp), "_T", rep(0:(Tp - 1), each = 2)),
    paste0("TD1_T", 0:(Tp - 1)), paste0("T", 0:(Tp - 1)))
  x <- matrix(c(1:(2 * Tp), 100 + 0:(Tp - 1), 0:(Tp - 1)), 1, length(nm),
    dimnames = list("1", nm))
  x[1, "T2"] <- NA
  x[1, "T9"] <- NA
  out <- suppressMessages(ctIntervalise(x, Tpoints = Tp, n.manifest = 2, n.TDpred = 1))
  gone <- c("Y1_T2", "Y2_T2", "TD1_T2", "Y1_T9", "Y2_T9", "TD1_T9")
  expect_true(all(is.na(out[1, gone])))
  kept <- setdiff(nm[1:(3 * Tp)], gone)
  expect_equal(out[1, kept], x[1, kept])
})
