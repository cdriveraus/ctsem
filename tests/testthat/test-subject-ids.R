# Functions that take subjects read them as the ids in the data, and fall back
# to ctsem's numbering 1..N for integers that are not ids, with a message.
# ctstantestfit's ids are 1..N, which cannot tell the two apart, so its id map
# is relabelled: the filter works by position, and reads ids only through it.

.relabelled_fit <- function() {
  fit <- ctstantestfit
  fit$standata$idmap[, 1] <- 100 + 10 * fit$standata$idmap[, 2]
  fit
}

test_that("subjects resolve as data ids first, then as ctsem's numbering", {
  fit <- .relabelled_fit()
  ids <- .ctResolveSubjects(fit, c(130, 110))
  expect_equal(as.vector(ids), c(3L, 1L))
  expect_true(attr(ids, "realid"))
  expect_equal(attr(ids, "ids"), c(130, 110))

  expect_message(numbered <- .ctResolveSubjects(fit, 1:2), "not ids in the data")
  expect_equal(as.vector(numbered), 1:2)
  expect_false(attr(numbered, "realid"))
  expect_no_message(.ctResolveSubjects(fit, 1:2, realid = FALSE))

  expect_error(.ctResolveSubjects(fit, 999), "Subjects not found: 999")
  expect_error(.ctResolveSubjects(fit, "abc"), "Subjects not found: abc")
  expect_error(.ctResolveSubjects(fit, 110, realid = FALSE), "Subjects not found: 110")
})

test_that("ctPredict and ctExtract give a subject the same by id or by number", {
  fit <- .relabelled_fit()
  byid <- ctPredict(fit, subjects = c(130, 120))
  bynumber <- suppressMessages(ctPredict(fit, subjects = 2:3))
  expect_setequal(as.character(unique(byid$Subject)), c("120", "130"))
  expect_setequal(as.character(unique(bynumber$Subject)), c("2", "3"))
  expect_equal(byid$value, bynumber$value)

  e_id <- ctExtract(fit, subjectMatrices = TRUE, subjects = 120)
  e_number <- ctExtract(fit, subjectMatrices = TRUE, subjects = 2, realid = FALSE)
  expect_equal(e_id$subj_DRIFT, e_number$subj_DRIFT)
})

test_that("ctDiscretePars labels each subject's result with its own id", {
  fit <- .relabelled_fit()
  set.seed(1)
  byid <- ctDiscretePars(fit, subjects = c(130, 120), times = 1, nsamples = 3, cores = 1)
  set.seed(1)
  bynumber <- ctDiscretePars(fit, subjects = 2:3, times = 1, nsamples = 3, cores = 1, realid = FALSE)
  expect_identical(dimnames(byid)$Subject, c("120", "130"))
  expect_identical(dimnames(bynumber)$Subject, c("2", "3"))
  expect_equal(unname(byid), unname(bynumber))
})
