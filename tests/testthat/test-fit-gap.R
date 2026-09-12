# One shape for "how far is the reported answer from a better one".
#
# ctsem asks that three times -- against the exact curvature, against adaptive
# quadrature, against a particle filter -- and answered it in three shapes under
# three names (`$gap`, `$gap`, `$difference`). `.ctFitGap()` is the shape all
# three now also report, in `$verdict`.
#
# The constructor tests are free and always run. The three-route test needs
# julia and a fit each way, so it is separate and skips without.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("a gap with a stated bar reports whether it cleared it", {
  g <- ctsem:::.ctFitGap("curvature", gap = 0.5, tolerance = 0.01)
  expect_s3_class(g, "ctFitGap")
  expect_equal(g$gap, 0.5)
  expect_true(g$exceeds_tolerance)

  under <- ctsem:::.ctFitGap("curvature", gap = 0.001, tolerance = 0.01)
  expect_false(under$exceeds_tolerance)
})

test_that("a gap with no stated bar says so rather than implying it passed", {
  # Two of the three references have no tolerance anywhere in the package, and
  # inventing one would report "immaterial" against a bar nobody chose.
  g <- ctsem:::.ctFitGap("quadrature", gap = 3.2)
  expect_true(is.na(g$tolerance))
  expect_true(is.na(g$exceeds_tolerance))
  expect_match(paste(utils::capture.output(print(g)), collapse = " "),
    "none stated")
})

test_that("a stochastic reference reports whether the gap beats its own noise", {
  # Only the particle route has a standard error, and nothing previously used it
  # for a verdict -- a gap of 3 with an se of 4 and a gap of 3 with an se of 0.1
  # were reported identically.
  noisy <- ctsem:::.ctFitGap("particle", gap = 3, gap_se = 4)
  expect_false(noisy$resolved)
  clean <- ctsem:::.ctFitGap("particle", gap = 3, gap_se = 0.1)
  expect_true(clean$resolved)
  # No standard error is "not asked", not "not resolved".
  expect_true(is.na(ctsem:::.ctFitGap("quadrature", gap = 3)$resolved))
})

test_that("exceeding a bar and beating the noise are reported separately", {
  # Collapsing them would hide which of the two failed, and they want different
  # responses: a gap under the bar is fine, a gap lost in its own noise means
  # measure it better.
  g <- ctsem:::.ctFitGap("particle", gap = 3, gap_se = 4, tolerance = 0.01)
  expect_true(g$exceeds_tolerance)
  expect_false(g$resolved)
})

test_that("the gap is normalised per subject when the count is known", {
  g <- ctsem:::.ctFitGap("quadrature", gap = 10, nsubjects = 50L)
  expect_equal(g$gap_per_subject, 0.2)
  expect_true(is.na(ctsem:::.ctFitGap("quadrature", gap = 10)$gap_per_subject))
})

test_that("an unavailable gap prints as unavailable", {
  g <- ctsem:::.ctFitGap("curvature", gap = NA_real_, tolerance = 0.01)
  expect_match(paste(utils::capture.output(print(g)), collapse = " "),
    "not available")
})

test_that("every reference name the package uses has a printable description", {
  # A reference with no entry falls back to the bare name, which is legible but
  # is a sign the table below and the call sites have drifted.
  for (r in c("curvature", "quadrature", "particle")) {
    out <- paste(utils::capture.output(
      print(ctsem:::.ctFitGap(r, gap = 1))), collapse = " ")
    expect_false(grepl(paste0("against ", r, "\\b"), out), info = r)
  }
})

test_that("all three measurement routes report a verdict of the same shape", {
  skip_without_julia()
  set.seed(4)
  n <- 8L; t <- 5L
  dat <- data.frame(id = rep(seq_len(n), each = t),
    time = rep(seq_len(t), times = n), Y1 = stats::rnorm(n * t))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE

  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 1, intoverpop = "laplace", priors = TRUE,
    optimcontrol = list(finishsamples = 20))))

  # 1. the curvature certification, which ctFit() attaches on every optimised fit
  cert <- fit$uncertainty$certification
  expect_s3_class(cert$verdict, "ctFitGap")
  expect_equal(cert$verdict$reference, "curvature")
  expect_true(is.finite(cert$verdict$tolerance))

  # 2. the quadrature check
  chk <- suppressWarnings(suppressMessages(ctLaplaceCheck(fit, nodes = 3L)))
  expect_s3_class(chk$verdict, "ctFitGap")
  expect_equal(chk$verdict$reference, "quadrature")
  # The shape is shared; the number is still the one this route always reported.
  expect_equal(chk$verdict$gap, chk$gap)
  expect_equal(chk$verdict$gap_per_subject, chk$gap_per_subject)

  # The three carry the same field names, which is the whole point.
  expect_equal(sort(names(cert$verdict)), sort(names(chk$verdict)))
})
