# ctReport(): the folder of output, on both backends.
#
# What is worth asserting here is not that each component runs -- each is an
# existing function with its own tests -- but the three things ctReport() itself
# is responsible for:
#
#   1. the files named in the index are the files on disk, and they are not
#      empty. A component that "succeeded" and wrote nothing is the failure
#      mode a status list alone would hide.
#   2. a component a backend does not support is recorded as skipped or failed
#      and the rest of the run still completes. That is the whole design, and
#      it is only actually exercised by handing it a fit that cannot do
#      everything.
#   3. the profile's reference curve. A quadratic log probability drops by
#      k^2/2 at k conditional standard errors, so a well-determined parameter
#      must come back at 2.0 -- which is a real check on .ctReportRawSE()
#      using the conditional (Hessian diagonal) width rather than the marginal
#      one. Scaling by the marginal SE instead is silent, plausible, and wrong,
#      and it made every parameter in a healthy fit look non-quadratic.

.ctReportFolder <- function() {
  d <- file.path(tempdir(), paste0("ctReportTest", as.integer(stats::runif(1, 0, 1e7))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Every file the index claims exists, is non-trivial, and nothing numbered is
# on disk without being in the record.
.ctReportCheckFolder <- function(res) {
  index <- file.path(res$folder, "00-index.md")
  testthat::expect_true(file.exists(index))
  testthat::expect_gt(file.size(index), 200)

  claimed <- unlist(lapply(res$components,
    function(x) if (identical(x$status, "ok")) x$files else NULL))
  testthat::expect_true(length(claimed) > 0)
  for (f in claimed) {
    p <- file.path(res$folder, f)
    testthat::expect_true(file.exists(p), info = paste("missing", f))
    testthat::expect_gt(file.size(p), 100)
  }
  ondisk <- setdiff(list.files(res$folder, pattern = "^[0-9][0-9]-"), "00-index.md")
  testthat::expect_setequal(ondisk, unique(claimed))
  claimed
}

# The profile's grid must be in conditional standard errors, not marginal ones.
# Under the marginal scaling each drop is multiplied by that parameter's ridge
# squared, so the signature of the error is drops that track the ridge and run
# to tens. Asserting the exact 2.0 instead would be asserting that the
# likelihood is quadratic, which is the thing the report exists to question.
.ctReportCheckCalibration <- function(pr) {
  d <- rbind(
    data.frame(drop = pr$table$drop_lo, ridge = pr$table$ridge),
    data.frame(drop = pr$table$drop_hi, ridge = pr$table$ridge))
  d <- d[is.finite(d$drop), ]
  testthat::expect_gt(nrow(d), 4)
  testthat::expect_lt(abs(stats::median(d$drop) - pr$expected), pr$expected / 2)
  testthat::expect_lt(max(d$drop), 4 * pr$expected)
  ridged <- d[is.finite(d$ridge) & d$ridge > 2, ]
  if (nrow(ridged)) {
    # These are the ones that would blow up: ridge 4 means a factor of 16.
    testthat::expect_lt(max(ridged$drop), 3 * pr$expected)
  }
  invisible(d)
}

test_that("ctReport writes the quick components and an index that matches the folder (stan)", {
  skip_on_cran()

  d <- .ctReportFolder()
  res <- suppressWarnings(suppressMessages(
    ctReport(ctstantestfit, folder = d, components = "quick", nsamples = 10,
      profilepoints = 5, verbose = FALSE)))

  expect_identical(normalizePath(res$folder, winslash = "/"),
    normalizePath(d, winslash = "/"))
  files <- .ctReportCheckFolder(res)

  # Every quick component, by name, so a component silently dropped from the
  # set is a failure rather than a shorter list.
  expect_setequal(names(res$components),
    c("summary", "parmatrices", "identification", "profile", "discretepars",
      "network", "predictions", "residuals"))
  expect_true(all(vapply(res$components, function(x) x$status, character(1)) == "ok"))

  # The numbering is the reading order: sorting the names must not reorder them.
  expect_identical(sort(files), files[order(files)])
  expect_true(any(grepl("^01-summary", files)))
  expect_true(any(grepl("^04-loglik-profile[.]pdf$", files)))

  # The index says how to read each file, not only what it is.
  index <- readLines(file.path(d, "00-index.md"))
  expect_true(any(grepl("Backend: stan", index, fixed = TRUE)))
  expect_true(any(grepl("Read this first", index, fixed = TRUE)))
})

test_that("a component the backend cannot do is recorded, and the run continues (stan)", {
  skip_on_cran()

  # No covariates, so tipredeffects has nothing to plot and must skip -- while
  # summary, before it, and covcheck, after it, both still run.
  m <- suppressWarnings(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1, 1, 1),
    MANIFESTVAR = matrix(0.2, 1, 1)))
  set.seed(2)
  dat <- data.frame(id = rep(1:8, each = 5), time = rep(0:4, 8),
    Y1 = stats::rnorm(40))
  fit <- suppressWarnings(suppressMessages(ctFit(dat, m, backend = "stan",
    cores = 1, priors = FALSE, optimcontrol = list(carefulfit = FALSE))))

  d <- .ctReportFolder()
  res <- suppressWarnings(suppressMessages(
    ctReport(fit, folder = d, components = c("summary", "tipredeffects", "covcheck"),
      nsamples = 10, verbose = FALSE)))

  expect_identical(res$components$tipredeffects$status, "skipped")
  expect_match(res$components$tipredeffects$note, "covariate")
  expect_identical(res$components$summary$status, "ok")
  expect_identical(res$components$covcheck$status, "ok")
  .ctReportCheckFolder(res)

  # And the index tells the user, rather than leaving a gap in the numbering.
  index <- readLines(file.path(d, "00-index.md"))
  expect_true(any(grepl("Not produced", index, fixed = TRUE)))
  expect_true(any(grepl("tipredeffects", index, fixed = TRUE)))
})

test_that("the profile's quadratic reference is calibrated: a determined parameter drops by 2 at 2 SE", {
  skip_on_cran()

  pr <- ctsem:::.ctReportProfile(ctstantestfit, npoints = 5L)
  expect_s3_class(pr, "ctReportProfile")
  expect_true(all(c("parameter", "estimate", "se_report", "se_slice", "ridge",
    "drop_lo", "drop_hi", "peak", "flag") %in% names(pr$table)))
  expect_identical(nrow(pr$table), length(ctsem:::.ctFitRawEstimate(ctstantestfit)))
  expect_true(pr$haveboth)

  .ctReportCheckCalibration(pr)
  # This fit has real ridges (DRIFT against CINT), so the check above has
  # something to bite on. If it ever stops being true of the fixture, the
  # calibration assertion has quietly lost its teeth.
  expect_gt(max(pr$table$ridge, na.rm = TRUE), 2)
  # And on the near-independent parameters, where the slice is the whole
  # story, the quadratic should be close to exact.
  clean <- pr$table[is.finite(pr$table$ridge) & pr$table$ridge < 1.2, ]
  expect_gt(nrow(clean), 2)
  expect_equal(clean$drop_lo, rep(pr$expected, nrow(clean)), tolerance = 0.25)
  expect_equal(clean$drop_hi, rep(pr$expected, nrow(clean)), tolerance = 0.25)

  # A converged fit peaks at the estimate in every coordinate.
  expect_true(all(abs(pr$table$peak) <= 0.5))
  expect_false(any(grepl("off-peak", pr$table$flag)))

  # And the ridge column is the ratio it claims to be.
  ses <- ctsem:::.ctReportRawSE(ctstantestfit)
  expect_equal(pr$table$ridge, ses$marginal / ses$conditional, tolerance = 1e-8)
  expect_true(all(ses$marginal >= ses$conditional * (1 - 1e-8)))
})

test_that("the profile degrades to the reported SE when the fit carries no Hessian", {
  skip_on_cran()

  # A fit optimised with uncertainty off has draws or a covariance but no
  # Hessian, so there is no conditional width and no ridge to report. It must
  # still produce a profile, and say which of its columns stopped meaning
  # what the header claims.
  f <- ctstantestfit
  f$stanfit$uncertainty$hessian <- NULL
  ses <- ctsem:::.ctReportRawSE(f)
  expect_null(ses$conditional)
  expect_false(is.null(ses$marginal))

  pr <- ctsem:::.ctReportProfile(f, npoints = 5L)
  expect_false(pr$haveboth)
  expect_true(all(is.na(pr$table$ridge)))
  expect_equal(pr$table$se_slice, pr$table$se_report)

  d <- .ctReportFolder()
  res <- suppressWarnings(suppressMessages(
    ctReport(f, folder = d, components = "profile", profilepoints = 5,
      verbose = FALSE)))
  expect_identical(res$components$profile$status, "ok")
  txt <- readLines(file.path(d, "04-loglik-profile.txt"))
  expect_true(any(grepl("No Hessian on this fit", txt, fixed = TRUE)))
  .ctReportCheckFolder(res)
})

test_that("ctReport removes only its own files from a reused folder", {
  skip_on_cran()

  d <- .ctReportFolder()
  writeLines("keep me", file.path(d, "notes.txt"))
  writeLines("stale", file.path(d, "99-stale.txt"))
  res <- suppressWarnings(suppressMessages(
    ctReport(ctstantestfit, folder = d, components = "summary", verbose = FALSE)))

  expect_true(file.exists(file.path(d, "notes.txt")))
  expect_false(file.exists(file.path(d, "99-stale.txt")))
  .ctReportCheckFolder(res)
})

test_that("the registry, the component functions and the filename prefixes agree", {
  reg <- ctsem:::.ctReportComponents()
  expect_setequal(names(reg), names(ctsem:::.ctReportDo))
  expect_setequal(unique(unname(reg)), c("quick", "default", "all"))

  # A component with no interpretation entry would write its heading into the
  # index and then two blank lines -- paste0(NULL) is character(0), which
  # writes nothing and reports no error. Both lines must exist and be prose.
  interp <- ctsem:::.ctReportInterpretation()
  expect_setequal(names(interp), names(reg))
  for (nm in names(interp)) {
    expect_length(interp[[nm]], 2)
    expect_true(all(nzchar(interp[[nm]])), info = nm)
    expect_true(all(nchar(interp[[nm]]) > 25), info = nm)
  }

  sets <- ctsem:::.ctReportSets()
  expect_true(all(sets$quick %in% sets$default))
  expect_true(all(sets$default %in% sets$all))
  expect_setequal(sets$all, names(reg))
  # 'default' and 'all' resolve in registry order, not set order, so a
  # component moved between sets does not silently move in the reading order.
  expect_identical(ctsem:::.ctReportResolve("all"), names(reg))
  expect_identical(ctsem:::.ctReportResolve(c("residuals", "summary")),
    c("summary", "residuals"))

  # The numeric prefixes must ascend in registry order: they are the reading
  # order, and they live in the component bodies rather than the registry.
  skip_on_cran()
  d <- .ctReportFolder()
  res <- suppressWarnings(suppressMessages(
    ctReport(ctstantestfit, folder = d, components = "quick", nsamples = 10,
      profilepoints = 3, verbose = FALSE)))
  prefix <- vapply(names(reg)[names(reg) %in% names(res$components)],
    function(nm) {
      f <- res$components[[nm]]$files
      if (is.null(f)) NA_character_ else substr(sort(f)[1], 1, 2)
    }, character(1))
  prefix <- prefix[!is.na(prefix)]
  expect_gt(length(prefix), 5)
  expect_identical(prefix, prefix[order(prefix)])
  expect_false(any(duplicated(prefix)))
})

test_that("ctReport rejects a non-fit and an unknown component", {
  expect_error(ctReport(1:10), "ctStanFit or ctJuliaFit")
  expect_error(ctReport(ctstantestfit, components = "nosuchthing"),
    "Unknown report component")
})

test_that("ctReport writes the same folder for a julia fit, and its own diagnostics", {
  skip_on_cran()
  skip_without_julia()

  set.seed(3)
  n <- 10; tp <- 5
  dat <- do.call(rbind, lapply(1:n, function(i) data.frame(
    id = i, time = 0:(tp - 1), Y1 = stats::rnorm(tp), Y2 = stats::rnorm(tp),
    TI1 = rep(stats::rnorm(1), tp))))
  m <- suppressWarnings(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), MANIFESTVAR = diag(.2, 2), MANIFESTMEANS = matrix(0, 2, 1),
    n.TIpred = 1, TIpredNames = "TI1"))
  # priors = FALSE: at HEAD this model shape (2 latents, 1 covariate, default
  # tipred effects) fails .ctBackendPriorSpec() with "the Stan layout accounts
  # for 21 of 24 free parameters". Nothing to do with ctReport, and the report
  # needs a fit rather than a prior.
  fit <- suppressWarnings(suppressMessages(ctFit(dat, m, backend = "julia",
    cores = 1, priors = FALSE)))

  d <- .ctReportFolder()
  res <- suppressWarnings(suppressMessages(
    ctReport(fit, folder = d, components = c("summary", "parmatrices",
      "identification", "profile", "discretepars", "predictions", "residuals",
      "covcheck", "postpred", "tipredeffects"),
      nsamples = 10, profilepoints = 5, verbose = FALSE)))

  .ctReportCheckFolder(res)
  expect_true(all(vapply(res$components, function(x) x$status, character(1)) == "ok"))

  index <- readLines(file.path(d, "00-index.md"))
  expect_true(any(grepl("Backend: julia", index, fixed = TRUE)))

  # The julia path is the one that carries an identifiability decomposition at
  # fit time, so the identification file must actually report it rather than
  # falling through to the stan branch.
  ident <- readLines(file.path(d, "03-identification.txt"))
  expect_true(any(grepl("condition number", ident, fixed = TRUE)))
  expect_true(any(grepl("weak directions", ident, fixed = TRUE)))
  expect_true(any(grepl("converged:", ident, fixed = TRUE)))

  # Same calibration check as for stan, on the other backend's Hessian and
  # its own log-probability function.
  pr <- ctsem:::.ctReportProfile(fit, npoints = 5L)
  expect_true(pr$haveboth)
  .ctReportCheckCalibration(pr)
})
