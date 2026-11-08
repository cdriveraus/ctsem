# anova() and ctChisqTest() on stan and julia fits, and across backends.
#
# Both read each fit's logLik() and parameter count through .ctAnova()
# (R/ctCompare.R), so a julia fit, a stan fit or one of each compare alike.
# ctChisqTest() returns the p-value of anova()'s second row, a plain number as
# in 3.11.1.

.chisq_stan_data <- function() {
  data.frame(id = 1, time = seq_along(sunspot.year[1:60]), Y1 = sunspot.year[1:60])
}

test_that("anova() and ctChisqTest() compare two nested stan fits", {
  skip_on_cran()

  dat <- .chisq_stan_data()
  m1 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0)
  m2 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0, DRIFT = .9)
  f1 <- suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0, backend = 'stan'))
  f2 <- suppressMessages(ctFit(dat, m2, cores = 1, verbose = 0, backend = 'stan'))

  p <- ctChisqTest(f1, f2)
  expect_true(is.numeric(p))
  expect_null(attributes(p))
  expect_true(!is.na(p))
  expect_true(p >= 0 && p <= 1)

  # The comparison is between models, not between the argument order they
  # were supplied in.
  expect_equal(p, ctChisqTest(f2, f1))

  tab <- anova(f1, f2)
  expect_s3_class(tab, "anova")
  expect_equal(rownames(tab), c("f2", "f1"))  # f2 fixes the drift
  expect_equal(tab$Df, c(NA, 1))
  expect_equal(tab$Chisq[2], 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f2))))
  expect_equal(tab[["Pr(>Chisq)"]][2], p)
})

.chisq_julia_model <- function(fixedDrift = NA) {
  drift <- if (is.na(fixedDrift)) "drift" else fixedDrift
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix(drift, 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.chisq_julia_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8; diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(20), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
  data
}

test_that("anova() compares nested julia fits, and a julia fit against a stan fit", {
  skip_without_julia()

  dat <- .chisq_julia_data()
  f_fixed <- suppressMessages(ctFit(dat, .chisq_julia_model(-0.5), backend = "julia", verbose = 0))
  f_free <- suppressMessages(ctFit(dat, .chisq_julia_model(), backend = "julia", verbose = 0))
  expect_s3_class(f_fixed, "ctJuliaFit")
  expect_s3_class(f_free, "ctJuliaFit")

  # The table: ordered by parameters whatever the argument order, rows named
  # by the arguments, the test on the second row.
  tab <- anova(f_free, f_fixed)
  expect_s3_class(tab, "anova")
  expect_equal(rownames(tab), c("f_fixed", "f_free"))
  expect_equal(names(tab), c("npar", "AIC", "BIC", "logLik", "-2LL", "Chisq",
    "Df", "Pr(>Chisq)"))
  # f_free has exactly one more free parameter (drift) than f_fixed.
  expect_equal(tab$npar, c(length(ctsem:::.ctFitRawEstimate(f_fixed)),
    length(ctsem:::.ctFitRawEstimate(f_free))))
  expect_equal(tab$Df, c(NA, 1))
  expect_equal(tab$logLik, c(as.numeric(logLik(f_fixed)), as.numeric(logLik(f_free))))
  expect_equal(tab$`-2LL`, -2 * tab$logLik)
  expect_equal(tab$AIC, c(AIC(f_fixed), AIC(f_free)))
  expect_equal(tab$BIC, c(BIC(f_fixed), BIC(f_free)))
  expect_equal(tab$Chisq, c(NA, 2 * diff(tab$logLik)))
  expect_equal(tab[["Pr(>Chisq)"]],
    c(NA, stats::pchisq(tab$Chisq[2], 1, lower.tail = FALSE)))
  out <- capture.output(print(tab))
  expect_match(out[1], "Likelihood-ratio tests", fixed = TRUE)
  expect_true(any(grepl("^f_free ", out)))

  # ctChisqTest() is that p-value, a plain number in either order.
  p <- ctChisqTest(f_fixed, f_free)
  expect_identical(p, tab[["Pr(>Chisq)"]][2])
  expect_identical(p, ctChisqTest(f_free, f_fixed))

  # Named arguments name the rows; fits handed over by value are numbered
  # rather than deparsed.
  expect_equal(rownames(anova(free = f_free, fixed = f_fixed)), c("fixed", "free"))
  expect_equal(rownames(do.call(anova, list(f_free, f_fixed))), c("fit2", "fit1"))

  # What anova() refuses or flags, on altered copies of the same fits.
  expect_error(anova(f_free), "two or more", fixed = TRUE)
  other <- f_free
  other$standata$ndatapoints <- other$standata$ndatapoints + 1L
  expect_error(anova(f_fixed, other), "same data", fixed = TRUE)
  worse <- f_free
  worse$estimate$loglik <- as.numeric(logLik(f_fixed)) - 1
  expect_warning(anova(f_fixed, worse), "lower log likelihood", fixed = TRUE)
  priored <- f_free
  priored$args$resolved$priors <- TRUE
  expect_true(any(grepl("priored: fitted with priors",
    attr(anova(f_fixed, priored), "heading"), fixed = TRUE)))
  expect_false(any(grepl("priors", attr(tab, "heading"), fixed = TRUE)))

  # Cross-backend: the same model fit with backend='stan' compares with the
  # julia fit like any other.
  f_free_stan <- suppressMessages(ctFit(dat, .chisq_julia_model(), cores = 1, verbose = 0, backend = 'stan'))
  cross <- anova(f_fixed, f_free_stan)
  expect_equal(cross$Df, c(NA, 1))
  expect_true(is.finite(cross$Chisq[2]))
  expect_identical(ctChisqTest(f_fixed, f_free_stan), cross[["Pr(>Chisq)"]][2])
})

test_that("anova() and ctChisqTest() refuse a sampled stan fit, naming ctLOO()", {
  skip_on_cran()
  # A sampled fit has a posterior, not a single log likelihood. Before the
  # check, the comparison failed with "arguments imply differing number of
  # rows", which names neither the sampled fit nor why. See
  # review/J7-sampled-fit-support.md.
  dat <- .chisq_stan_data()
  m1 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0)
  f_opt <- suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0, backend = 'stan'))
  f_samp <- suppressWarnings(suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0,
    optimize = FALSE, chains = 1, iter = 60, backend = 'stan')))
  expect_gt(length(f_samp$stanfit$stanfit@sim), 0)  # genuinely sampled

  err <- tryCatch({
    ctChisqTest(f_opt, f_samp)
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_false(is.na(err))
  expect_match(err, "sampled", fixed = TRUE)
  expect_match(err, "ctLOO", fixed = TRUE)

  # Order shouldn't matter for the refusal either, and anova() names the fit.
  expect_error(ctChisqTest(f_samp, f_opt), "sampled", fixed = TRUE)
  expect_error(anova(f_opt, f_samp), "f_samp is sampled", fixed = TRUE)

  # Two optimized fits with the same parameters: nothing to test, so no
  # p-value, and no error.
  expect_true(is.na(ctChisqTest(f_opt, f_opt)))
  same <- anova(f_opt, f_opt)
  expect_equal(rownames(same), c("fit1", "fit2"))
  expect_equal(same$Df, c(NA, 0))
})
