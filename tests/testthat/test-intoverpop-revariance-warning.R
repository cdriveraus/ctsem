# Individual variation on a variance cell (DIFFUSION, MANIFESTVAR) is only
# partially identified under `intoverpop='augmented'`: the random effect is a
# static latent state that never reaches the observation mean, so the Kalman
# update cannot move it and only its covariance with the mean-affecting effects
# is determined. `ctFit()` warns and points at `intoverpop='laplace'`. The check
# sits ahead of the backend dispatch because the filter is the same on both, so
# these run with `fit=FALSE` and need no fit.
#
# See review/RANDOMEFFECTS-partial-identification-2026-09-07.md.

.revar_data <- function() {
  set.seed(1)
  data.frame(id = rep(1:6, each = 4), time = rep(0:3, 6),
    y1 = stats::rnorm(24), y2 = stats::rnorm(24))
}

.revar_model <- function() {
  suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("y1", "y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), CINT = matrix(0, 2, 1), T0MEANS = matrix(0, 2, 1)))
}

# Returns the partial-identification warnings only. `expect_silent()` is no use
# here: ctFit emits unrelated messages and warnings on these models, and the
# question is whether this one fired.
.revar_warnings <- function(model, ...) {
  w <- character(0)
  withCallingHandlers(
    suppressMessages(ctFit(.revar_data(), model, fit = FALSE, verbose = 0, ...)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  grep("only partially identified", w, value = TRUE)
}

test_that("indvarying DIFFUSION warns under every spelling of augmented", {
  m <- .revar_model()
  m$pars$indvarying[m$pars$matrix %in% "DIFFUSION"] <- TRUE

  # 'auto' resolves to augmented when optimising with free indvarying pars,
  # TRUE means augmented, and 'augmented' says so. All three must warn.
  for (iop in list("auto", TRUE, "augmented")) {
    w <- .revar_warnings(m, intoverpop = iop)
    expect_length(w, 1)                       # once per fit, not once per par
    expect_match(w, "DIFFUSION or MANIFESTVAR")
    expect_match(w, "intoverpop='laplace'")
    expect_match(w, "diff_eta1")              # names the parameters found
    expect_match(w, "diff_eta2")
  }
})

test_that("indvarying MANIFESTVAR warns under augmented", {
  m <- .revar_model()
  m$pars$indvarying[m$pars$matrix %in% "MANIFESTVAR" &
      is.na(m$pars$value)] <- TRUE
  w <- .revar_warnings(m, intoverpop = "augmented")
  expect_length(w, 1)
  expect_match(w, "mvary1")
  expect_match(w, "mvary2")
})

test_that("indvarying DIFFUSION is silent under laplace", {
  skip_without_julia()
  m <- .revar_model()
  m$pars$indvarying[m$pars$matrix %in% "DIFFUSION"] <- TRUE
  expect_length(.revar_warnings(m, intoverpop = "laplace", backend = "julia"), 0)
})

test_that("mean-affecting indvarying parameters are silent under augmented", {
  m <- .revar_model()
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$matrix %in%
      c("MANIFESTMEANS", "CINT", "T0MEANS") & is.na(m$pars$value)] <- TRUE
  expect_true(any(m$pars$indvarying))
  expect_length(.revar_warnings(m, intoverpop = "augmented"), 0)
})

test_that("no indvarying parameters is silent", {
  m <- .revar_model()
  m$pars$indvarying <- FALSE
  expect_length(.revar_warnings(m, intoverpop = "auto"), 0)
})

test_that("indvarying MANIFESTVAR on a binary indicator is silent", {
  # ctFit fixes free MANIFESTVAR for a non-Gaussian indicator to a
  # deterministic calculation further down, so warning about it here would
  # contradict that message. The check excludes those rows.
  skip_without_julia()
  m <- .revar_model()
  m$manifesttype[] <- 1L
  m$pars$indvarying[m$pars$matrix %in% "MANIFESTVAR" &
      is.na(m$pars$value)] <- TRUE
  set.seed(2)
  d <- .revar_data()
  d$y1 <- stats::rbinom(24, 1, .5)
  d$y2 <- stats::rbinom(24, 1, .5)
  w <- character(0)
  withCallingHandlers(
    suppressMessages(ctFit(d, m, fit = FALSE, verbose = 0, backend = "julia")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(grep("only partially identified", w), 0)
})
