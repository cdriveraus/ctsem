# Workspace slots are handed down as a band, and every route budgets its own.
#
# The Laplace unit functions can divide a unit's members across tasks, and each
# task needs its own adjoint workspace. Those come from `laplace.workspaces`,
# which is sized by whichever route is driving: the fit objective sizes it to
# chunks times member width, `quadrature.jl` to its chunk count, and
# `sample_density.jl` gives each chain a contiguous block so no two chains can
# reach the same workspace.
#
# So a slot may only ever be taken from the band the caller hands down. An
# earlier version derived slots from a global chunk count instead, which the
# fit tuner *pins* for the rest of the session -- so a fit followed by a
# `ctLaplaceCheck()` asked for slots the quadrature route had never budgeted.
# Nothing in the suite caught it, because the quadrature tests fit on one core
# and so never set a width above one. Hence this file.

.band_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), T0MEANS = matrix(paste0("t0_", 1:2), 2, 1),
    T0VAR = diag(2), MANIFESTMEANS = matrix(paste0("mm_", 1:2), 2, 1),
    MANIFESTVAR = "diag", DRIFT = "auto", CINT = matrix(0, 2, 1),
    DIFFUSION = "auto", id = c("subject", "grp"), time = "time",
    tipredDefault = FALSE)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying_grp <- FALSE
  m$pars$indvarying[m$pars$matrix %in% c("T0MEANS", "MANIFESTMEANS")] <- TRUE
  m$pars$indvarying_grp[m$pars$matrix == "MANIFESTMEANS"] <- TRUE
  m
}

.band_data <- function(ngroup = 3L, nperson = 12L) {
  set.seed(4)
  d <- expand.grid(time = c(0, .2, .45, .7), person = seq_len(nperson),
    grp = seq_len(ngroup))
  d <- d[order(d$grp, d$person, d$time), ]
  d$subject <- paste0(d$grp, "_", d$person)
  d$Y1 <- rnorm(nrow(d))
  d$Y2 <- rnorm(nrow(d))
  d
}

test_that("a fit that pins a member width leaves other routes working", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctFit(.band_data(), .band_model(),
    backend = "julia", intoverpop = "laplace", cores = 4, verbose = 0,
    optimcontrol = list(estonly = TRUE, maxiter = 3))))
  expect_true(is.finite(as.numeric(fit$estimate$loglik)))

  # The fit's tuner pins both numbers for the session. This is the state the
  # quadrature route below then inherits, and the whole point of the test.
  width <- as.integer(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM._LAPLACE_MEMBER_WIDTH[]"))
  expect_gte(width, 1L)

  # The quadrature route sizes `laplace.workspaces` to its own chunk count and
  # budgets no member tasks, so it must run serially whatever the fit pinned.
  # Before the band was passed explicitly this either indexed past the vector
  # or took a slot another chunk was using.
  chk <- suppressWarnings(suppressMessages(
    ctLaplaceCheck(fit, nodes = 3L, correction = FALSE)))
  expect_s3_class(chk, "ctLaplaceCheck")
  expect_true(all(is.finite(unlist(chk[vapply(chk, is.numeric, logical(1))]))))
})

test_that("asking for slots the caller never budgeted is refused, loudly", {
  skip_without_julia()
  # The guard exists because the failure it replaces is a `BoundsError` inside a
  # spawned task, which `.ctBackendLpgFunc` turns into a large finite penalty --
  # the optimiser then sees a bad point rather than a bug.
  msg <- tryCatch({
    JuliaConnectoR::juliaEval(
      "let lp = nothing
         try
           ContinuousTimeSEM._laplace_check_band((workspaces = [Dict{Any,Any}()],), 1, 4)
           \"no error\"
         catch e
           sprint(showerror, e)
         end
       end")
  }, error = function(e) conditionMessage(e))
  expect_match(as.character(msg), "slots 1\\.\\.4 asked for, 1 exist",
    fixed = FALSE)
})
