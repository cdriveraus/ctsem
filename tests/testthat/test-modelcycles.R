if(identical(Sys.getenv("NOT_CRAN"), "true")) {

  # The T0MEANS / state / PARS loop is expressible in a model specification and
  # the generated Stan program breaks it by evaluation order rather than by a
  # rule, so a cell caught in a real loop is read before anything writes it and
  # comes back as a plausible number. These pin the refusal, and just as
  # importantly they pin the shapes that must keep working: a PARS cell reading a
  # state is a supported pattern and only becomes a loop when the state's
  # T0MEANS entry depends back on that PARS cell.

  mn2 <- c("Y1", "Y2")
  ln2 <- c("eta1", "eta2")
  m2 <- function(...) ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = mn2, latentNames = ln2, LAMBDA = diag(2), ...)

  test_that("a t0 matrix referencing a latent state is refused", {
    expect_error(suppressMessages(m2(T0MEANS = c("0.5 * eta2", "t0m2"))),
      "t0 matrix cannot reference a latent state")
    # The message names the offending cell and the state, because the user's
    # only way back is to find the cell they wrote.
    expect_error(suppressMessages(m2(T0MEANS = c("0.5 * eta2", "t0m2"))),
      "T0MEANS\\[1,1\\]")
    expect_error(suppressMessages(m2(T0MEANS = c("0.5 * eta2", "t0m2"))), "eta2")
  })

  test_that("a T0MEANS to PARS to state loop is refused, and names the path", {
    e <- tryCatch(suppressMessages(m2(
      PARS = matrix("0.3 + 0.05 * eta1", 1, 1),
      T0MEANS = c("PARS[1,1]", "t0m2"))), error = function(e) e)
    expect_s3_class(e, "error")
    expect_match(conditionMessage(e), "Circular dependency")
    expect_match(conditionMessage(e), "PARS\\[1,1\\]")
    expect_match(conditionMessage(e), "T0MEANS\\[1,1\\]")
  })

  test_that("a loop through two PARS cells is refused", {
    expect_error(suppressMessages(m2(
      PARS = matrix(c("PARS[2,1]", "0.1 * eta1"), 2, 1),
      T0MEANS = c("PARS[1,1]", "t0m2"))),
      "Circular dependency")
  })

  test_that("state dependent specifications that are not loops still build", {
    # A PARS cell reading a state, referenced by DRIFT. The tested pattern.
    expect_s3_class(suppressMessages(m2(
      PARS = matrix("0.3 + 0.05 * eta1", 1, 1),
      DRIFT = matrix(c("PARS[1,1]", 0, 0, "d22"), 2, 2))), "ctStanModel")

    # Direct state references in the measurement and dynamics matrices.
    expect_s3_class(suppressMessages(ctModel(type = "ct", n.latent = 2,
      n.manifest = 2, manifestNames = mn2, latentNames = ln2,
      LAMBDA = matrix(c("lbystate * eta2 + 1", 0, 0, 1), 2, 2),
      DRIFT = matrix(c("-0.5 * eta2", 0, 0, "d22"), 2, 2))), "ctStanModel")

    # PARS as an ordinary free parameter with a transform.
    expect_s3_class(suppressMessages(m2(PARS = c("dr11|-log1p_exp(param)"),
      DRIFT = matrix(c("PARS[1,1]", 0, 0, "d22"), 2, 2))), "ctStanModel")

    # A T0MEANS transform specification is not a state reference.
    expect_s3_class(suppressMessages(m2(
      T0MEANS = c("t0m1", "t0m2|log1p_exp(param)"))), "ctStanModel")
  })

  test_that("a latent name is not matched inside a longer latent name", {
    # eta1 must not match inside eta10, or a ten latent model would report a
    # dependency on the wrong state.
    d <- diag(-0.5, 10)
    d[1, 2] <- "0.1 * eta10"
    expect_s3_class(suppressMessages(ctModel(type = "ct", n.latent = 10,
      n.manifest = 10, manifestNames = paste0("Y", 1:10),
      latentNames = paste0("eta", 1:10), LAMBDA = diag(10),
      DRIFT = d)), "ctStanModel")
  })

  test_that("the reference extractors read names and bracketed forms alike", {
    expect_equal(ctsem:::.ctCycleStateRefs("0.3 + 0.05 * eta2", ln2), 2L)
    expect_equal(ctsem:::.ctCycleStateRefs("state[2] * 2", ln2), 2L)
    expect_equal(ctsem:::.ctCycleStateRefs("eta1 + eta2", ln2), c(1L, 2L))
    expect_equal(ctsem:::.ctCycleStateRefs("d22", ln2), integer(0))
    # Only the first element of a "|" separated specification carries the value.
    expect_equal(ctsem:::.ctCycleStateRefs("t0m2|log1p_exp(param)", ln2),
      integer(0))
    expect_equal(ctsem:::.ctCycleParsRefs("PARS[1,1] * (1 + PARS[2,1])"),
      c("1,1", "2,1"))
    expect_equal(ctsem:::.ctCycleParsRefs("PARS[ 1 , 2 ]"), "1,2")
    expect_equal(ctsem:::.ctCycleParsRefs("d22"), character(0))
  })
}
