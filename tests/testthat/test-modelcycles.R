skip_on_cran()
{  # body of the guard this replaced; indentation unchanged

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

  test_that("a t0 matrix referencing a latent state is accepted and evaluates correctly", {
    # T0MEANS[1,1] reads eta2 (state 2) with a multiplier, T0MEANS[2,1] is the
    # ordinary free parameter that state 2 itself resolves to at t0. This is
    # acyclic -- state 2 does not depend back on state 1 -- and is the standard
    # idiom for a stable latent intercept, ctsem's state-space replacement for
    # MANIFESTTRAITVAR. It must build, not be refused.
    mod <- suppressMessages(m2(T0MEANS = c("0.5 * eta2", "t0m2")))
    expect_s3_class(mod, "ctStanModel")

    # Stan side: the referenced-state cell gets its own `stateref` column
    # (matsetup column 10) rather than the `param` column, with `param` left
    # at 0 so nothing can mistake it for a parameter reference, and `when = 1`
    # so it materialises in the t0 pass, after the state it reads is known.
    dat <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2), Y1 = 0, Y2 = 0)
    prepared <- suppressMessages(ctFit(dat, mod, backend = "stan", fit = FALSE))
    ms <- prepared$setup$matsetup
    row1 <- ms[ms$matrix == 1 & ms$row == 1 & ms$col == 1, ]
    row2 <- ms[ms$matrix == 1 & ms$row == 2 & ms$col == 1, ]
    expect_equal(row1$param, 0L)
    expect_equal(row1$stateref, 2L)
    expect_equal(row1$when, 1L)
    expect_true(row2$param > 0L)
    expect_equal(row2$stateref, 0L)

    # Julia side: T0MEANS[1,1] resolves, at model-build time, to the same free
    # parameter as T0MEANS[2,1] with its own transform composed on top --
    # not to whichever parameter happens to share state 2's index.
    table <- ctsem:::.ctJuliaParameterTable(mod)
    t1 <- table[table$matrix == "T0MEANS" & table$row == 1 & table$col == 1, ]
    t2 <- table[table$matrix == "T0MEANS" & table$row == 2 & table$col == 1, ]
    expect_false(is.na(t1$parnumber))
    expect_equal(t1$parnumber, t2$parnumber)
    eval_transform <- function(text, value) {
      eval(parse(text = gsub("param\\[\\d+\\]", value, text)))
    }
    v2 <- eval_transform(t2$transform, 0.37)
    v1 <- eval_transform(t1$transform, 0.37)
    expect_equal(v1, 0.5 * v2)
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
