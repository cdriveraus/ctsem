# ctModelTransformsToNum() reads a transform text as `offset + multiplier *
# g(param * meanscale + inneroffset)` exactly where one of the templates
# reproduces it (.ctTransformMatch), and leaves the least-squares search for the
# rest. Fit-free.

test_that("transform texts are read exactly, floors included", {
  cases <- list(
    list("param", 0, 1, 0, 1, 0),
    list("10 * param", 0, 10, 0, 1, 0),
    list("0.5 + 5 * param", 0, 5, 0.5, 1, 0),
    list("log1p_exp(param)", 1, 1, 0, 1, 0),
    list("-log1p_exp(-param)", 1, -1, 0, -1, 0),
    list("1e-10 + 5 * log1p_exp(2 * param)", 1, 5, 1e-10, 2, 0),
    list("10 * log1p_exp(2 * param) + 1e-10", 1, 10, 1e-10, 2, 0),
    list("-(1e-06 + 2 * log1p_exp(-(2 * param)))", 1, -2, -1e-6, -2, 0),
    list("log1p_exp(2*(param-1))*0.7", 1, 0.7, 0, 2, -2),
    list("-exp(param)", 2, -1, 0, 1, 0),
    list("1/(1+exp(-param))", 3, 1, 0, 1, 0),
    list("2/(1+exp(-param))-1", 3, 2, -1, 1, 0),
    list("exp(param)/(1+exp(param))", 3, 1, 0, 1, 0),
    list("(param)^3 * 0.5", 4, 0.5, 0, 1, 0),
    list("log1p(param)", 5, 1, 0, 1, 0))
  for (case in cases) {
    r <- ctsem:::.ctTransformMatch(case[[1]])
    expect_false(is.null(r), info = case[[1]])
    expect_equal(c(r$type, r$multiplier, r$offset, r$meanscale, r$inneroffset),
      unlist(case[-1]), tolerance = 0, info = case[[1]])
  }
})

test_that("an exact reading reproduces its text, and one no template fits stays text", {
  ctm <- list(pars = data.frame(matrix = "DRIFT", row = 1:4, col = 1L,
    param = paste0("p", 1:4), value = NA_real_,
    transform = c("1e-10 + 5 * log1p_exp(2 * param - 1)", "abs(param)",
      "2/(1+exp(-param))-1", "param^2"),
    indvarying = FALSE, sdscale = 1, stringsAsFactors = FALSE))
  pars <- ctModelTransformsToNum(ctm)$pars
  expect_equal(pars$transform[c(2, 4)], c("abs(param)", "param^2"))
  p <- seq(-3, 3, by = 0.25)
  for (r in c(1, 3)) {
    expected <- eval(str2lang(ctm$pars$transform[r]), list(param = p),
      asNamespace("ctsem"))
    read <- tform(p, as.integer(pars$transform[r]), pars$multiplier[r],
      pars$meanscale[r], pars$offset[r], pars$inneroffset[r])
    expect_equal(read, expected, tolerance = 1e-13, info = ctm$pars$transform[r])
  }
})

# A cell reading one state is first offered to the template reading above. A
# long one -- the Wiener first-passage series written as cells -- once stopped
# setup twice: the search's start grid, the count of the expression's numbers
# to the fourth (1458 GB at 200 terms), and an environment keyed by the cell's
# text, which R limits to 10000 bytes.
test_that("a long cell reading one state is set up as a calculation", {
  bsum <- function(x) if (length(x) == 1L) x else {
    h <- length(x) %/% 2L
    paste0("(", bsum(x[1:h]), " + ", bsum(x[-(1:h)]), ")")
  }
  digits <- paste(rep("1234567890", 25), collapse = "")
  long <- paste0("eta1 * ", bsum(sprintf("exp(-0.%d%s * eta1)", 1:40, digits)))
  expect_gt(nchar(long), 10000)
  read <- ctsem:::ctModelTransformsToNum(list(pars = data.frame(matrix = "PARS",
    row = 1, col = 1, param = "p1", value = NA_real_,
    transform = gsub("eta1", "param", long), indvarying = FALSE, sdscale = 1,
    stringsAsFactors = FALSE)))$pars
  expect_identical(read$transform, gsub("eta1", "param", long))
  m <- suppressMessages(ctModel(type = "ct", manifestNames = "y1",
    latentNames = "eta1", LAMBDA = matrix(1), MANIFESTMEANS = matrix(long),
    DRIFT = matrix("drift"), CINT = matrix(0), T0MEANS = matrix("t0m")))
  dat <- data.frame(id = rep(1:3, each = 4), time = rep(0:3, 3),
    y1 = c(0.1, -0.3, 0.2, 0.5, -0.1, 0.4, 0.0, -0.2, 0.3, 0.1, -0.4, 0.2))
  expect_s3_class(suppressMessages(ctFit(dat, m, fit = FALSE)), "ctFitModel")
})
