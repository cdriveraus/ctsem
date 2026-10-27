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
