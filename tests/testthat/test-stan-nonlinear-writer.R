# The stan model writer simplifies a state-dependent model's calculations with
# Deriv, which inlines a function it has no rule for. log1p_exp's body became
# overflow-free in 4.0.0 (pmax), Deriv failed on it with "invalid .Internal()
# argument", and no stan fit of a model with log1p_exp in a state-dependent
# cell could be prepared. Fit-free: prepared and checked by stanc, not compiled.

test_that("a stan model with log1p_exp in a state-dependent drift is written", {
  ln <- c("eta1", "eta2")
  DRIFT <- matrix("0", 2, 2)
  for (i in 1:2) DRIFT[i, i] <- paste0("-log1p_exp(dself", i, " + dbylevel * ",
    ln[i %% 2 + 1], ")")
  model <- suppressMessages(ctModel(LAMBDA = diag(2), manifestNames = c("Y1", "Y2"),
    latentNames = ln, DRIFT = DRIFT, PARS = c("dself1", "dself2", "dbylevel"),
    type = "ct", silent = TRUE))
  prepared <- suppressMessages(ctFit(ctstantestdat, model, backend = "stan",
    fit = FALSE))
  expect_match(prepared$stanmodeltext, "log1p_exp", fixed = TRUE)
  parsed <- rstan::stanc(model_code = prepared$stanmodeltext, verbose = FALSE)
  expect_true(isTRUE(parsed$status))
})
