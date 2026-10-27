# The symbolic Jacobian ctJacobian() writes for a state-dependent DRIFT cell,
# for each function of the engine's cell vocabulary, against central
# differences of the drift itself. Fit-free. min and max once came out 0 here
# (Deriv read an escaped comma as part of a name), and inv stopped the model.

test_that("each vocabulary function's drift Jacobian matches its differences", {
  cells <- c("-exp(eta2)", "-log(2 + eta2^2)", "-sqrt(1 + eta2^2)", "-abs(eta2)",
    "-log1p_exp(eta2)", "-inv_logit(eta2)", "-tanh(eta2)", "-log1p(eta2^2)",
    "-expm1(eta2)", "-sin(eta2)", "-cos(eta2)", "-inv(2 + eta2^2)",
    "-min(eta2, 1)", "-max(eta2, -1)", "-min(eta2, 0.2, 2*eta2)")
  env <- new.env(parent = asNamespace("ctsem"))
  for (cell in cells) {
    DRIFT <- matrix(c(cell, "0", "0", "-1"), 2, 2, byrow = TRUE)
    m <- suppressMessages(ctModel(type = "ct", LAMBDA = diag(2),
      manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
      DRIFT = DRIFT, silent = TRUE))
    m$pars <- ctsem:::ctModelStatesAndPARS(m$pars, statenames = m$latentNames,
      tdprednames = m$TDpredNames)
    J <- ctsem:::ctJacobian(m, types = "JAx")$JAx
    drift <- gsub("eta2", "state[2]", cell, fixed = TRUE)
    f <- function(s) eval(str2lang(paste0("(", drift, ") * state[1]")),
      list(state = s), env)
    for (s2 in c(-1.7, -0.3, 0.6, 1.4)) {
      s <- c(0.8, s2)
      h <- 1e-6
      numeric <- vapply(1:2, function(j) {
        e <- replace(numeric(2), j, h)
        (f(s + e) - f(s - e)) / (2 * h)
      }, numeric(1))
      written <- vapply(J[1, ], function(x) eval(str2lang(x), list(state = s), env),
        numeric(1))
      expect_equal(unname(written), numeric, tolerance = 1e-6,
        info = paste(cell, "at state", s2))
    }
  }
})
