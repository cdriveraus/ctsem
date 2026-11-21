# Rows where nothing is observed and no TD predictor acts leave the likelihood
# where the dynamics are linear and the same at every row: the prediction
# across such a row is the prediction over the whole interval
# (`_ctsem_rows_mergeable`, `_ctsem_likelihood_rows` in ctsem_backend.jl).
# The objective without them must be the objective with every row, and what is
# reported by row must still have every row. Where rows do not merge --
# discrete time, where each row is a step, and dynamics that read the state --
# every row stays. `ctsem_set_skip_empty_rows!(false)` is the oracle.

.emptyrows_data <- function(type = "ct") {
  set.seed(3)
  n <- 12
  rows <- 22
  out <- do.call(rbind, lapply(seq_len(n), function(i) {
    time <- if (type == "dt") seq_len(rows) else cumsum(stats::runif(rows, .2, 1.5))
    data.frame(id = i, time = time, Y1 = stats::rnorm(rows), Y2 = stats::rnorm(rows),
      td = 0)
  }))
  empty <- sample(nrow(out), round(.4 * nrow(out)))
  out[empty, c("Y1", "Y2")] <- NA
  out$Y2[sample(nrow(out), 30)] <- NA
  # A TD input on empty rows and on observed ones.
  out$td[c(sample(empty, 6), sample(setdiff(seq_len(nrow(out)), empty), 6))] <- 1
  # Leading and trailing empty rows.
  for (i in 1:3) {
    r <- which(out$id == i)
    out[r[c(1:2, rows - 1, rows)], c("Y1", "Y2")] <- NA
    out$td[r[c(1:2, rows - 1, rows)]] <- 0
  }
  out
}

# Value, gradient and per-row log likelihoods with rows skipped or not.
.emptyrows_evaluate <- function(asmodel, raw, skip) {
  ctsem:::.ctJuliaModule()
  old <- ctsem:::.ctJuliaEval(sprintf(
    "ContinuousTimeSEM.ctsem_set_skip_empty_rows!(%s)", if (skip) "true" else "false"))
  on.exit({
    ctsem:::.ctJuliaEval(sprintf("ContinuousTimeSEM.ctsem_set_skip_empty_rows!(%s)",
      if (isTRUE(old)) "true" else "false"))
    ctsem:::.ctJuliaForgetObjective(asmodel)
  })
  ctsem:::.ctJuliaForgetObjective(asmodel)
  out <- ctJuliaEvaluate(asmodel, raw, gradient = TRUE)
  out$llrow <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)$llrow
  out
}

.emptyrows_compare <- function(model, data, tolerance = 1e-9, fixed = NULL) {
  spec <- suppressWarnings(suppressMessages(ctFit(data, model, fit = FALSE)))
  npar <- ctsem:::.ctBackendNpar(spec)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  # Raw values by parameter name, where a random one could make the process
  # explode.
  for (p in names(fixed)) {
    raw[spec$parameter_table$parnumber[match(p, spec$parameter_table$param)]] <- fixed[[p]]
  }
  asmodel <- ctsem:::.ctBackendAsModel(spec)
  on <- .emptyrows_evaluate(asmodel, raw, TRUE)
  off <- .emptyrows_evaluate(asmodel, raw, FALSE)
  expect_true(is.finite(off$value))
  expect_equal(on$value, off$value, tolerance = tolerance)
  expect_equal(unlist(on$gradient), unlist(off$gradient), tolerance = 1e-6)
  # Every row reported, and the rows sum to the likelihood.
  expect_length(on$llrow, nrow(data))
  expect_equal(on$llrow, off$llrow, tolerance = 1e-8)
  invisible(list(on = on, off = off))
}

test_that("rows with nothing to filter leave the likelihood unchanged", {
  skip_without_julia()
  data <- .emptyrows_data()
  model <- ctModel(type = "ct", LAMBDA = diag(2), n.TDpred = 1, TDpredNames = "td",
    manifestNames = c("Y1", "Y2"))
  .emptyrows_compare(model, data)
})

test_that("discrete time and state-dependent dynamics filter every row", {
  skip_without_julia()
  # A row of a discrete-time model is a step of the process, so leaving one
  # out would change the model.
  dt <- ctModel(type = "dt", LAMBDA = diag(2), n.TDpred = 1, TDpredNames = "td",
    manifestNames = c("Y1", "Y2"))
  .emptyrows_compare(dt, .emptyrows_data("dt"))
  # Where the drift reads the state, the filter linearises at every row.
  nl <- suppressWarnings(ctModel(type = "ct", LAMBDA = diag(1), manifestNames = "Y1",
    PARS = matrix("drift", 1, 1), DRIFT = matrix("PARS[1,1] * (1 + 0.5 * eta1)", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(2, 1, 1)))
  .emptyrows_compare(nl, .emptyrows_data()[, c("id", "time", "Y1")], tolerance = 1e-12,
    fixed = list(drift = -0.5))
})
