# Persons nested in studies, as in a meta-analysis of panels: two processes
# measured directly, each person with their own intercepts and a starting
# point, each study with its own drift. On the layer a study's unit holds every
# person's effects, so these are the wide nested units whose curvature must
# stay in its blocks. Shared by the mode-search, soft-quadrature and map tests
# (test-julia-laplace-modes.R, test-julia-sample.R).
.nested_panel_data <- function(nstudy = 2, nperson = 12, ntime = 6) {
  set.seed(5)
  A0 <- matrix(c(-0.6, 0.15, 0.1, -0.5), 2, 2)
  G <- t(chol(matrix(c(0.6, 0.1, 0.1, 0.5), 2, 2)))
  out <- list()
  for (s in seq_len(nstudy)) {
    A <- A0 + matrix(rnorm(4, 0, 0.08), 2, 2)
    for (p in seq_len(nperson)) {
      mm <- rnorm(2); x <- rnorm(2)
      times <- cumsum(c(0, runif(ntime - 1, 0.5, 1.5)))
      X <- matrix(NA, ntime, 2); X[1, ] <- x
      for (t in 2:ntime) {
        dt <- times[t] - times[t - 1]; n <- ceiling(dt / 0.05); h <- dt / n
        for (j in seq_len(n)) x <- x + A %*% x * h + G %*% rnorm(2, 0, sqrt(h))
        X[t, ] <- x
      }
      Y <- X + matrix(mm, ntime, 2, byrow = TRUE) + matrix(rnorm(2 * ntime, 0, 0.4), ntime, 2)
      out[[length(out) + 1]] <- data.frame(id = (s - 1) * nperson + p, study = s,
        time = times, Y1 = Y[, 1], Y2 = Y[, 2])
    }
  }
  do.call(rbind, out)
}

# The model for `.nested_panel_data()`: "intercepts", person-level random
# intercepts and initial means under a study-level random drift; "drift", each
# person's own random drift diagonal and intercepts under a study-level random
# drift diagonal, whose short records leave a unit's curvature soft.
.nested_panel_model <- function(kind = c("intercepts", "drift")) {
  kind <- match.arg(kind)
  m <- suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2),
    MANIFESTMEANS = matrix(c("mm1", "mm2"), 2, 1), CINT = matrix(0, 2, 1),
    id = c("id", "study"), time = "time", tipredDefault = FALSE))
  free <- is.na(m$pars$value)
  if (kind == "intercepts") {
    m$pars$indvarying <- ifelse(m$pars$matrix %in% c("MANIFESTMEANS", "T0MEANS") & free,
      "id", ifelse(m$pars$matrix == "DRIFT" & free, "study", FALSE))
  } else {
    diagonal <- m$pars$matrix == "DRIFT" & m$pars$row == m$pars$col
    m$pars$indvarying <- ifelse(m$pars$matrix == "MANIFESTMEANS" & free, "id",
      ifelse(diagonal & free, "id,study", FALSE))
  }
  m
}
