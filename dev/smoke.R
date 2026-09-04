# Smoke test. Run this before merging ANYTHING that touches the julia engine,
# the R backend glue, or the build. It is deliberately tiny: it proves a fit
# completes and returns, not that the statistics are right.
#
#   Rscript dev/smoke.R
#
# Why it exists. On 2026-09-04 a change was merged that returned an empty vector
# from the engine in the normal case. A zero-length vector deadlocks the
# JuliaConnectoR bridge in both directions, so every healthy julia fit hung the
# moment its result crossed back to R. The engine's own 1400-assertion suite
# passed, because it runs inside julia and never crosses the bridge, and the
# change was merged on inspection after two attempts at a live fit were lost to
# machine contention. This file is the thing that would have caught it in under
# a minute.
#
# So: engine tests are not a substitute for one real fit through R.

Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("CTSEM_JULIA_AGREE"))) Sys.setenv(CTSEM_JULIA_AGREE = "yes")

t_all <- Sys.time()
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
cat("loaded\n"); flush.console()

set.seed(1)
n <- 5; tp <- 5
dat <- data.frame(id = rep(1:n, each = tp), time = rep(0:(tp - 1), n),
  Y1 = rnorm(n * tp))

ok <- TRUE
step <- function(label, expr, check) {
  t0 <- Sys.time()
  r <- tryCatch(eval(expr), error = function(e) e)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (inherits(r, "error")) {
    cat(sprintf("FAIL %-34s %5.1fs  %s\n", label, el, conditionMessage(r)))
    ok <<- FALSE
    return(invisible(NULL))
  }
  good <- isTRUE(try(check(r), silent = TRUE))
  cat(sprintf("%-4s %-34s %5.1fs\n", if (good) "ok" else "FAIL", label, el))
  if (!good) ok <<- FALSE
  invisible(r)
}

m <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, manifestNames = "Y1",
  latentNames = "eta1", LAMBDA = matrix(1, 1, 1),
  MANIFESTVAR = matrix(0.2, 1, 1))

# A fit that RETURNS is the point. The bridge deadlock this file exists for
# produced a fit that computed correctly and never came back.
fj <- step("julia optimise", quote(suppressMessages(suppressWarnings(
  ctFit(dat, m, backend = "julia", cores = 1, priors = FALSE)))),
  function(r) is.finite(r$estimate$loglik))

# Individually varying parameters take the augmented-state path, which is where
# most of the model-writer machinery lives and where several 2026-09 defects sat.
mi <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, manifestNames = "Y1",
  latentNames = "eta1", LAMBDA = matrix(1, 1, 1),
  MANIFESTVAR = matrix(0.2, 1, 1), CINT = "cint|param|TRUE",
  MANIFESTMEANS = 0)
step("julia optimise, indvarying", quote(suppressMessages(suppressWarnings(
  ctFit(dat, mi, backend = "julia", cores = 1, priors = TRUE,
    intoverpop = "augmented")))),
  function(r) is.finite(r$estimate$loglik))

# Post-fit accessors cross the bridge separately from the fit itself.
if (!is.null(fj)) {
  step("summary", quote(suppressMessages(summary(fj))),
    function(r) !is.null(r$popmeans))
  step("ctKalman", quote(suppressMessages(ctKalman(fj))),
    function(r) nrow(r) > 0)
}

# Stan, so a build or model-writer change that breaks it is caught here too.
step("stan optimise", quote(suppressMessages(suppressWarnings(
  ctFit(dat, m, backend = "stan", cores = 1, priors = FALSE,
    optimcontrol = list(carefulfit = FALSE))))),
  function(r) is.finite(summary(r)$loglik))

cat(sprintf("\n%s in %.0fs total\n", if (ok) "SMOKE PASSED" else "SMOKE FAILED",
  as.numeric(difftime(Sys.time(), t_all, units = "secs"))))
if (!ok) quit(status = 1)
