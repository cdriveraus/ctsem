# The quadrature correction's gain against what it predicts before any round
# runs, and its cost under the stopping rules before and after 2026-09-27, on
# the optimiser bench's default Laplace cells (one start each; the AnomAuth
# cells from their stored spurious maxima).
#
#   Rscript dev/lapcontinue/calibrate-cost.R <tree> <outdir> [case ...]
#
# Per case: one Laplace fit (laplace_correct = FALSE, cores = 1, the bench's
# start); the penalised exact log posterior at its optimum (the bench's
# reference: a trapezoid over the soft directions for units of dimension <= 3,
# importance sampling for the nested ones); then the correction applied post
# hoc twice to that same fit, timed, with the skip off in both:
#   old  stop_gain = 0: rounds stop on the fit's certification tolerance, or
#        on running out of rounds or attempts, as at juliaFit c74b500e;
#   new  the stopping rule on predicted gain in nats.
# Writes <outdir>/calibrate.tsv, one line per case, and <outdir>/<case>.rds.
# One Julia thread; the first case is run once untimed so that no timed
# correction pays compilation.
args <- commandArgs(TRUE)
TREE <- normalizePath(args[1], winslash = "/")
OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NOT_CRAN = "true", JULIA_NUM_THREADS = "1")
if (!nzchar(Sys.getenv("CTSEM_JULIA_AGREE"))) Sys.setenv(CTSEM_JULIA_AGREE = "yes")
suppressMessages(devtools::load_all(TREE, compile = FALSE, quiet = TRUE))
stamp <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
stamp("LOADED", TREE)
STORE <- Sys.getenv("BENCH_DATA", path.expand("~/dev/ctsem-bench-data"))
source(file.path(TREE, "dev/optimbench/cells.R"))
source(file.path(TREE, "dev/optimbench/starts.R"))
invisible(.ctJuliaModule())
JuliaConnectoR::juliaEval(sprintf('Core.eval(ContinuousTimeSEM, :(include(%s))); nothing',
  deparse(file.path(TREE, "dev/optimbench/reference.jl"))))
jc <- function(fn, ...) JuliaConnectoR::juliaCall(paste0("ContinuousTimeSEM.", fn), ...)
jg <- function(x) if (inherits(x, "JuliaProxy")) JuliaConnectoR::juliaGet(x) else x
jv <- function(x) .ctJuliaNumericVector(as.numeric(x))
now <- function() proc.time()[["elapsed"]]

CASES <- c(acnonlin = "default:1", anomS1 = "stored:anomS1_spurious",
  anomS2 = "stored:anomS2_spurious", cf_binary = "default:1", cf_mixed = "default:1",
  cf_ordinal = "default:1", gA1 = "default:1", gA14 = "default:1", gB2 = "default:1",
  gB8 = "default:1", gC2 = "default:1", gC8 = "default:1", gD1 = "default:1",
  gD3 = "default:1", gN1 = "default:1", gN3 = "default:1", jflat = "default:1",
  mvmix = "default:1", ordinal = "default:1")
which <- if (length(args) > 2) args[-(1:2)] else names(CASES)

reference <- function(fit, x) {
  obj <- .ctJuliaObjective(fit)
  ev <- jg(jc("ctsem_laplace_evaluate", obj, jv(x), gradient = FALSE))
  prior <- as.numeric(ev$value) - sum(as.numeric(ev$unit_loglik))
  PU <- as.matrix(jg(jc("bench_probe_units", obj, jv(x))))
  dmax <- ncol(PU) - 6L
  R <- if (dmax <= 3L) {
    vapply(seq_len(nrow(PU)), function(U) as.numeric(jg(jc("bench_probe_reference",
      obj, jv(x), as.integer(U), softcut = 3.5, nstiff = 9L, h = 0.1,
      halfwidth = 6.0))$value), 0)
  } else {
    vapply(seq_len(nrow(PU)), function(U) as.numeric(jg(jc("bench_probe_is", obj,
      jv(x), as.integer(U), n = 100000L, seed = 1L, inflate = 1.5))$value), 0)
  }
  sum(R) + prior
}
controls <- list(
  old = utils::modifyList(.ctLaplaceContinueDefaults, list(skip_gain = 0, stop_gain = 0)),
  new = utils::modifyList(.ctLaplaceContinueDefaults, list(skip_gain = 0)))
warmed <- FALSE
for (case in which) {
  f <- file.path(OUT, paste0(case, ".rds"))
  if (file.exists(f)) { stamp("have", case); next }
  stamp("CASE", case)
  res <- tryCatch({
    P <- bench_problem(case)
    D <- bench_data(case, if (isTRUE(P$datafixed)) "cfg" else "1", STORE)
    d <- D$data; m <- P$model()
    spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = "laplace", fit = FALSE, cores = 1)))
    npar <- .ctBackendNpar(spec)
    rawnames <- suppressWarnings(.ctBackendRawParameterNames(list(model_spec = spec), npar))
    st <- bench_start(CASES[[case]], npar, rawnames)
    t0 <- now(); set.seed(st$seed)
    off <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = "laplace", cores = 1, inits = st$inits,
      optimcontrol = list(laplace_correct = FALSE, finishsamples = 100))))
    fit_secs <- now() - t0
    x0 <- as.numeric(off$estimate$raw)
    exact0 <- reference(off, x0)
    if (!warmed) {
      for (ctl in controls) invisible(suppressWarnings(suppressMessages(
        .ctLaplaceContinue(off, control = ctl))))
      warmed <- TRUE
    }
    rec <- list(case = case, start = CASES[[case]], data_md5 = D$md5, npar = npar,
      fit_secs = fit_secs, laplace = x0, exact_laplace = exact0)
    line <- list(case = case, npar = npar, units = NA, screen = NA, predicted = NA,
      exact_laplace = exact0)
    for (nm in names(controls)) {
      t0 <- now()
      out <- suppressWarnings(suppressMessages(.ctLaplaceContinue(off,
        control = controls[[nm]])))
      secs <- now() - t0
      cr <- out$laplace$correction
      x <- as.numeric(out$estimate$raw)
      ex <- if (max(abs(x - x0)) == 0) exact0 else reference(off, x)
      rec[[nm]] <- list(record = cr, estimate = x, exact = ex, secs = secs)
      line$units <- cr$units; line$screen <- cr$screen
      if (identical(nm, "new")) line$predicted <- .ctJuliaOr(cr$predicted_gain, NA)
      line[[paste0("gain_", nm)]] <- ex - exact0
      line[[paste0("move_se_", nm)]] <- if (isTRUE(cr$applied))
        max(abs(cr$delta_se), na.rm = TRUE) else 0
      line[[paste0("secs_", nm)]] <- secs
      line[[paste0("grads_", nm)]] <- as.numeric(.ctJuliaOr(cr$evaluations[["gradients"]], NA))
      line[[paste0("status_", nm)]] <- paste0(cr$status, "/", .ctJuliaOr(cr$continuation, ""))
      stamp(sprintf("  %-4s %6.1fs %4s gradients, %s, gain %.4f nats, move %.3f se",
        nm, secs, format(line[[paste0("grads_", nm)]]), line[[paste0("status_", nm)]],
        ex - exact0, line[[paste0("move_se_", nm)]]))
    }
    stamp(sprintf("  screen %.4g, predicted gain %.4g", line$screen, line$predicted))
    utils::write.table(as.data.frame(line, stringsAsFactors = FALSE),
      file.path(OUT, "calibrate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
      append = file.exists(file.path(OUT, "calibrate.tsv")),
      col.names = !file.exists(file.path(OUT, "calibrate.tsv")))
    saveRDS(rec, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    "ok"
  }, error = function(e) { stamp("ERROR", conditionMessage(e)); conditionMessage(e) })
  stamp("DONE", case, res)
}
stamp("ALLDONE")
