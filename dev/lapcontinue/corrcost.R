# Where the quadrature correction's time goes, stage by stage, on the
# optimiser bench's Laplace cells (review/OPTIM-next-2026-09-27.md, job F).
#
#   Rscript dev/lapcontinue/corrcost.R <outdir> [case[:start] ...]
#
# The build: CC_LIB, an installed library (dev1: ship.sh's
# ~/dev/ctsemlib-bench-<label>), or CC_TREE, a source tree to load_all. Per
# case: one Laplace fit (laplace_correct = FALSE, cores = 1, the bench's
# start), saved to and on a rerun read from CC_FITS, so that two builds
# correct the same fit; then the default correction applied to it CC_REPS
# times (default 2), the first untimed so that nothing timed compiles. Writes
# <outdir>/<case>.rds (the records, estimates, per-subject values and the
# continuation's Hessian) and one line per timed run to <outdir>/corrcost.tsv:
# the record's wall seconds by stage and the engine's seconds by kind of call,
# with the counts. CC_VARIANT names an entry of `variants` below.
args <- commandArgs(TRUE)
OUT <- args[1]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NOT_CRAN = "true", CTSEM_JULIA_AGREE = "yes", JULIA_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
stamp <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
LIB <- Sys.getenv("CC_LIB"); TREE <- Sys.getenv("CC_TREE")
if (nzchar(LIB)) {
  .libPaths(c(LIB, .libPaths()))
  suppressMessages(library(ctsem))
  HARNESS <- Sys.getenv("CC_HARNESS", file.path(TREE, "dev/optimbench"))
} else {
  suppressMessages(devtools::load_all(TREE, compile = FALSE, quiet = TRUE))
  HARNESS <- file.path(TREE, "dev/optimbench")
}
ns <- asNamespace("ctsem")
for (nm in c(".ctLaplaceContinue", ".ctLaplaceContinueDefaults", ".ctBackendNpar",
  ".ctBackendRawParameterNames", ".ctJuliaOr")) assign(nm, get(nm, ns))
stamp("LOADED", if (nzchar(LIB)) LIB else TREE)
STORE <- Sys.getenv("BENCH_DATA", path.expand("~/dev/ctsem-bench-data"))
FITS <- Sys.getenv("CC_FITS", file.path(OUT, "fits"))
dir.create(FITS, showWarnings = FALSE, recursive = TRUE)
REPS <- as.integer(Sys.getenv("CC_REPS", "2"))
source(file.path(HARNESS, "cells.R"))
source(file.path(HARNESS, "starts.R"))
jbin <- Sys.getenv("BENCH_JULIA_BIN")
if (!nzchar(jbin)) {
  cand <- c(Sys.glob(file.path(Sys.getenv("HOME"), ".julia", "juliaup", "julia-1.12*", "bin")),
    Sys.glob(file.path(Sys.getenv("USERPROFILE"), ".julia", "juliaup", "julia-1.12*", "bin")))
  cand <- cand[file.exists(cand)]
  jbin <- if (length(cand)) sort(cand, decreasing = TRUE)[1L] else ""
}
if (nzchar(jbin)) suppressMessages(ctJuliaSetup(julia_bin = jbin, threads = 1L)) else
  suppressMessages(ctJuliaSetup(threads = 1L))
now <- function() proc.time()[["elapsed"]]

variants <- list(default = list())
VARIANT <- Sys.getenv("CC_VARIANT", "default")
control <- utils::modifyList(.ctLaplaceContinueDefaults, variants[[VARIANT]])

cases <- if (length(args) > 1) args[-1] else "ord4"
for (spec_arg in cases) {
  parts <- strsplit(spec_arg, ":", fixed = TRUE)[[1]]
  case <- parts[1]
  start <- if (length(parts) > 1) paste(parts[-1], collapse = ":") else "default:1"
  tag <- gsub("[^A-Za-z0-9_.-]", "_", paste0(case, ".", start))
  stamp("CASE", case, start)
  res <- tryCatch({
    P <- bench_problem(case)
    D <- bench_data(case, if (isTRUE(P$datafixed)) "cfg" else "1", STORE)
    d <- D$data; m <- P$model()
    spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = "laplace", fit = FALSE, cores = 1)))
    npar <- .ctBackendNpar(spec)
    rawnames <- suppressWarnings(.ctBackendRawParameterNames(list(model_spec = spec), npar))
    st <- bench_start(start, npar, rawnames)
    saved <- file.path(FITS, paste0(tag, "_fit.rds"))
    if (file.exists(saved)) {
      off <- readRDS(saved)$fit
      stamp("  fit read from", saved)
    } else {
      t0 <- now(); set.seed(st$seed)
      off <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
        intoverpop = "laplace", cores = 1, inits = st$inits,
        optimcontrol = list(laplace_correct = FALSE))))
      saveRDS(list(fit = off, secs = now() - t0, data_md5 = D$md5), saved)
      stamp(sprintf("  fit %.1fs", now() - t0))
    }
    runs <- list()
    for (r in seq_len(REPS)) {
      t0 <- now()
      out <- suppressWarnings(suppressMessages(.ctLaplaceContinue(off, control = control)))
      wall <- now() - t0
      cr <- out$laplace$correction
      runs[[r]] <- list(wall = wall, record = cr, estimate = as.numeric(out$estimate$raw),
        se = as.numeric(out$estimate$se), loglik = out$estimate$loglik,
        subject_loglik = out$estimate$subject_loglik)
      s <- cr$seconds; e <- cr$engine_seconds; ev <- cr$evaluations
      g <- function(v, k) if (!is.null(v) && k %in% names(v)) unname(v[[k]]) else NA_real_
      line <- data.frame(case = case, start = start, variant = VARIANT, rep = r,
        wall = wall, status = paste0(cr$status, "/", .ctJuliaOr(cr$continuation, "")),
        flagged = .ctJuliaOr(cr$flagged, NA), units = .ctJuliaOr(cr$units, NA),
        screen = g(s, "screen"), rounds = g(s, "rounds"), hessian = g(s, "hessian"),
        certification = g(s, "certification"), uncertainty = g(s, "uncertainty"),
        total = g(s, "total"), e_values = g(e, "values"), e_gradients = g(e, "gradients"),
        e_placements = g(e, "placements"), e_hessian = g(e, "hessian"),
        n_values = g(ev, "values"), n_gradients = g(ev, "gradients"),
        n_placements = g(ev, "placements"), member_values = g(ev, "member_values"),
        member_sweeps = g(ev, "member_sweeps"),
        loglik = .ctJuliaOr(out$estimate$loglik, NA),
        stringsAsFactors = FALSE)
      stamp(sprintf(paste0("  rep %d %.1fs %s: screen %.1f rounds %.1f hessian %.1f ",
        "cert %.1f unc %.1f | engine values %.1f (%s) gradients %.1f (%s) ",
        "placements %.1f (%s) hessian %.1f | loglik %.4f"), r, wall, line$status,
        line$screen, line$rounds, line$hessian, line$certification,
        line$uncertainty, line$e_values, line$n_values, line$e_gradients,
        line$n_gradients, line$e_placements, line$n_placements, line$e_hessian,
        line$loglik))
      tsv <- file.path(OUT, "corrcost.tsv")
      utils::write.table(line, tsv, sep = "\t", quote = FALSE, row.names = FALSE,
        append = file.exists(tsv), col.names = !file.exists(tsv))
    }
    saveRDS(list(case = case, start = start, variant = VARIANT, data_md5 = D$md5,
      laplace = as.numeric(off$estimate$raw), runs = runs),
      file.path(OUT, paste0(tag, ".", VARIANT, ".rds")))
    "ok"
  }, error = function(e) { stamp("ERROR", conditionMessage(e)); conditionMessage(e) })
  stamp("DONE", case, res)
}
stamp("ALLDONE")
