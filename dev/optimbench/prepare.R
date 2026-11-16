# Before a grid runs its cells in parallel: store every cell's data, check each
# cell's fields, and load the build's Julia engine once.
#
#   Rscript prepare.R <grid.txt>
#
# run_grid.sh calls this with the same environment the cells get (BENCH_LIB,
# BENCH_DATA, BENCH_DIR). Simulating here, serially, means no two cells ever
# write the same data file, and loading the engine here means a build with a
# new engine precompiles once rather than in sixteen processes at the same time.
args <- commandArgs(TRUE)
grid <- args[1]
HERE <- Sys.getenv("BENCH_DIR", dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))))
STORE <- Sys.getenv("BENCH_DATA", file.path(Sys.getenv("HOME"), "dev", "ctsem-bench-data"))
Sys.setenv(NOT_CRAN = "true", CTSEM_JULIA_AGREE = "yes", JULIA_NUM_THREADS = "1")
lib <- Sys.getenv("BENCH_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(ctsem))
cat("LOADED ctsem", as.character(utils::packageVersion("ctsem")), "from",
  find.package("ctsem"), "\n")
source(file.path(HERE, "cells.R"))
source(file.path(HERE, "starts.R"))
source(file.path(HERE, "variants.R"))

lines <- readLines(grid, warn = FALSE)
lines <- trimws(lines[!grepl("^\\s*(#|$)", lines)])
cells <- do.call(rbind, lapply(strsplit(lines, "\\s+"), function(f) {
  if (length(f) != 6L) stop("a grid line has ", length(f), " fields, not 6: ",
    paste(f, collapse = " "))
  data.frame(id = f[1], model = f[2], data = f[3], route = f[4], start = f[5],
    variant = f[6], stringsAsFactors = FALSE)
}))
if (anyDuplicated(cells$id)) stop("duplicate cell ids: ",
  paste(unique(cells$id[duplicated(cells$id)]), collapse = ", "))
for (i in seq_len(nrow(cells))) {
  P <- bench_problem(cells$model[i])
  if (!cells$route[i] %in% P$routes) stop(cells$id[i], ": route ", cells$route[i],
    " not allowed for ", cells$model[i])
  invisible(bench_variant(cells$variant[i]))
  if (!grepl("^(default:[0-9]+|zeros|stored:.+)$", cells$start[i]))
    stop(cells$id[i], ": bad start ", cells$start[i])
  if (grepl("^stored:", cells$start[i]) && is.null(BENCH_STARTS[[sub("^stored:", "",
    cells$start[i])]])) stop(cells$id[i], ": no stored start ", cells$start[i])
}
cat(nrow(cells), "cells checked\n")
keys <- unique(cells[, c("model", "data")])
for (i in seq_len(nrow(keys))) {
  t0 <- proc.time()[["elapsed"]]
  D <- bench_data(keys$model[i], keys$data[i], STORE)
  cat(sprintf("data %-10s %-4s %6d rows  md5 %s  (%.1fs)\n", keys$model[i],
    keys$data[i], nrow(D$data), D$md5, proc.time()[["elapsed"]] - t0))
}

# The Julia the cells start (harness.R): the build's own pin. An image is a
# Julia version's, so one built here under another is no use to them: this
# looked for 1.12 after the pin moved to 1.13, and the first cells of every
# grid then built the engine again, side by side, 12-15 minutes on dev1.
jbin <- Sys.getenv("BENCH_JULIA_BIN")
if (!nzchar(jbin)) {
  pin <- paste0("julia-", get0(".ct_julia_version", envir = asNamespace("ctsem"),
    ifnotfound = "1.12"), "*")
  cand <- c(Sys.glob(file.path(Sys.getenv("HOME"), ".julia", "juliaup", pin, "bin")),
    Sys.glob(file.path(Sys.getenv("USERPROFILE"), ".julia", "juliaup", pin, "bin")))
  cand <- cand[file.exists(cand)]
  jbin <- if (length(cand)) sort(cand, decreasing = TRUE)[1L] else ""
}
t0 <- proc.time()[["elapsed"]]
if (nzchar(jbin)) suppressMessages(ctJuliaSetup(julia_bin = jbin, threads = 1L)) else
  suppressMessages(ctJuliaSetup(threads = 1L))
st <- ctJuliaStatus()
cat(sprintf("engine ready in %.1fs: Julia %s, %s\n", proc.time()[["elapsed"]] - t0,
  as.character(JuliaConnectoR::juliaEval("string(VERSION)")),
  paste(format(st$engine %||% st), collapse = " ")))
# One small fit, so every path a cell will take has been through this build's
# engine at least once and a broken build fails here rather than in 200 cells.
D <- bench_data("small", "1", STORE)
t0 <- proc.time()[["elapsed"]]
f <- suppressWarnings(suppressMessages(ctFit(D$data[D$data$id <= 10, ],
  bench_problem("small")$model(), backend = "julia", cores = 1L,
  optimcontrol = list(maxiter = 20L))))
cat(sprintf("smoke fit %.1fs: loglik %.4f, %d iterations\n",
  proc.time()[["elapsed"]] - t0, as.numeric(f$estimate$loglik)[1],
  as.integer(f$optim$iterations)))
cat("PREPARED\n")
