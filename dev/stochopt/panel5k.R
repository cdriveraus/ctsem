# Batching and the Newton finish, separately, on a 5000-subject panel.
#   Rscript panel5k.R <lib> <seed> <outdir>
#
# One process per dataset. It first fits every configuration, with and without
# uncertainty, on a 500-subject slice -- enough subjects to batch -- so that
# every code path the timed fits use is compiled (compiled code does not depend
# on the number of subjects). Then it times all eight on the full data, the
# order rotated by seed so no configuration always runs first.
#
# Configurations, all on the new optimiser: none (own L-BFGS alone), batch,
# newton, both (the default).
args <- commandArgs(TRUE)
LIB <- args[1]; SEED <- as.integer(args[2]); OUT <- args[3]
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v3")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem("panel5k", SEED)
configs <- c("none", "batch", "newton", "both")
control <- function(config, estonly) {
  oc <- list(estonly = estonly)
  if (config %in% c("none", "newton")) oc$batch <- FALSE
  if (config %in% c("none", "batch")) oc$newton <- FALSE
  oc
}
fitit <- function(d, config, estonly) {
  set.seed(SEED)
  suppressWarnings(suppressMessages(ctFit(d, P$m, backend = "julia",
    intoverpop = P$route, cores = 1, optimcontrol = control(config, estonly))))
}
ids <- unique(P$d$id)
slice <- P$d[P$d$id %in% ids[1:500], ]
for (config in configs) for (estonly in c(TRUE, FALSE)) {
  w <- fitit(slice, config, estonly)
  cat("WARM", config, estonly, "batch", paste(w$optim$batch_sizes, collapse = ","),
    "newton", w$optim$newton_steps, "\n"); flush(stdout())
}
g <- function(x) if (is.null(x)) NA else x
jobs <- expand.grid(config = configs, estonly = c(TRUE, FALSE),
  stringsAsFactors = FALSE)
jobs <- jobs[((seq_len(nrow(jobs)) + 3 * (SEED - 1) - 1) %% nrow(jobs)) + 1, ]
rows <- list()
for (j in seq_len(nrow(jobs))) {
  config <- jobs$config[j]; estonly <- jobs$estonly[j]
  t0 <- proc.time()[["elapsed"]]
  f <- fitit(P$d, config, estonly)
  secs <- proc.time()[["elapsed"]] - t0
  o <- f$optim
  row <- data.frame(config = config, estonly = estonly, seed = SEED,
    position = j, secs = secs, loglik = as.numeric(f$estimate$loglik)[1],
    converged = isTRUE(o$converged), iterations = as.integer(g(o$iterations)),
    f_calls = as.integer(g(o$f_calls)), g_calls = as.integer(g(o$g_calls)),
    newton = as.integer(g(o$newton_steps)),
    batch = paste(g(o$batch_sizes), collapse = ","),
    batch_at = paste(g(o$batch_iterations), collapse = ","),
    status = as.character(g(f$uncertainty$certification$status)),
    stringsAsFactors = FALSE)
  print(row); flush(stdout())
  rows[[j]] <- row
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(do.call(rbind, rows), file.path(OUT, sprintf("seed-%d.rds", SEED)))
cat("SAVED\n")
