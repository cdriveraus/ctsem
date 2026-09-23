# End to end: ctFit(backend = 'julia') with its defaults, on one build.
#   Rscript e2e.R <lib> <label> <model> <seed> <outdir>
args <- commandArgs(TRUE)
LIB <- args[1]; LABEL <- args[2]; MODEL <- args[3]; SEED <- as.integer(args[4])
OUT <- args[5]
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v3")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem(MODEL, SEED)
fitit <- function(d) {
  set.seed(SEED)
  suppressWarnings(suppressMessages(ctFit(d, P$m, backend = "julia",
    intoverpop = P$route, cores = 1, optimcontrol = list(estonly = FALSE))))
}
# warm-up on a small slice so compilation is not in the timing
ids <- unique(P$d$id)
invisible(try(fitit(P$d[P$d$id %in% ids[seq_len(min(length(ids), 10))], ]), silent = TRUE))
cat("WARM\n"); flush(stdout())
t0 <- proc.time()[["elapsed"]]
f <- try(fitit(P$d), silent = TRUE)
secs <- proc.time()[["elapsed"]] - t0
if (inherits(f, "try-error")) {
  cat("FAILED", as.character(f), "\n"); quit(status = 1)
}
cert <- f$uncertainty$certification
row <- data.frame(label = LABEL, model = MODEL, seed = SEED, secs = secs,
  loglik = as.numeric(f$estimate$loglik)[1],
  converged = isTRUE(f$optim$converged),
  iterations = as.integer(.subset2(f$optim, "iterations")),
  newton = if (is.null(f$optim$newton_steps)) NA else as.integer(f$optim$newton_steps),
  batch = paste(f$optim$batch_sizes, collapse = ","),
  stopped_by_gap = isTRUE(f$optim$stopped_by_gap),
  predicted_gain = if (is.null(f$optim$predicted_gain)) NA else as.numeric(f$optim$predicted_gain)[1],
  gradient_norm = if (is.null(f$optim$gradient_norm)) NA else as.numeric(f$optim$gradient_norm)[1],
  corrections = length(f$optim$corrections),
  status = if (is.null(cert$status)) NA else as.character(cert$status),
  gap = if (is.null(cert$gap)) NA else as.numeric(cert$gap),
  stringsAsFactors = FALSE)
print(row)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(row = row, raw = f$estimate$raw, se = f$estimate$se),
  file.path(OUT, sprintf("%s-%s-%d.rds", LABEL, MODEL, SEED)))
cat("SAVED\n")
