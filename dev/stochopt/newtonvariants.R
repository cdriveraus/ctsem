# The Newton finish's curvature: exact (refresh on slow contraction), chord
# (keep the hand-over Hessian), subset (a scaled subset Hessian for the steps).
#   Rscript newtonvariants.R <lib> <model> <seed> <outdir>
# Warm first on a slice with every variant, then time each on the full data,
# order rotated by seed. Everything else is the default fit.
args <- commandArgs(TRUE)
LIB <- args[1]; MODEL <- args[2]; SEED <- as.integer(args[3]); OUT <- args[4]
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v3")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem(MODEL, SEED)
variants <- c("exact", "chord", "subset")
fitit <- function(d, v) {
  set.seed(SEED)
  suppressWarnings(suppressMessages(ctFit(d, P$m, backend = "julia",
    intoverpop = P$route, cores = 1,
    optimcontrol = list(estonly = FALSE, newton = v))))
}
ids <- unique(P$d$id)
slice <- P$d[P$d$id %in% ids[seq_len(min(length(ids), 500))], ]
for (v in variants) invisible(fitit(slice, v))
cat("WARM\n"); flush(stdout())
g <- function(x) if (is.null(x)) NA else x
order <- variants[((seq_along(variants) + SEED - 2) %% length(variants)) + 1]
rows <- list()
for (j in seq_along(order)) {
  v <- order[j]
  t0 <- proc.time()[["elapsed"]]
  f <- fitit(P$d, v)
  secs <- proc.time()[["elapsed"]] - t0
  o <- f$optim
  row <- data.frame(model = MODEL, variant = v, seed = SEED, position = j,
    secs = secs, loglik = as.numeric(f$estimate$loglik)[1],
    converged = isTRUE(o$converged), iterations = as.integer(g(o$iterations)),
    newton = as.integer(g(o$newton_steps)),
    full_hessians = as.integer(g(o$newton_hessians)),
    subset_hessians = as.integer(g(o$newton_subset_hessians)),
    status = as.character(g(f$uncertainty$certification$status)),
    gap = as.numeric(g(f$uncertainty$certification$gap)),
    se_max = max(abs(g(f$estimate$se))), stringsAsFactors = FALSE)
  print(row); flush(stdout())
  rows[[j]] <- row
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(do.call(rbind, rows), file.path(OUT, sprintf("%s-%d.rds", MODEL, SEED)))
cat("SAVED\n")
