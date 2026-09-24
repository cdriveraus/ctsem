# Batching and the Newton finish, separately, on a 5000-subject panel.
#   Rscript panel5k.R <lib> <config> <seed> <outdir>
# config: old (Optim build) | none | batch | newton | both
args <- commandArgs(TRUE)
LIB <- args[1]; CONFIG <- args[2]; SEED <- as.integer(args[3]); OUT <- args[4]
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v3")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem("panel5k", SEED)
oc <- list(estonly = FALSE)
if (CONFIG %in% c("none", "newton")) oc$batch <- FALSE
if (CONFIG %in% c("none", "batch")) oc$newton <- FALSE
fitit <- function(d) {
  set.seed(SEED)
  suppressWarnings(suppressMessages(ctFit(d, P$m, backend = "julia",
    intoverpop = P$route, cores = 1, optimcontrol = oc)))
}
# Warm-up large enough to batch, so compilation is not in the timing.
ids <- unique(P$d$id)
invisible(try(fitit(P$d[P$d$id %in% ids[1:200], ]), silent = TRUE))
cat("WARM\n"); flush(stdout())
t0 <- proc.time()[["elapsed"]]
f <- fitit(P$d)
secs <- proc.time()[["elapsed"]] - t0
o <- f$optim
g <- function(x) if (is.null(x)) NA else x
row <- data.frame(config = CONFIG, seed = SEED, secs = secs,
  loglik = as.numeric(f$estimate$loglik)[1], converged = isTRUE(o$converged),
  iterations = as.integer(g(o$iterations)), f_calls = as.integer(g(o$f_calls)),
  g_calls = as.integer(g(o$g_calls)), newton = as.integer(g(o$newton_steps)),
  batch = paste(g(o$batch_sizes), collapse = ","),
  status = as.character(g(f$uncertainty$certification$status)),
  gap = as.numeric(g(f$uncertainty$certification$gap)),
  stringsAsFactors = FALSE)
print(row)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(row = row, raw = f$estimate$raw), file.path(OUT,
  sprintf("%s-%d.rds", CONFIG, SEED)))
cat("SAVED\n")
