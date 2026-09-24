# Summarise a sweep: Rscript summarise.R <sweepdir>
args <- commandArgs(TRUE)
dir <- args[1]
files <- list.files(dir, "\\.rds$", full.names = TRUE)
rows <- list()
for (f in files) {
  r <- readRDS(f)
  b <- r$baseline
  basetot <- b$secs + b$hess_secs
  common <- data.frame(model = r$model, seed = r$seed, npar = r$npar,
    N = r$nsubjects, base_it = b$iterations, base_secs = b$secs,
    hess_secs = b$hess_secs, base_total = basetot,
    noP_it = b$noP$iterations, noP_total = b$noP$secs + b$hess_secs,
    noP_df = b$noP$f - b$f, stringsAsFactors = FALSE)
  if (!is.null(r$pbatch)) for (p in r$pbatch) {
    rows[[length(rows) + 1L]] <- cbind(common, method = sprintf("PB metric=%s grow=%s th=%g",
      p$metric, p$grow, p$theta), total = p$pb_secs + p$eg_secs,
      df = p$f - b$f, status = paste(p$pb_status, p$eg_status),
      detail = sprintf("it=%d+%d sizes=%s", p$pb_iterations, p$eg_iterations,
        paste(p$sizes, collapse = ",")), stringsAsFactors = FALSE)
  }
  if (!is.null(r$endgame)) for (e in r$endgame) {
    rows[[length(rows) + 1L]] <- cbind(common, method = sprintf("EG tau=%g %s",
      e$tau, e$curvature), total = e$stage1_secs + e$endgame_secs,
      df = e$f - b$f, status = e$status,
      detail = sprintf("it=%d+%d refresh=%d H=%.2f", e$stage1_iterations,
        e$iterations, e$refreshes, e$hess), stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
d$speedup <- d$base_total / d$total
d$speedup_vs_noP <- d$noP_total / d$total
saveRDS(d, file.path(dir, "summary.rds"))
options(width = 250)
cat("== baselines (per model, seed)\n")
bl <- unique(d[, c("model", "seed", "npar", "N", "base_it", "base_total", "hess_secs",
  "noP_it", "noP_total", "noP_df")])
print(bl[order(bl$model, bl$seed), ], row.names = FALSE, digits = 3)
cat("\n== median speedup over seeds (vs engine default incl. certification Hessian); worst df\n")
agg <- do.call(rbind, lapply(split(d, list(d$model, d$method), drop = TRUE), function(x)
  data.frame(model = x$model[1], method = x$method[1],
    speedup = median(x$speedup), speedup_noP = median(x$speedup_vs_noP),
    min_speedup = min(x$speedup), worst_df = min(x$df),
    fails = sum(x$df < -1e-3), stringsAsFactors = FALSE)))
agg <- agg[order(agg$model, -agg$speedup), ]
print(agg, row.names = FALSE, digits = 3)
cat("\n== runs landing > 1e-3 below the baseline optimum\n")
bad <- d[d$df < -1e-3, c("model", "seed", "method", "df", "status", "detail")]
print(bad, row.names = FALSE, digits = 4)
