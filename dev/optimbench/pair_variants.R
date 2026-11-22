# Pair variants of the same cell within one summary CSV (summarise.R's output):
# for each variant against `default`, the score change, gradient and time
# ratios, and which cells moved.
#   Rscript pair_variants.R <summary.csv> [out.csv]
a <- commandArgs(trailingOnly = TRUE)
d <- read.csv(a[1], stringsAsFactors = FALSE)
d <- d[d$status %in% "ok", ]
d$base <- sub("\\.[^.]+$", "", d$id)
grads <- function(x) ifelse(is.na(x$engine_g_calls), x$g_calls, x$engine_g_calls)
ref <- d[d$variant %in% "default", ]
rows <- list()
for (v in setdiff(unique(d$variant), "default")) {
  x <- d[d$variant %in% v, ]
  m <- merge(ref, x, by = "base", suffixes = c(".0", ".1"))
  if (!nrow(m)) next
  ds <- m$score.1 - m$score.0
  g0 <- ifelse(is.na(m$engine_g_calls.0), m$g_calls.0, m$engine_g_calls.0)
  g1 <- ifelse(is.na(m$engine_g_calls.1), m$g_calls.1, m$engine_g_calls.1)
  rows[[v]] <- data.frame(variant = v, pairs = nrow(m),
    better = sum(ds > 0.01, na.rm = TRUE), worse = sum(ds < -0.01, na.rm = TRUE),
    net = round(sum(ds, na.rm = TRUE), 2),
    worst = round(min(ds, na.rm = TRUE), 3), best = round(max(ds, na.rm = TRUE), 3),
    grad_ratio_sum = round(sum(g1, na.rm = TRUE) / sum(g0, na.rm = TRUE), 3),
    grad_ratio_med = round(stats::median(g1 / g0, na.rm = TRUE), 3),
    time_ratio_sum = round(sum(m$secs_fit.1, na.rm = TRUE) /
      sum(m$secs_fit.0, na.rm = TRUE), 3),
    certified = paste0(sum(m$cert_status.1 %in% "certified"), "/",
      sum(m$cert_status.0 %in% "certified")))
  moved <- m[abs(ds) > 0.01 & !is.na(ds), ]
  if (nrow(moved)) {
    cat("\n", v, ": cells that moved\n", sep = "")
    print(data.frame(cell = moved$base, d_score = round(moved$score.1 - moved$score.0, 3),
      grads0 = g0[abs(ds) > 0.01 & !is.na(ds)], grads1 = g1[abs(ds) > 0.01 & !is.na(ds)],
      cert0 = moved$cert_status.0, cert1 = moved$cert_status.1), row.names = FALSE)
  }
}
out <- do.call(rbind, rows)
cat("\n"); print(out, row.names = FALSE)
if (length(a) > 1) write.csv(out, a[2], row.names = FALSE)
