# Write the bench's grids as text, one cell per line, so that what a grid holds
# is reviewable and a job can copy one and edit it.
#
#   Rscript make_grid.R            writes grids/*.txt next to this script
#
# Line format (whitespace-separated, # comments):
#   id  model  data  route  start  variant
# `data` is the simulator's seed, or `cfg` for a problem that fixes its own
# data. `start` is default:<seed> (the user's path), zeros or stored:<name>.
HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
dir.create(file.path(HERE, "grids"), showWarnings = FALSE)

abbr <- c(augmented = "aug", laplace = "lap", auto = "auto")
cell <- function(model, data, route, start, variant) {
  s <- if (grepl("^default:", start)) paste0("s", sub("^default:", "", start)) else
    if (identical(start, "zeros")) "z" else sub("^stored:", "", start)
  data.frame(id = sprintf("%s.%s.%s.%s.%s", model, abbr[[route]],
    if (identical(data, "cfg")) "cfg" else paste0("d", data), s, variant),
    model = model, data = data, route = route, start = start, variant = variant,
    stringsAsFactors = FALSE)
}
# Every combination, the start varying fastest.
cells <- function(models, data, routes, starts, variants) {
  g <- expand.grid(start = starts, variant = variants, route = routes,
    data = data, model = models, stringsAsFactors = FALSE)
  out <- do.call(rbind, Map(cell, g$model, g$data, g$route, g$start, g$variant))
  rownames(out) <- NULL
  out
}
seeds <- function(n) paste0("default:", seq_len(n))
control <- function(k) data.frame(id = paste0("control.", k), model = "panel",
  data = "1", route = "augmented", start = "zeros", variant = "control",
  stringsAsFactors = FALSE)
write_grid <- function(name, parts, header) {
  g <- do.call(rbind, parts)
  stopifnot(!anyDuplicated(g$id))
  f <- file.path(HERE, "grids", paste0(name, ".txt"))
  w <- max(nchar(g$id))
  body <- sprintf(paste0("%-", w, "s  %-12s %-4s %-10s %-24s %s"), g$id, g$model,
    g$data, g$route, g$start, g$variant)
  writeLines(c(paste0("# ", header), "#",
    sprintf(paste0("# %-", w - 2, "s  %-12s %-4s %-10s %-24s %s"), "id", "model",
      "data", "route", "start", "variant"), body), f)
  cat(name, ":", nrow(g), "cells ->", f, "\n")
}

gated <- c("gA1", "gA10", "gA14", "gB1", "gB2", "gB8", "gC1", "gC2", "gC8")
binary <- c("gD1", "gD3")
nested <- paste0("gN", 1:4)
anom <- c("anomS1", "anomS2")
hist <- c(gated, binary, nested, anom)

# The baseline: every cell family once, at the build's defaults. Fast families
# first, so a first summary exists while the slow ones run; the three control
# cells sit at the start, the middle and the end of the launch order.
write_grid("baseline", list(
  control(1),
  # carefulfit / gaptol design: two data sets, both routes, 3 default starts,
  # with and without the priors that switch the warm-up off at this build.
  cells(paste0("cf_", c("gaussian", "binary", "ordinal", "mixed")), c("1", "2"),
    c("augmented", "laplace"), seeds(3), c("default", "priorsFALSE")),
  # gated-gaps configs: 5 default starts each (the note's default path).
  cells(gated, "cfg", "laplace", seeds(5), "default"),
  cells(binary, "cfg", "laplace", seeds(5), "default"),
  # test fixtures
  cells("acnonlin", "cfg", "laplace", seeds(3), "default"),
  cells("mvmix", "cfg", c("laplace", "augmented"), seeds(3), "default"),
  control(2),
  cells("jflat", "cfg", "laplace", c("stored:flatdrift8", seeds(3)), "default"),
  cells(nested, "cfg", "laplace", seeds(5), "default"),
  # stochopt regimes, each on its own route
  cells(c("small", "nonlin", "long", "bigp", "panel"), "1", "augmented", seeds(3), "default"),
  cells("ordinal", "1", "laplace", seeds(3), "default"),
  # AnomAuth S1 and S2 from the default path and from the stored spurious maxima
  cells(anom, "cfg", "laplace", seeds(3), "default"),
  cell("anomS1", "cfg", "laplace", "stored:anomS1_spurious", "default"),
  cell("anomS2", "cfg", "laplace", "stored:anomS2_spurious", "default"),
  # the reference at each config's best-known point, which must reproduce
  # references.csv (a check of data, model, prior term and reference code)
  do.call(rbind, lapply(hist, function(m) cell(m, "cfg", "laplace",
    paste0("stored:hist_", m), "evalonly"))),
  control(3)),
  "baseline grid: every family at the build's defaults (make_grid.R)")

# A quick grid to check a build and the harness end to end in a few minutes.
write_grid("smoke", list(control(1),
  cells("cf_gaussian", "1", c("augmented", "laplace"), seeds(1), "default"),
  cells("small", "1", "augmented", seeds(1), "default"),
  cell("gA1", "cfg", "laplace", "default:1", "estonly")),
  "smoke grid: minutes, not a measurement (make_grid.R)")
