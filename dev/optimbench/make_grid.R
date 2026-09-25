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

# The baseline covers every family at the build's defaults, with few seeds:
# 2 default starts per cell, 3 for the gated-gaps configs, and of those the
# gaps note's worst and one typical config per family (A14/A1, B8/B2, C2/C8,
# D3/D1, N1/N3). Cut from the first design to keep dev1 time down: the second
# cf data set, the third cf seed, seeds 4-5 and configs A10, B1, C1, N2, N4,
# the third fixture and regime seeds, the second AnomAuth default seed, and
# the reference checks outside one config per family. Fast families first, so
# a first summary exists while the slow ones run; the three control cells sit
# at the start, the middle and the end of the launch order.
gated <- c("gA1", "gA14", "gB2", "gB8", "gC2", "gC8")
binary <- c("gD1", "gD3")
nested <- c("gN1", "gN3")
anom <- c("anomS1", "anomS2")
hist <- c("gA1", "gB2", "gC2", "gD1", "gN1", "anomS1")

write_grid("baseline", list(
  control(1),
  # carefulfit / gaptol design: both routes, with and without the priors that
  # switch the warm-up off at this build.
  cells(paste0("cf_", c("gaussian", "binary", "ordinal", "mixed")), "1",
    c("augmented", "laplace"), seeds(2), c("default", "priorsFALSE")),
  cells(gated, "cfg", "laplace", seeds(3), "default"),
  cells(binary, "cfg", "laplace", seeds(3), "default"),
  # test fixtures
  cells("acnonlin", "cfg", "laplace", seeds(2), "default"),
  cells("mvmix", "cfg", c("laplace", "augmented"), seeds(2), "default"),
  control(2),
  cells("jflat", "cfg", "laplace", c("stored:flatdrift8", seeds(1)), "default"),
  cells(nested, "cfg", "laplace", seeds(3), "default"),
  # stochopt regimes, each on its own route
  cells(c("small", "nonlin", "long", "bigp", "panel"), "1", "augmented", seeds(2), "default"),
  cells("ordinal", "1", "laplace", seeds(2), "default"),
  # AnomAuth S1 and S2 from the default path and from the stored spurious maxima
  cells(anom, "cfg", "laplace", seeds(1), "default"),
  cell("anomS1", "cfg", "laplace", "stored:anomS1_spurious", "default"),
  cell("anomS2", "cfg", "laplace", "stored:anomS2_spurious", "default"),
  # the reference at a config's best-known point, which must reproduce
  # references.csv (a check of data, model, prior term and reference code)
  do.call(rbind, lapply(hist, function(m) cell(m, "cfg", "laplace",
    paste0("stored:hist_", m), "evalonly"))),
  control(3)),
  "baseline grid: every family at the build's defaults, few seeds (make_grid.R)")

# A quick grid to check a build and the harness end to end in a few minutes.
write_grid("smoke", list(control(1),
  cells("cf_gaussian", "1", c("augmented", "laplace"), seeds(1), "default"),
  cells("small", "1", "augmented", seeds(1), "default"),
  cell("gA1", "cfg", "laplace", "default:1", "estonly")),
  "smoke grid: minutes, not a measurement (make_grid.R)")
