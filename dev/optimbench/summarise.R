# Summarise bench results: a markdown report and a CSV.
#
#   Rscript summarise.R [--out <prefix>] [--refs <references.csv>] <results dir> [<results dir> ...]
#
# Each results dir is one grid run against one build (run_grid.sh's
# <results>/<label>/<grid>). With several, the first is the baseline and every
# cell id found in more than one is compared pairwise against it.
#
# Scoring. A cell's score is, in order of preference, its penalised exact log
# likelihood (Laplace cells whose units allow the reference; see reference.jl),
# else its objective re-scored on a freshly built objective, else the fit's
# log posterior. Cells are compared only with cells of the same model, data,
# route, objective tag (variants.R) and score kind; `best` is the highest score
# in that group across every results dir given, and for exact scores also the
# best-known value in references.csv. `loss` = score - best, so 0 is the best
# known and a loss below -0.05 means the fit ended materially short of it.
# `dx` is the largest raw-scale distance from the group's best point.
#
# Printed first: the contamination check (a quantity no optimiser change can
# affect) and the check that each fit's warm-up flag says what ran.
args <- commandArgs(TRUE)
HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
opt <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) return(default)
  v <- args[i + 1L]; args <<- args[-c(i, i + 1L)]; v
}
outprefix <- opt("--out", "bench-summary")
reffile <- opt("--refs", file.path(HERE, "references.csv"))
title <- opt("--title", "Optimiser bench")
dirs <- args
if (!length(dirs)) stop("usage: Rscript summarise.R [--out prefix] <results dir> ...")
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
source(file.path(HERE, "starts.R"))

bind <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (k in setdiff(cols, names(r))) r[[k]] <- NA
    r[, cols, drop = FALSE]
  }))
}

raws <- list(); records <- list(); rows <- list(); order_in_grid <- list()
for (dir in dirs) {
  files <- list.files(dir, "\\.rds$", full.names = TRUE)
  grid <- file.path(dir, "grid.txt")
  gl <- if (file.exists(grid)) {
    l <- trimws(readLines(grid, warn = FALSE)); l <- l[!grepl("^(#|$)", l)]
    do.call(rbind, lapply(strsplit(l, "\\s+"), function(f) data.frame(id = f[1],
      model = f[2], data = f[3], route = f[4], start = f[5], variant = f[6],
      stringsAsFactors = FALSE)))
  } else NULL
  status <- if (file.exists(file.path(dir, "status.tsv"))) utils::read.delim(
    file.path(dir, "status.tsv"), header = FALSE, col.names = c("id", "exit",
      "secs", "load1"), stringsAsFactors = FALSE) else NULL
  seen <- character()
  label <- NA_character_
  for (f in files) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r) || is.null(r$cell)) next
    label <- r$label
    key <- paste(r$label, r$cell$id)
    records[[key]] <- r
    row <- r$row
    if (is.null(row)) row <- data.frame(id = r$cell$id, model = r$cell$model,
      data = r$cell$data, route = r$cell$route, start = r$cell$start,
      variant = r$cell$variant, label = r$label, host = r$host,
      status = r$status, error = r$error, stringsAsFactors = FALSE)
    row$dir <- dir
    rows[[key]] <- row
    raws[[key]] <- if (!is.null(r$estimate$raw)) r$estimate$raw else
      if (identical(r$status, "evalonly")) r$start_vector else NULL
    seen <- c(seen, r$cell$id)
  }
  if (!is.null(gl)) {
    lab <- if (is.na(label)) basename(dirname(normalizePath(dir))) else label
    for (i in seq_len(nrow(gl))) if (!gl$id[i] %in% seen) {
      st <- status[status$id == gl$id[i], , drop = FALSE]
      rows[[paste(lab, gl$id[i])]] <- data.frame(gl[i, ], label = lab,
        status = if (nrow(st)) paste0("killed(exit ", st$exit[1], ")") else "pending",
        dir = dir, stringsAsFactors = FALSE)
    }
    order_in_grid[[dir]] <- gl$id
  }
}
if (!length(rows)) stop("no results in ", paste(dirs, collapse = ", "))
d <- bind(unname(rows))
d$key <- paste(d$label, d$id)
num <- function(x) suppressWarnings(as.numeric(x))
for (k in intersect(c("exact", "rescored", "logposterior", "loglik", "secs_fit",
  "iterations", "engine_iterations", "engine_g_calls", "engine_f_calls", "ctrl1",
  "ctrl2", "ctrl3", "cert_gap", "quad5", "quad5_gap", "newton_hessians",
  "newton_subset_hessians", "hessians_computed", "warm_iterations",
  "main_iterations", "resume_iterations", "load1_fit_start", "load1_fit_end",
  "is_check_maxdiff", "unit_lambda_min", "warm_dx", "warm_dll", "runs_overshot",
  "overshot_maxabs_min", "runs_maxabs_lt2", "maxabs_raw_max", "probe_secs",
  "engine_runs", "probe_runs", "lapcorr_secs", "lapcorr_gradients", "cont_delta_se_max",
  "hess_secs_central", "hess_secs_forward", "hess_secs_laplace_x", "se_rel_forward",
  "hess_secs_exact", "se_rel_exact",
  "se_rel_laplace_x", "se_rel_laplace_est", "se_rel_reported", "npar",
  "secs_lapcorrect"), names(d))) d[[k]] <- num(d[[k]])
# Written only by the harness since the quadhess job; absent from older results.
for (k in c("hess_secs_exact", "se_rel_exact")) if (is.null(d[[k]])) d[[k]] <- NA_real_
done <- d$status %in% c("ok", "evalonly")

# ---- scores ----------------------------------------------------------------------
d$score_kind <- ifelse(is.finite(d$exact), d$exact_kind %||% "exact",
  ifelse(is.finite(d$rescored), "objective", ifelse(is.finite(d$logposterior),
    "logposterior", NA)))
d$score <- ifelse(is.finite(d$exact), d$exact, ifelse(is.finite(d$rescored),
  d$rescored, d$logposterior))
d$objective[is.na(d$objective)] <- "default"
d$group <- paste(d$model, d$data, d$route, d$objective, d$score_kind)
refs <- if (file.exists(reffile)) utils::read.csv(reffile, stringsAsFactors = FALSE) else NULL
d$ref_known <- NA_real_
if (!is.null(refs)) {
  m <- match(d$model, refs$model)
  exactish <- d$score_kind %in% c("exact-softcut3.5", "is-t4-1.5") & d$objective == "default"
  d$ref_known[exactish & !is.na(m)] <- num(refs$value[m[exactish & !is.na(m)]])
}
d$best <- NA_real_; d$loss <- NA_real_; d$dx <- NA_real_; d$best_from <- NA_character_
for (g in unique(d$group[done & is.finite(d$score)])) {
  i <- which(d$group == g & done & is.finite(d$score))
  b <- max(d$score[i])
  from <- d$key[i][which.max(d$score[i])]
  known <- suppressWarnings(max(d$ref_known[i], na.rm = TRUE))
  if (is.finite(known) && known > b) { b <- known; from <- "references.csv" }
  d$best[i] <- b; d$loss[i] <- d$score[i] - b; d$best_from[i] <- from
  braw <- if (identical(from, "references.csv")) {
    f <- BENCH_STARTS[[paste0("hist_", d$model[i[1]])]]
    tryCatch(f(length(raws[[d$key[i[1]]]] %||% 0), NULL), error = function(e) NULL)
  } else raws[[from]]
  for (j in i) {
    r <- raws[[d$key[j]]]
    if (!is.null(r) && length(r) == length(braw)) d$dx[j] <- max(abs(r - braw))
  }
}
d$hessians_total <- rowSums(cbind(d$newton_hessians, d$newton_subset_hessians,
  d$hessians_computed), na.rm = TRUE)
d$ctrl_min <- suppressWarnings(apply(cbind(d$ctrl1, d$ctrl2, d$ctrl3), 1, min, na.rm = TRUE))
d$ctrl_min[!is.finite(d$ctrl_min)] <- NA

# ---- report ------------------------------------------------------------------------
md <- character()
add <- function(...) md <<- c(md, paste0(...))
# Significant digits for small quantities (gaps, distances, seconds per call),
# fixed decimals for scores and losses, whole numbers for counts.
fmt <- function(x, digits = 3) ifelse(!is.finite(x), "", formatC(x, digits = digits, format = "g"))
f2 <- function(x, d = 2) ifelse(!is.finite(x), "", sprintf(paste0("%.", d, "f"), x))
f0 <- function(x) ifelse(!is.finite(x), "", sprintf("%.0f", x))
tab <- function(df) {
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  add("| ", paste(names(df), collapse = " | "), " |")
  add("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  for (i in seq_len(nrow(df))) add("| ", paste(unlist(df[i, ]), collapse = " | "), " |")
  add("")
}
hosts <- unique(stats::na.omit(d$host))
add("# ", title, "\n")
add("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), " by `dev/optimbench/summarise.R` from ",
  paste0("`", dirs, "`", collapse = ", "), ".\n")
builds <- unique(d[, c("label", "build")])
add("Builds: ", paste(sprintf("`%s` (%s)", builds$label, ifelse(is.na(builds$build), "?",
  builds$build)), collapse = "; "), ". Machine(s): ", paste(hosts, collapse = ", "),
  ". Every timing below is from the machine named; counts are machine-independent.\n")
st <- table(d$label, d$status)
add("Cells by status:\n")
tab(data.frame(label = rownames(st), unclass(st), check.names = FALSE, row.names = NULL))

add("## Contamination check\n")
add("Seconds per value-and-gradient evaluation at raw zeros, timed inside Julia, ",
  "three batches per cell. No optimiser change can move it, so where it moves the ",
  "machine moved. First the control cells (panel regime, 1000 subjects), which run ",
  "at the start, middle and end of the grid:\n")
cc <- d[d$variant %in% "control", c("label", "id", "host", "ctrl1", "ctrl2", "ctrl3",
  "load1_start"), drop = FALSE]
if (nrow(cc)) { cc[, 4:6] <- lapply(cc[, 4:6], fmt, 4); tab(cc) } else add("(no control cells)\n")
add("Then every cell's own control, per model: the spread of its minimum over the ",
  "cells of that model and build. A max/min ratio near 1 says the cells ran on a ",
  "comparable machine; much above 1.3 says some ran contended.\n")
k <- d[is.finite(d$ctrl_min) & !d$variant %in% "control", ]
if (nrow(k)) {
  # By route as well: the augmented and Laplace objectives of one model are
  # different computations and cost different amounts.
  agg <- do.call(rbind, lapply(split(k, list(k$label, k$model, k$route), drop = TRUE), function(x)
    data.frame(label = x$label[1], model = x$model[1], route = x$route[1], cells = nrow(x),
      min = fmt(min(x$ctrl_min), 3), median = fmt(stats::median(x$ctrl_min), 3),
      max = fmt(max(x$ctrl_min), 3), ratio = f2(max(x$ctrl_min) / min(x$ctrl_min)),
      load_fit_start = if (any(is.finite(x$load1_fit_start))) paste(f2(range(
        x$load1_fit_start, na.rm = TRUE), 1), collapse = " to ") else "",
      stringsAsFactors = FALSE)))
  tab(agg[order(agg$label, agg$model), ])
}

add("## Does the warm-up flag say what ran\n")
add("`claimed` is `fit$optim$carefulfit`; `ran` is whether a prior warm-up stage ",
  "actually went through `.ctJuliaOptimise` (the harness's stage record).\n")
w <- d[d$status %in% "ok", ]
if (nrow(w)) {
  wt <- as.data.frame(table(label = w$label, variant = w$variant,
    claimed = w$warmup_claimed, ran = w$warmup_ran), stringsAsFactors = FALSE)
  tab(wt[wt$Freq > 0, ])
}

if ("warm_dx" %in% names(d)) {
  add("## Did the warm-up fit and the timed fit agree\n")
  add("Every fit cell fits twice in one session, the first time as its warm-up, so ",
    "that the timed fit pays no first-call compilation. The two are the same fit ",
    "and must end at the same point: `dx` is the largest raw difference between ",
    "their estimates, `dll` the timed fit's log posterior minus the warm-up's, and ",
    "`warm_s` / `s` the two fits' seconds.\n")
  wf <- d[d$status %in% "ok" & is.finite(d$warm_dx), ]
  if (nrow(wf)) {
    add(sprintf("%d fit cells; %d end at the same point (dx < 1e-8), %d within 1e-4, %d further apart.\n",
      nrow(wf), sum(wf$warm_dx < 1e-8), sum(wf$warm_dx >= 1e-8 & wf$warm_dx < 1e-4),
      sum(wf$warm_dx >= 1e-4)))
    x <- wf[wf$warm_dx >= 1e-8, ]
    if (nrow(x)) tab(data.frame(label = x$label, id = x$id, dx = fmt(x$warm_dx, 2),
      dll = fmt(x$warm_dll, 2), warm_s = f2(num(x$secs_warm), 1), s = f2(x$secs_fit, 1),
      stringsAsFactors = FALSE))
  }
}

add("## By model, route and variant\n")
add("Medians over the cells of each group (seeds); `best` counts cells within 0.05 of ",
  "the group's best score; `worst` is the largest loss. `it` is iterations summed ",
  "over every optimiser stage the fit ran (engine runs, including stall escapes), ",
  "`g` gradient calls likewise, `H` Hessians formed (Newton finish plus certification), ",
  "`s` wall seconds of the fit alone.\n")
ok <- d[d$status %in% c("ok"), ]
if (nrow(ok)) {
  agg <- do.call(rbind, lapply(split(ok, list(ok$label, ok$model, ok$route, ok$variant),
    drop = TRUE), function(x) data.frame(label = x$label[1], model = x$model[1],
      route = x$route[1], variant = x$variant[1], n = nrow(x),
      conv = sum(x$converged %in% TRUE), cert = sum(x$cert_status %in% "certified"),
      it = f0(stats::median(x$engine_iterations)),
      g = f0(stats::median(x$engine_g_calls)),
      H = f0(stats::median(x$hessians_total)),
      s = f2(stats::median(x$secs_fit), 1),
      best = sum(is.finite(x$loss) & x$loss > -0.05),
      worst = f2(suppressWarnings(min(x$loss, na.rm = TRUE)), 3),
      warm = sum(x$warmup_ran %in% TRUE), stringsAsFactors = FALSE)))
  agg$worst[agg$worst %in% c("Inf", "-Inf")] <- ""
  tab(agg[order(agg$label, agg$model, agg$route, agg$variant), ])
}

if ("runs_overshot" %in% names(d)) {
  add("## The overshoot probe\n")
  add("Every engine optimisation run ends on the overshoot probe (`_ctsem_overshot`), ",
    "at the run's end point. Per cell: `runs` engine runs, `found` how many of their ",
    "probes found a gain, `at` the smallest `max |raw|` of an end point where one did, ",
    "`lt2` runs whose end point had every `|raw| < 2`, `maxraw` the largest `max |raw|` ",
    "over the runs, `probed` the runs whose probe ran (the prior warm-up's is off), ",
    "`probe_s` one probe at the estimate (timed after the fit), and `share` = ",
    "probed x probe_s / fit seconds, the probe's estimated share of the fit.\n")
  pc <- d[d$status %in% "ok" & !d$variant %in% "noprobe", ]
  if (nrow(pc)) {
    add(sprintf(paste0("%d fits, %d engine runs; %d runs found a gain, in %d cells; ",
      "%d runs ended with every |raw| < 2, and %d of those found one.\n"), nrow(pc),
      sum(pc$engine_runs, na.rm = TRUE), sum(pc$runs_overshot, na.rm = TRUE),
      sum(pc$runs_overshot > 0, na.rm = TRUE), sum(pc$runs_maxabs_lt2, na.rm = TRUE),
      sum(pc$runs_overshot > 0 & is.finite(pc$overshot_maxabs_min) & pc$overshot_maxabs_min < 2)))
    pc <- pc[order(pc$model, pc$route, pc$variant, pc$id), ]
    tab(data.frame(label = pc$label, id = pc$id, runs = f0(pc$engine_runs),
      found = f0(pc$runs_overshot), at = f2(pc$overshot_maxabs_min, 2),
      lt2 = f0(pc$runs_maxabs_lt2), maxraw = f2(pc$maxabs_raw_max, 2),
      probed = f0(pc$probe_runs), probe_s = f2(pc$probe_secs, 2), s = f2(pc$secs_fit, 1),
      share = f2(pc$probe_runs * pc$probe_secs / pc$secs_fit, 3), stringsAsFactors = FALSE))
  }
  np <- d[d$status %in% "ok" & d$variant %in% "noprobe", ]
  if (nrow(np)) {
    add("Paired with the probe off (`noprobe`, same cell otherwise): `ratio` is the ",
      "fit seconds off over on, `d_score` the score off minus on, `dx` the largest raw ",
      "difference between the two estimates.\n")
    rows <- lapply(seq_len(nrow(np)), function(i) {
      b <- d[d$status %in% "ok" & d$variant %in% "default" & d$label == np$label[i] &
        d$model == np$model[i] & d$data == np$data[i] & d$route == np$route[i] &
        d$start == np$start[i], ]
      if (!nrow(b)) return(NULL)
      rb <- raws[[b$key[1]]]; rn <- raws[[np$key[i]]]
      data.frame(id = np$id[i], s_on = f2(b$secs_fit[1], 1), s_off = f2(np$secs_fit[i], 1),
        ratio = f2(np$secs_fit[i] / b$secs_fit[1], 3),
        d_score = fmt(np$score[i] - b$score[1], 2),
        dx = if (length(rb) && length(rb) == length(rn)) fmt(max(abs(rb - rn)), 2) else "",
        found_on = f0(b$runs_overshot[1]), stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) tab(do.call(rbind, rows))
  }
}

if ("cont_hessian" %in% names(d)) {
  add("## The quadrature continuation's Hessian\n")
  add("Laplace fits whose quadrature continuation reached its fixed point and reported ",
    "the Hessian of its own objective. Against that Hessian by central differences of ",
    "the hybrid gradient (2 npar gradients), the largest relative difference in any ",
    "standard error from: `exact` the same Hessian by the exact scheme (on a build ",
    "that has it); `fwd` by forward differences (npar + 1 gradients); ",
    "`lap_x` the Laplace Hessian at the continuation's estimate; `lap_est` the Laplace ",
    "Hessian the fit already held, at the Laplace optimum. `rep` checks the rule: the ",
    "fit's reported standard errors against the Hessian it reported. `moved` is the ",
    "continuation's largest move in Laplace standard errors; seconds are dev1's, ",
    "taken after the fit, one after another; `corr_s` is the whole correction in the fit. ",
    "Before the quadhess harness the reference was the Hessian the fit reported, ",
    "whichever scheme its build took it by.\n")
  ch <- d[d$status %in% "ok" & d$cont_hessian %in% TRUE, ]
  cs <- d[d$status %in% "ok" & !is.na(d$lapcorr_status), ]
  if (nrow(cs)) {
    st <- as.data.frame(table(status = cs$lapcorr_status,
      continuation = ifelse(is.na(cs$lapcorr_continuation), "", cs$lapcorr_continuation),
      hessian = cs$cont_hessian), stringsAsFactors = FALSE)
    tab(st[st$Freq > 0, ])
  }
  if (nrow(ch)) tab(data.frame(label = ch$label, id = ch$id, npar = f0(ch$npar),
    moved = f2(ch$cont_delta_se_max, 2), exact = fmt(ch$se_rel_exact, 2),
    fwd = fmt(ch$se_rel_forward, 2),
    lap_x = fmt(ch$se_rel_laplace_x, 2), lap_est = fmt(ch$se_rel_laplace_est, 2),
    rep = fmt(ch$se_rel_reported, 2), exact_s = f2(ch$hess_secs_exact, 1),
    central_s = f2(ch$hess_secs_central, 1),
    fwd_s = f2(ch$hess_secs_forward, 1), lap_x_s = f2(ch$hess_secs_laplace_x, 1),
    corr_s = f2(ch$lapcorr_secs, 1), s = f2(ch$secs_fit, 1), stringsAsFactors = FALSE))
}

add("## Reference checks\n")
add("Each `evalonly` cell evaluates the reference at a config's best-known point ",
  "(starts.R `hist_*`); `diff` is this build's value there minus references.csv. ",
  "`IS` is the largest difference between the reference and importance sampling ",
  "on the three units where the reference departs most from Laplace.\n")
ev <- d[d$status %in% "evalonly", ]
if (nrow(ev)) {
  m <- match(ev$model, refs$model)
  tab(data.frame(label = ev$label, id = ev$id, kind = ev$exact_kind,
    value = f2(ev$exact, 4), known = f2(num(refs$value[m]), 4),
    diff = fmt(ev$exact - num(refs$value[m]), 3), IS = fmt(ev$is_check_maxdiff, 2),
    quad5_gap = f2(ev$quad5_gap, 3), secs = f2(num(ev$secs_ref), 0),
    stringsAsFactors = FALSE))
}

add("## Every cell\n")
add("`it` warm/main/resume iterations; `g` engine gradient calls in total; `H` ",
  "Hessians; `cert` the certification status and gap; `score` as above with its ",
  "kind in `kind`; `loss` against the group best; `dx` max raw distance from the ",
  "group's best point; `quad` the 5-node quadrature gap (quadrature minus Laplace) at ",
  "the end point; `s` fit seconds; `load` 1-minute load at the fit's start.\n")
fams <- c("control", "cf_", "gA", "gB", "gC", "gD", "gN", "acnonlin", "mvmix", "jflat",
  "small", "nonlin", "long", "bigp", "panel", "ordinal", "anomS")
famof <- function(m) { f <- fams[vapply(fams, function(p) startsWith(m, p), TRUE)]
  if (length(f)) f[1] else "other" }
d$family <- vapply(d$model, famof, "")
for (fam in unique(d$family[!d$family %in% "control"])) {
  x <- d[d$family == fam, ]
  x <- x[order(x$model, x$route, x$variant, x$label, x$id), ]
  add("### ", fam, "\n")
  tab(data.frame(label = x$label, id = x$id, status = x$status,
    conv = x$converged, cert = ifelse(is.na(x$cert_status), "",
      paste0(x$cert_status, " ", fmt(x$cert_gap, 2))),
    it = ifelse(is.na(x$engine_iterations), "", paste0(x$warm_iterations %||% "", "/",
      x$main_iterations, "/", x$resume_iterations)),
    g = f0(x$engine_g_calls), H = f0(x$hessians_total),
    warm = ifelse(is.na(x$warmup_ran), "", paste0(substr(x$warmup_claimed, 1, 1), "/",
      substr(x$warmup_ran, 1, 1))),
    score = f2(x$score, 3), kind = x$score_kind, loss = f2(x$loss, 3),
    dx = fmt(x$dx, 2), quad = f2(x$quad5_gap, 2), s = f2(x$secs_fit, 1),
    load = f2(x$load1_fit_start, 1), stringsAsFactors = FALSE))
}

if (length(unique(d$label)) > 1L) {
  add("## Paired against the first results dir\n")
  base <- d[d$dir == dirs[1] & d$status %in% c("ok", "evalonly"), ]
  other <- d[d$dir != dirs[1] & d$status %in% c("ok", "evalonly"), ]
  p <- merge(base, other, by = "id", suffixes = c(".base", ".new"))
  if (nrow(p)) tab(data.frame(id = p$id, label = p$label.new,
    d_loss = f2(p$loss.new - p$loss.base, 3),
    d_it = p$engine_iterations.new - p$engine_iterations.base,
    d_g = p$engine_g_calls.new - p$engine_g_calls.base,
    d_H = p$hessians_total.new - p$hessians_total.base,
    ratio_s = f2(p$secs_fit.new / p$secs_fit.base, 2),
    cert = paste(p$cert_status.base, "->", p$cert_status.new),
    stringsAsFactors = FALSE))
}

errs <- d[d$status %in% c("error", "timeout") | grepl("^killed", d$status), ]
if (nrow(errs)) {
  add("## Errors, timeouts and killed cells\n")
  tab(data.frame(label = errs$label, id = errs$id, status = errs$status,
    error = substr(gsub("[|\n]", " ", errs$error %||% ""), 1, 200), stringsAsFactors = FALSE))
}

writeLines(md, paste0(outprefix, ".md"))
keep <- setdiff(names(d), c("dir", "key", "family"))
utils::write.csv(d[, keep], paste0(outprefix, ".csv"), row.names = FALSE)
cat("wrote", paste0(outprefix, ".md"), "and", paste0(outprefix, ".csv"), ":", nrow(d), "cells,",
  sum(done), "finished\n")
