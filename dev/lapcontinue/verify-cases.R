# The quadrature continuation on the known-broken Laplace cases.
#
#   Rscript dev/lapcontinue/verify-cases.R <tree> <outdir> [case ...]
#
# For each case one uncorrected Laplace fit (laplace_correct = FALSE, cores 1)
# from the case's start; the step correction and the continuation are then
# applied to that same fit post hoc, so all three share the Laplace optimum,
# its Hessian and its draws. At each estimate: the reported log likelihood,
# the 5-node quadrature, and the penalised exact reference of the optimiser
# bench (dev/optimbench/reference.jl: per unit a trapezoid over the soft
# directions for units of dimension <= 3, importance sampling for the nested
# ones), the same quantity as the bench's references.csv best-known values.
#
# Cases (review/OPTIM-consolidation-plan-2026-09-25.md section 10, and the
# coordinator's list of 2026-09-26): AnomAuth S1 and S2 from their stored
# spurious maxima; the gated-gaps shortfall configs B8 (seeds 1-3), A14 and
# D3; N1; the carefulfit study's binary data sets 16 and 27 (the study's own
# simulator, ctGenerate(backend = 'r'), copied from
# dev/simstudies/simstudy-carefulfit.R); the 40-subject nonlinear fixture.
# Data for the bench cells come from the bench's store (BENCH_DATA, default
# ~/dev/ctsem-bench-data), so they are the rows the bench saw.
#
# Writes <outdir>/<case>.rds and appends one line per method to
# <outdir>/summary.tsv. VERIFY_REUSE names an earlier run's <outdir>: a point
# whose estimate is identical to one it scored reuses that reference, so a
# rerun after a change to the continuation pays only for its new points.
# VERIFY_FITS names a directory of the uncorrected Laplace fits: a fit found
# there is loaded rather than refitted, and one fitted is saved there. Neither
# the fit nor the step correction runs any continuation code, so a run testing
# a change to the continuation can reuse both. VERIFY_SKIP_CONTINUE=true leaves
# the continuation out (a run that only makes the fits and their references),
# VERIFY_SKIP_PROBE=true the gB8 line probe. Correctness, not timing: seconds are recorded but are
# only comparable within one run on one machine.
args <- commandArgs(TRUE)
TREE <- normalizePath(args[1], winslash = "/")
OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("CTSEM_JULIA_AGREE"))) Sys.setenv(CTSEM_JULIA_AGREE = "yes")
suppressMessages(devtools::load_all(TREE, compile = FALSE, quiet = TRUE))
stamp <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
stamp("LOADED", TREE)
STORE <- Sys.getenv("BENCH_DATA", path.expand("~/dev/ctsem-bench-data"))
source(file.path(TREE, "dev/optimbench/cells.R"))
source(file.path(TREE, "dev/optimbench/starts.R"))
invisible(.ctJuliaModule())
JuliaConnectoR::juliaEval(sprintf('Core.eval(ContinuousTimeSEM, :(include(%s))); nothing',
  deparse(file.path(TREE, "dev/optimbench/reference.jl"))))
jc <- function(fn, ...) JuliaConnectoR::juliaCall(paste0("ContinuousTimeSEM.", fn), ...)
jg <- function(x) if (inherits(x, "JuliaProxy")) JuliaConnectoR::juliaGet(x) else x
jv <- function(x) .ctJuliaNumericVector(as.numeric(x))
now <- function() proc.time()[["elapsed"]]
known <- utils::read.csv(file.path(TREE, "dev/optimbench/references.csv"),
  stringsAsFactors = FALSE)

# The carefulfit study's binary design (dev/simstudies/simstudy-carefulfit.R).
cfs_data <- function(seed) {
  set.seed(seed)
  cints <- stats::rnorm(60, 0, 0.5)
  d <- do.call(rbind, lapply(seq_len(60), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
      manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
      DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(1e-6),
      T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(cints[i]),
      MANIFESTMEANS = matrix(0), Tpoints = 10))
    one <- data.frame(ctGenerate(gen, n = 1, Tpoints = 10, backend = "r"))
    one$id <- i
    one
  }))
  eta <- d$eta
  for (nm in c("b1", "b2", "b3")) d[[nm]] <- stats::rbinom(length(eta), 1,
    1 / (1 + exp(-eta)))
  d$eta <- NULL
  d
}
cfs_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 3, manifestNames = c("b1", "b2", "b3"), latentNames = "eta1",
    LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 3),
    manifesttype = c(1L, 1L, 1L))))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  m
}

CASES <- list(
  anomS1 = list(model = "anomS1", start = "stored:anomS1_spurious", best = "hist_anomS1"),
  anomS2 = list(model = "anomS2", start = "stored:anomS2_spurious", best = "hist_anomS2"),
  gB8s1 = list(model = "gB8", start = "default:1", best = "hist_gB8"),
  gB8s2 = list(model = "gB8", start = "default:2", best = "hist_gB8"),
  gB8s3 = list(model = "gB8", start = "default:3", best = "hist_gB8"),
  gA14 = list(model = "gA14", start = "default:1", best = "hist_gA14"),
  gD3 = list(model = "gD3", start = "default:1", best = "hist_gD3"),
  gN1 = list(model = "gN1", start = "default:1", best = "hist_gN1"),
  bin16 = list(model = "cfs_binary", seed = 16L),
  bin27 = list(model = "cfs_binary", seed = 27L),
  acnonlin = list(model = "acnonlin", start = "default:1"))
which <- if (length(args) > 2) args[-(1:2)] else names(CASES)

# The penalised exact log likelihood at `x`, as the bench computes it.
reference <- function(fit, x) {
  obj <- .ctJuliaObjective(fit)
  ev <- jg(jc("ctsem_laplace_evaluate", obj, jv(x), gradient = FALSE))
  prior <- as.numeric(ev$value) - sum(as.numeric(ev$unit_loglik))
  PU <- as.matrix(jg(jc("bench_probe_units", obj, jv(x))))
  dmax <- ncol(PU) - 6L
  R <- if (dmax <= 3L) {
    vapply(seq_len(nrow(PU)), function(U) as.numeric(jg(jc("bench_probe_reference",
      obj, jv(x), as.integer(U), softcut = 3.5, nstiff = 9L, h = 0.1,
      halfwidth = 6.0))$value), 0)
  } else {
    vapply(seq_len(nrow(PU)), function(U) as.numeric(jg(jc("bench_probe_is", obj,
      jv(x), as.integer(U), n = 100000L, seed = 1L, inflate = 1.5))$value), 0)
  }
  quad <- tryCatch(as.numeric(jg(jc("ctsem_laplace_quadrature", obj, jv(x),
    nodes = 5L))$value), error = function(e) NA_real_)
  list(exact = sum(R) + prior, kind = if (dmax <= 3L) "exact-softcut3.5" else "is-t4-1.5",
    laplace = as.numeric(ev$value), quad5 = quad, prior = prior, units = R,
    below_one = sum(apply(PU[, 7:ncol(PU), drop = FALSE], 1,
      function(e) any(e < 1, na.rm = TRUE))))
}

# The Laplace objective along a direction through `x`: its value, its
# directional derivative from the gradient, and how many units the gated rule
# scored, at each step -- what a line search that finds no decrease meets.
JuliaConnectoR::juliaEval("lc_gated_units(o) = o.gated_units")
line_probe <- function(fit, x, direction, ts, floor = NULL) {
  obj <- .ctJuliaObjective(fit)
  if (!is.null(floor)) {
    own <- as.character(.ctJuliaOr(fit$model_spec$laplace$inner$floor, "total"))
    jc("ctsem_set_laplace_floor!", obj, floor)
    on.exit(jc("ctsem_set_laplace_floor!", obj, own), add = TRUE)
  }
  do.call(rbind, lapply(ts, function(t) {
    y <- x + t * direction
    ev <- jg(jc("ctsem_laplace_evaluate", obj, jv(y), gradient = TRUE))
    gated <- as.integer(JuliaConnectoR::juliaCall("lc_gated_units", obj))
    data.frame(t = t, value = as.numeric(ev$value),
      slope = sum(as.numeric(ev$gradient) * direction),
      converged = isTRUE(ev$converged), gated = gated)
  }))
}

REUSE <- Sys.getenv("VERIFY_REUSE", "")
FITS <- Sys.getenv("VERIFY_FITS", "")
if (nzchar(FITS)) dir.create(FITS, showWarnings = FALSE, recursive = TRUE)
SKIP_CONTINUE <- identical(tolower(Sys.getenv("VERIFY_SKIP_CONTINUE", "")), "true")
SKIP_PROBE <- identical(tolower(Sys.getenv("VERIFY_SKIP_PROBE", "")), "true")
reused_reference <- function(case, x) {
  if (!nzchar(REUSE)) return(NULL)
  f <- file.path(REUSE, paste0(case, ".rds"))
  if (!file.exists(f)) return(NULL)
  old <- readRDS(f)
  for (nm in names(old$estimates)) {
    y <- old$estimates[[nm]]
    if (length(y) == length(x) && max(abs(y - x)) == 0 &&
        !is.null(old$reference[[nm]])) return(old$reference[[nm]])
  }
  NULL
}

for (case in which) {
  cs <- CASES[[case]]
  if (is.null(cs)) { stamp("unknown case", case); next }
  f <- file.path(OUT, paste0(case, ".rds"))
  if (file.exists(f)) { stamp("have", case); next }
  stamp("CASE", case)
  rec <- list(case = case, spec = cs, host = Sys.info()[["nodename"]],
    tree = TREE, sha = tryCatch(system2("git", c("-C", TREE, "rev-parse", "HEAD"),
      stdout = TRUE), error = function(e) NA_character_))
  res <- tryCatch({
    if (identical(cs$model, "cfs_binary")) {
      d <- cfs_data(cs$seed); m <- cfs_model()
      # As the study does: the default start, drawn from the RNG state its
      # simulator left, so the same fit it made.
      inits <- NULL
      seed <- NULL
    } else {
      P <- bench_problem(cs$model)
      D <- bench_data(cs$model, "cfg", STORE)
      d <- D$data; m <- P$model()
      rec$data_md5 <- D$md5
      spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
        intoverpop = "laplace", fit = FALSE, cores = 1)))
      npar <- .ctBackendNpar(spec)
      rawnames <- suppressWarnings(.ctBackendRawParameterNames(list(model_spec = spec),
        npar))
      st <- bench_start(cs$start, npar, rawnames)
      inits <- st$inits; seed <- st$seed
    }
    saved <- if (nzchar(FITS)) file.path(FITS, paste0(case, "_fit.rds")) else ""
    if (nzchar(saved) && file.exists(saved)) {
      stored <- readRDS(saved)
      off <- stored$fit
      rec$fit_seconds <- stored$seconds
      rec$fit_loaded <- saved
      stamp("fit loaded from", saved)
    } else {
      t0 <- now()
      if (!is.null(seed)) set.seed(seed)
      off <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
        intoverpop = "laplace", cores = 1, inits = inits,
        optimcontrol = list(laplace_correct = FALSE, finishsamples = 100))))
      rec$fit_seconds <- now() - t0
      if (nzchar(saved)) {
        saveRDS(list(fit = off, seconds = rec$fit_seconds), paste0(saved, ".tmp"))
        file.rename(paste0(saved, ".tmp"), saved)
      }
    }
    rec$optim <- off$optim[intersect(names(off$optim), c("converged", "iterations",
      "f_calls", "g_calls", "stop_reason", "corrections", "hessian_evaluations",
      "stage_iterations", "saturated", "overshot"))]
    rec$certification <- off$uncertainty$certification[c("status", "certified", "gap")]
    stamp(sprintf("fit %.0fs, %s, stop %s, certification %s %.3g", rec$fit_seconds,
      if (isTRUE(off$optim$converged)) "converged" else "not converged",
      paste(off$optim$stop_reason, collapse = ","), rec$certification$status,
      as.numeric(rec$certification$gap)))
    t0 <- now()
    step <- suppressWarnings(.ctLaplaceAutoCorrect(off))
    rec$step_seconds <- now() - t0
    cont <- NULL
    if (!SKIP_CONTINUE) {
      t0 <- now()
      cont <- suppressWarnings(.ctLaplaceContinue(off, verbose = 1L))
      rec$continue_seconds <- now() - t0
      rec$continue_record <- cont$laplace$correction
    }
    rec$step_record <- step$laplace$correction
    stamp(sprintf("step %.0fs (%s)%s", rec$step_seconds,
      step$laplace$correction$status, if (is.null(cont)) "" else
        sprintf(", continue %.0fs (%s, %s, %d rounds)", rec$continue_seconds,
          cont$laplace$correction$status,
          cont$laplace$correction$continuation %||% "",
          as.integer(cont$laplace$correction$rounds %||% 0L))))
    se <- as.numeric(off$estimate$se)
    best <- if (!is.null(cs$best)) bench_start(paste0("stored:", cs$best),
      length(off$estimate$raw), NULL)$inits else NULL
    points <- list(laplace = as.numeric(off$estimate$raw),
      step = as.numeric(step$estimate$raw))
    fits <- list(laplace = off, step = step)
    if (!is.null(cont)) {
      points$continue <- as.numeric(cont$estimate$raw)
      fits$continue <- cont
    }
    rec$estimates <- points
    rec$se <- list(laplace = se,
      continue = if (!is.null(cont)) as.numeric(cont$estimate$se))
    # Distances in the Laplace fit's standard errors, over the coordinates that
    # have one: a saturated parameter reports an se near zero (5e-8 on gA14),
    # and a distance divided by it measures nothing.
    usable <- is.finite(se) & se > 1e-6
    dist_se <- function(a, b) if (any(usable) && length(a) == length(b))
      max(abs((a - b)[usable] / se[usable])) else NA_real_
    rec$best <- best
    rec$reference <- list()
    for (nm in names(points)) {
      x <- points[[nm]]
      same <- Filter(function(k) max(abs(points[[k]] - x)) == 0,
        names(rec$reference))
      again <- reused_reference(case, x)
      rec$reference[[nm]] <- if (length(same)) rec$reference[[same[1]]] else
        if (!is.null(again)) again else {
        t0 <- now(); r <- reference(off, x); r$seconds <- now() - t0; r }
      r <- rec$reference[[nm]]
      line <- data.frame(case = case, method = nm,
        reported = as.numeric(fits[[nm]]$estimate$loglik),
        logposterior = as.numeric(fits[[nm]]$estimate$logposterior),
        laplace = r$laplace, quad5 = r$quad5, exact = r$exact, kind = r$kind,
        best_known = if (!is.null(cs$best)) known$value[known$model == cs$model] else NA_real_,
        dist_best_se = if (!is.null(best)) dist_se(x, best) else NA_real_,
        dist_laplace_se = dist_se(x, points$laplace),
        below_one = r$below_one, stringsAsFactors = FALSE)
      utils::write.table(line, file.path(OUT, "summary.tsv"), sep = "\t",
        append = file.exists(file.path(OUT, "summary.tsv")),
        col.names = !file.exists(file.path(OUT, "summary.tsv")), row.names = FALSE,
        quote = FALSE)
      stamp(sprintf("  %-8s exact %.4f  quad5 %.4f  reported %.4f  dist to best %.3f se",
        nm, r$exact, r$quad5, line$reported, line$dist_best_se))
    }
    # gB8: why a resume stops. Along the Newton step the certification
    # predicts, the gated objective and the total-floor one side by side.
    if (!SKIP_PROBE && grepl("^gB8", case) &&
        !is.null(off$uncertainty$certification$step)) {
      direction <- as.numeric(off$uncertainty$certification$step)
      if (length(direction) == length(points$laplace) && any(direction != 0)) {
        ts <- seq(-1, 1, length.out = 41)
        rec$line_gated <- line_probe(off, points$laplace, direction, ts)
        rec$line_total <- line_probe(off, points$laplace, direction, ts, floor = "total")
        stamp("  line probe done")
      }
    }
    rec$trace <- off$optim$trace
    "ok"
  }, error = function(e) { stamp("ERROR", conditionMessage(e)); conditionMessage(e) })
  rec$status <- res
  saveRDS(rec, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  stamp("DONE", case)
}
stamp("ALLDONE")
