# One bench cell, one R process.
#
#   Rscript harness.R <id> <model> <data> <route> <start> <variant> <outdir>
#
# A cell is (id, model, data, route, start, variant) plus the build it runs
# against, which comes from the environment:
#
#   BENCH_LIB      library holding the build under test (dev1: ship.sh puts it
#                  at ~/dev/ctsemlib-bench-<label>, with a BENCH_BUILD file);
#   BENCH_TREE     or a source tree to devtools::load_all(compile = FALSE),
#                  for trying a cell locally (correctness only, never timing);
#   BENCH_LABEL    the build's label, recorded on every result;
#   BENCH_DATA     the data store (default ~/dev/ctsem-bench-data);
#   BENCH_TIMEOUT  cap on the fit, in seconds (default 7200);
#   BENCH_REFCAP   cap on the reference computations (default 5400);
#   BENCH_WARM     0 to skip the warm-up fit (default 1);
#   BENCH_REF      0 to skip the reference integrals (default 1);
#   BENCH_JULIA_BIN  the Julia bin directory, if not the versioned juliaup one.
#
# See README.md for the contract. What happens, in order:
#   1. the data, read from the store (simulated and stored on first use);
#   2. the contamination control: one gradient of this cell's objective at a
#      fixed point (raw zeros), after one untimed evaluation, three times;
#   3. a warm-up fit, the timed fit itself run once first, so compilation is
#      not in the timed fit (BENCH_WARM);
#   4. the fit, from a freshly built objective, under set.seed(seed) when the
#      start is default:<seed>, with every call to the optimiser, the Hessian
#      and the post-fit stages recorded by thin wrappers in the namespace;
#   5. the end point re-scored on a freshly built objective (cold inner modes);
#      for Laplace cells the 5-node quadrature there (ctLaplaceCheck) and the
#      exact reference where the cell's units allow one (BENCH_REF); where
#      the quadrature continuation reported its own Hessian, that Hessian
#      against cheaper ones (5b); and one overshoot probe at the end point,
#      timed (5c);
#   6. one .rds per cell in <outdir>, with a one-row data frame `row` of
#      scalars for summarise.R and the full record beside it.
#
# Nothing here changes what the build does: the wrappers forward every argument
# untouched and only time and count.

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
.num1 <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v)) v[1L] else NA_real_ }
.int1 <- function(x) { v <- suppressWarnings(as.integer(x)); if (length(v)) v[1L] else NA_integer_ }
.lgl1 <- function(x) if (is.null(x) || !length(x)) NA else isTRUE(as.logical(x[1L]))
.chr1 <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[1L])
.fin <- function(x) if (length(x) == 1L && is.finite(x)) x else NA_real_
.now <- function() proc.time()[["elapsed"]]
.loadavg <- function() {
  if (!file.exists("/proc/loadavg")) return(rep(NA_real_, 3))
  tryCatch(as.numeric(strsplit(readLines("/proc/loadavg", warn = FALSE)[1L],
    " ")[[1L]][1:3]), error = function(e) rep(NA_real_, 3))
}
.stamp <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }

args <- commandArgs(TRUE)
if (length(args) < 7L) stop("usage: Rscript harness.R <id> <model> <data> <route> <start> <variant> <outdir>")
CELL <- list(id = args[1], model = args[2], data = args[3], route = args[4],
  start = args[5], variant = args[6])
OUT <- args[7]
HERE <- Sys.getenv("BENCH_DIR", "")
if (!nzchar(HERE)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  HERE <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1L]))) else "."
}
LABEL <- Sys.getenv("BENCH_LABEL", "unlabelled")
STORE <- Sys.getenv("BENCH_DATA", file.path(Sys.getenv("HOME"), "dev", "ctsem-bench-data"))
CAP <- as.numeric(Sys.getenv("BENCH_TIMEOUT", "7200"))
REFCAP <- as.numeric(Sys.getenv("BENCH_REFCAP", "5400"))
WARM <- !identical(Sys.getenv("BENCH_WARM", "1"), "0")
DOREF <- !identical(Sys.getenv("BENCH_REF", "1"), "0")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "pids"), showWarnings = FALSE)
writeLines(as.character(Sys.getpid()), file.path(OUT, "pids", paste0(CELL$id, ".R")))
Sys.setenv(NOT_CRAN = "true", CTSEM_JULIA_AGREE = "yes", JULIA_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

rec <- list(cell = CELL, label = LABEL, host = Sys.info()[["nodename"]],
  pid = Sys.getpid(), started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  load_start = .loadavg(), ncpu = parallel::detectCores(), status = "started",
  error = NA_character_,
  harness = if (file.exists(file.path(HERE, "HARNESS_SHA")))
    readLines(file.path(HERE, "HARNESS_SHA"), warn = FALSE)[1L] else "working tree")
t_process <- .now()
save_rec <- function() {
  rec$load_end <<- .loadavg()
  rec$secs_process <<- .now() - t_process
  f <- file.path(OUT, paste0(CELL$id, ".rds"))
  tmp <- paste0(f, ".tmp")
  saveRDS(rec, tmp)
  invisible(file.rename(tmp, f))
}
fail <- function(stage, e) {
  rec$status <<- if (grepl("elapsed time limit", conditionMessage(e))) "timeout" else "error"
  rec$error <<- paste0(stage, ": ", conditionMessage(e))
  .stamp("FAILED", rec$status, rec$error)
  rec$row <<- NULL
  save_rec()
  quit(save = "no", status = 0)
}

# ---- the build ---------------------------------------------------------------
lib <- Sys.getenv("BENCH_LIB"); tree <- Sys.getenv("BENCH_TREE")
tryCatch({
  if (nzchar(lib)) {
    .libPaths(c(lib, .libPaths()))
    suppressMessages(library(ctsem))
    bf <- file.path(lib, "BENCH_BUILD")
    rec$build <- if (file.exists(bf)) readLines(bf, warn = FALSE) else NA_character_
  } else if (nzchar(tree)) {
    suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
    rec$build <- paste0("tree=", tree, " sha=", tryCatch(system2("git",
      c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA))
  } else stop("set BENCH_LIB (an installed build) or BENCH_TREE (a source tree)")
}, error = function(e) fail("load", e))
rec$ctsem_version <- as.character(utils::packageVersion("ctsem"))
ns <- asNamespace("ctsem")
.stamp("LOADED", CELL$id, "ctsem", rec$ctsem_version, "label", LABEL)

source(file.path(HERE, "cells.R"))
source(file.path(HERE, "starts.R"))
source(file.path(HERE, "variants.R"))

jbin <- Sys.getenv("BENCH_JULIA_BIN")
if (!nzchar(jbin)) {
  cand <- c(Sys.glob(file.path(Sys.getenv("HOME"), ".julia", "juliaup", "julia-1.12*", "bin")),
    Sys.glob(file.path(Sys.getenv("USERPROFILE"), ".julia", "juliaup", "julia-1.12*", "bin")))
  cand <- cand[file.exists(cand)]
  jbin <- if (length(cand)) sort(cand, decreasing = TRUE)[1L] else ""
}
tryCatch({
  if (nzchar(jbin)) suppressMessages(ctJuliaSetup(julia_bin = jbin, threads = 1L))
  else suppressMessages(ctJuliaSetup(threads = 1L))
  rec$julia_pid <- .int1(JuliaConnectoR::juliaEval("getpid()"))
  writeLines(as.character(rec$julia_pid), file.path(OUT, "pids", paste0(CELL$id, ".julia")))
  rec$julia_version <- as.character(JuliaConnectoR::juliaEval("string(VERSION)"))
}, error = function(e) fail("julia", e))
jeval <- function(txt) JuliaConnectoR::juliaEval(txt)
jcall <- function(fn, ...) JuliaConnectoR::juliaCall(paste0("ContinuousTimeSEM.", fn), ...)
# The bench's own Julia (reference integrals and the control timer), included
# into the engine module for this session only.
REFJL_ERR <- ""
REFJL <- tryCatch({
  jeval(sprintf("Core.eval(ContinuousTimeSEM, :(include(%s)))",
    deparse(normalizePath(file.path(HERE, "reference.jl"), winslash = "/"))))
  TRUE
}, error = function(e) { REFJL_ERR <<- conditionMessage(e); FALSE })
rec$reference_jl <- if (REFJL) "loaded" else REFJL_ERR
# JuliaConnectoR translates plain values (a Matrix, a Vector) on return and
# hands back a proxy for anything else; juliaGet only applies to the latter.
jget <- function(x) if (inherits(x, "JuliaProxy")) JuliaConnectoR::juliaGet(x) else x
jvec <- function(x) ctsem:::.ctJuliaNumericVector(as.numeric(x))
fresh_objectives <- function() {
  ce <- get(".ct_julia_cache", envir = ns)
  ce$objectives <- new.env(parent = emptyenv())
  invisible(NULL)
}
opcounts <- function() tryCatch({
  v <- as.numeric(jeval("collect(values(ContinuousTimeSEM.ctsem_opcounts()))"))
  names(v) <- as.character(jeval("string.(collect(keys(ContinuousTimeSEM.ctsem_opcounts())))"))
  v
}, error = function(e) NULL)

# ---- the problem ---------------------------------------------------------------
tryCatch({
  P <- bench_problem(CELL$model)
  if (!CELL$route %in% P$routes) stop("route ", CELL$route, " is not one of ",
    paste(P$routes, collapse = ", "), " for ", CELL$model)
  D <- bench_data(CELL$model, CELL$data, STORE)
  dat <- D$data
  rec$data_file <- D$file; rec$data_md5 <- D$md5
  model <- P$model()
  V <- bench_variant(CELL$variant)
  rec$variant <- V[c("name", "optimcontrol", "args", "objective", "fit")]
  # Only what the cell sets: an argument left at its default is not passed, so
  # nothing that asks missing() inside ctFit sees a different call from the
  # one a user would write.
  fitargs <- function(d, inits, oc_extra = list()) {
    a <- list(datalong = d, model = model, backend = "julia",
      intoverpop = CELL$route, cores = 1L)
    if (!is.null(inits)) a$inits <- inits
    oc <- utils::modifyList(V$optimcontrol, oc_extra)
    if (length(oc)) a$optimcontrol <- oc
    utils::modifyList(a, V$args)
  }
  spec <- suppressWarnings(suppressMessages(do.call(ctFit,
    c(fitargs(dat, NULL), list(fit = FALSE)))))
  npar <- ctsem:::.ctBackendNpar(spec)
  rawnames <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
  S <- bench_start(CELL$start, npar, rawnames)
  rec$npar <- npar; rec$rawnames <- rawnames; rec$start_vector <- S$inits
  rec$nsubjects <- length(unique(dat[[P$idcol]])); rec$nrows <- nrow(dat)
  laplace <- !is.null(spec$laplace)
  rec$laplace <- laplace
}, error = function(e) fail("problem", e))
.stamp("DATA", nrow(dat), "rows", rec$nsubjects, "groups; npar", npar, "; md5", D$md5)

# ---- 2. contamination control ----------------------------------------------------
# One gradient at raw zeros, the same point for every variant and start of this
# model, after one untimed evaluation that pays the compilation. Nothing an
# optimiser change does can move it, so if it moves between cells the machine
# moved.
# Timed inside Julia (bench_time_gradient in reference.jl) in batches of about
# a quarter second, so neither the bridge's round trip nor the timer's
# resolution is in it; the seconds are per evaluation.
ctrl <- tryCatch({
  at <- rep(0, npar)
  t0 <- .now(); invisible(ctJuliaEvaluate(spec, at, gradient = TRUE)); first <- .now() - t0
  if (!isTRUE(REFJL)) stop("reference.jl did not load: ", REFJL_ERR)
  tg <- jget(jcall("bench_time_gradient", ctsem:::.ctJuliaObjective(spec), jvec(at)))
  list(first = first, reps = as.numeric(tg$seconds), batch = .int1(tg$batch))
}, error = function(e) list(first = NA_real_, reps = rep(NA_real_, 3),
  error = conditionMessage(e)))
rec$ctrl <- ctrl
.stamp(sprintf("CONTROL gradient at zeros: first call %.3fs; per evaluation %s (batches of %s)",
  ctrl$first, paste(sprintf("%.3g", ctrl$reps), collapse = " "), format(ctrl$batch)),
  if (!is.null(ctrl$error)) paste("ERROR", ctrl$error) else "")
fresh_objectives()

if (isTRUE(V$control_only)) {
  rec$status <- "control"
  rec$row <- data.frame(id = CELL$id, model = CELL$model, data = CELL$data,
    route = CELL$route, start = CELL$start, variant = CELL$variant,
    label = LABEL, host = rec$host, status = "control",
    ctrl1 = ctrl$reps[1], ctrl2 = ctrl$reps[2], ctrl3 = ctrl$reps[3],
    load1_start = rec$load_start[1], stringsAsFactors = FALSE)
  save_rec(); .stamp("DONE", CELL$id, "control"); quit(save = "no", status = 0)
}

# ---- instrumentation ---------------------------------------------------------------
# Thin wrappers in the namespace: each forwards its arguments untouched and
# records what ran. Names that a build no longer has are skipped and the record
# says which were in place, so a missing stage table is never mistaken for a
# stage that did not run.
REC <- new.env()
REC$stages <- list(); REC$runs <- list(); REC$hessians <- list(); REC$post <- list()
REC$cont <- NULL; REC$continue <- NULL
REC$current <- 0L
wrap <- function(name, factory) {
  if (!exists(name, envir = ns, inherits = FALSE)) return(FALSE)
  orig <- get(name, envir = ns)
  w <- factory(orig)
  unlockBinding(name, ns); assign(name, w, envir = ns); lockBinding(name, ns)
  TRUE
}
stage_kind <- function(label, caller) {
  if (identical(label, "prior warm-up")) return("warmup")
  if (identical(caller, "optimise")) return("resume")
  if (grepl("Restart", caller)) return("restart")
  if (grepl("Sample", caller)) return("sample-placement")
  "main"
}
scalars <- function(r, keep) {
  out <- list()
  for (k in keep) if (!is.null(r[[k]]) && length(r[[k]]) >= 1L && length(r[[k]]) <= 20L &&
    (is.numeric(r[[k]]) || is.logical(r[[k]]) || is.character(r[[k]]))) out[[k]] <- r[[k]]
  out
}
RESULT_KEYS <- c("iterations", "f_calls", "g_calls", "maximum_loglik", "converged",
  "stopped_by_gap", "stopped_by_stall", "stalled", "overshot", "saturated",
  "predicted_gain", "last_gain", "gradient_norm", "newton_steps", "newton_hessians",
  "newton_subset_hessians", "batch_sizes", "batch_iterations", "stall_escapes",
  "stall_triggers", "gap_tol", "linesearch", "chunks", "overshoot_gain",
  "stopped_by_stall", "stall_gain")
rec$instrumented <- c(
  .ctJuliaOptimise = wrap(".ctJuliaOptimise", function(orig) function(...) {
    a <- list(...)
    # The function whose frame made the call: sys.parent(), not sys.call(-1),
    # because a call made inside suppressWarnings() is evaluated as a promise
    # several frames below its caller.
    pf <- sys.parent()
    caller <- if (pf > 0L) substr(paste(deparse(sys.call(pf)[[1L]]),
      collapse = ""), 1L, 60L) else "top"
    label <- if (!is.null(a$progress_label)) as.character(a$progress_label) else ""
    oc <- a$optimcontrol
    k <- length(REC$stages) + 1L
    REC$stages[[k]] <- list(kind = stage_kind(label, caller), caller = caller,
      label = label, maxiter = .int1(a$maxiter %||% oc$maxiter),
      innergaptol = .num1(oc$innergaptol), g_tol = .num1(oc$g_tol), secs = NA_real_)
    outer <- REC$current; REC$current <- k
    on.exit(REC$current <- outer, add = TRUE)
    t0 <- .now()
    r <- orig(...)
    REC$stages[[k]]$secs <- .now() - t0
    REC$stages[[k]]$result <- scalars(r, RESULT_KEYS)
    r
  }),
  # Every engine optimisation passes through here, including the stall-escape
  # stages inside one .ctJuliaOptimise call, whose counts that call's own result
  # does not sum. Recorded only when the value looks like an optimiser result.
  .ctBackendWithMaxChunks = wrap(".ctBackendWithMaxChunks", function(orig) function(chunks, expr) {
    t0 <- .now()
    r <- orig(chunks, expr)
    if (is.list(r) && !is.null(r$minimizer) && !is.null(r$iterations)) {
      # Where the run's overshoot probe ran (it runs at the run's end point,
      # the minimizer), and what it pulled back when it found a gain.
      m <- suppressWarnings(as.numeric(r$minimizer))
      op <- suppressWarnings(as.integer(r$overshoot_parameters))
      op <- op[!is.na(op) & op >= 1L & op <= length(m)]
      REC$runs[[length(REC$runs) + 1L]] <- c(list(stage = REC$current,
        secs = .now() - t0, npar = length(m),
        maxabs_raw = if (length(m) && all(is.finite(m))) max(abs(m)) else NA_real_,
        overshoot_maxabs = if (length(op)) max(abs(m[op])) else NA_real_,
        overshoot_n = length(op)), scalars(r, RESULT_KEYS))
    }
    r
  }),
  .ctBackendHessianAt = wrap(".ctBackendHessianAt", function(orig) function(...) {
    t0 <- .now(); r <- orig(...)
    REC$hessians[[length(REC$hessians) + 1L]] <- list(where = "certification",
      stage = REC$current, secs = .now() - t0, ok = is.matrix(r), computed = TRUE)
    r
  }),
  # This one returns the certification's matrix unchanged when it was
  # evaluated at the same point, which costs nothing and is not a Hessian
  # computed; `computed` tells the two apart.
  .ctBackendHessian = wrap(".ctBackendHessian", function(orig) function(...) {
    t0 <- .now(); r <- orig(...)
    first <- list(...)[[1L]]
    stored <- tryCatch(first$uncertainty$hessian, error = function(e) NULL)
    REC$hessians[[length(REC$hessians) + 1L]] <- list(where = "uncertainty",
      stage = REC$current, secs = .now() - t0, ok = !is.null(r),
      computed = !(is.matrix(r) && is.matrix(stored) && identical(r, stored)))
    r
  }),
  .ctLaplaceAutoCorrect = wrap(".ctLaplaceAutoCorrect", function(orig) function(...) {
    t0 <- .now(); r <- orig(...)
    REC$post$laplace_correct <- (REC$post$laplace_correct %||% 0) + .now() - t0
    r
  }),
  # The fit calls ctFitUncertainty since juliaFit 908b068d, and
  # ctOptimUncertainty is then an alias that calls it, so the old name is
  # wrapped only on a build without the new one (else its time counts twice).
  # NA in `instrumented` is "deliberately not wrapped".
  ctFitUncertainty = wrap("ctFitUncertainty", function(orig) function(...) {
    t0 <- .now(); r <- orig(...)
    REC$post$uncertainty <- (REC$post$uncertainty %||% 0) + .now() - t0
    r
  }),
  ctOptimUncertainty = if (exists("ctFitUncertainty", envir = ns, inherits = FALSE)) NA else
    wrap("ctOptimUncertainty", function(orig) function(...) {
      t0 <- .now(); r <- orig(...)
      REC$post$uncertainty <- (REC$post$uncertainty %||% 0) + .now() - t0
      r
    }),
  .ctBackendCorrectResult = wrap(".ctBackendCorrectResult", function(orig) function(...) {
    t0 <- .now(); r <- orig(...)
    REC$post$certification <- (REC$post$certification %||% 0) + .now() - t0
    r
  }),
  # The quadrature continuation replaces the fit's uncertainty with the Hessian
  # of its own objective. What it replaced is kept here -- the Laplace Hessian
  # at the Laplace optimum and the standard errors from it -- so section 5b can
  # say whether the replacement changed the standard errors.
  # Timed too: it is the default correction since juliaFit 55caeeb7, and the
  # fit calls it directly rather than through .ctLaplaceAutoCorrect.
  .ctLaplaceContinue = wrap(".ctLaplaceContinue", function(orig) function(fit, ...) {
    REC$continue <- list(hessian = fit$uncertainty$hessian,
      evaluated_at = as.numeric(fit$uncertainty$evaluated_at %||% NA_real_),
      raw = as.numeric(fit$estimate$raw), se = as.numeric(fit$estimate$se %||% NA_real_))
    t0 <- .now(); r <- orig(fit, ...)
    REC$post$laplace_correct <- (REC$post$laplace_correct %||% 0) + .now() - t0
    r
  }),
  # Called only when the continuation's own Hessian becomes the fit's; its
  # `cont` is the continuation, placed at the final estimate, which 5b uses to
  # time that Hessian and its cheaper alternatives.
  .ctLaplaceContinueLpg = wrap(".ctLaplaceContinueLpg", function(orig) function(module, cont) {
    REC$cont <- cont
    orig(module, cont)
  })
)

# ---- 3. warm-up ----------------------------------------------------------------------
# The timed fit itself, run once first: same data, arguments, start and seed.
# Until 2026-09-26 this was a fit on a slice of the subjects capped at
# `maxiter = 5`, which reaches only the stages a five-iteration fit reaches, so
# every later one -- the Newton finish, the certification, the quadrature
# continuation -- compiled inside the timed fit, and which of them did changed
# with the build: gA1 read 6.9 s at one build and 113 s at the next with
# identical counts, and the same fit run twice in one session took 207.9 s and
# then 7.4 s. Only the fit itself is guaranteed to reach exactly the stages the
# timed fit will. It doubles a cell's fitting time. The two fits must also
# agree, since nothing is meant to carry from one fit to the next, and the row
# says whether they did (`warm_dx`, `warm_dll`). BENCH_WARM=0 skips it.
if (WARM && isTRUE(V$fit)) {
  t0 <- .now()
  w <- tryCatch({
    set.seed(S$seed)
    setTimeLimit(elapsed = CAP, transient = TRUE)
    wf <- suppressWarnings(suppressMessages(do.call(ctFit, fitargs(dat, S$inits))))
    setTimeLimit()
    list(result = "ok", raw = as.numeric(wf$estimate$raw),
      loglik = .num1(wf$estimate$loglik), logposterior = .num1(wf$estimate$logposterior),
      iterations = .int1(wf$optim$iterations), g_calls = .int1(wf$optim$g_calls))
  }, error = function(e) { setTimeLimit(); list(result = conditionMessage(e)) })
  rec$warm <- c(list(secs = .now() - t0), w)
  REC$stages <- list(); REC$runs <- list(); REC$hessians <- list(); REC$post <- list()
  REC$cont <- NULL; REC$continue <- NULL
  .stamp(sprintf("WARM %.1fs %s", rec$warm$secs, substr(w$result, 1, 120)))
  # A warm-up that ran out of time says the fit will too, and its Julia may
  # still be computing with nobody to answer: stop it, as a timed-out fit is.
  if (grepl("elapsed time limit", w$result)) {
    rec$status <- "timeout"
    rec$error <- paste0("warm-up fit: ", w$result)
    save_rec()
    if (is.finite(rec$julia_pid %||% NA)) try(tools::pskill(rec$julia_pid), silent = TRUE)
    .stamp("DONE", CELL$id, rec$status)
    quit(save = "no", status = 0)
  }
}
fresh_objectives()

# ---- 4. the fit ------------------------------------------------------------------------
warns <- character(); msgs <- character()
fit <- NULL
if (isTRUE(V$fit)) {
  try(jeval("ContinuousTimeSEM.ctsem_reset_opcounts!()"), silent = TRUE)
  set.seed(S$seed)
  rec$load_fit_start <- .loadavg()
  t0 <- .now()
  fit <- tryCatch({
    setTimeLimit(elapsed = CAP, transient = TRUE)
    withCallingHandlers(do.call(ctFit, fitargs(dat, S$inits)),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) {
        if (length(msgs) < 200L) msgs <<- c(msgs, substr(conditionMessage(m), 1L, 400L))
        invokeRestart("muffleMessage")
      })
  }, error = function(e) e)
  setTimeLimit()
  rec$secs_fit <- .now() - t0
  rec$load_fit_end <- .loadavg()
  rec$opcounts <- opcounts()
  rec$warnings <- warns; rec$messages <- msgs
  rec$stages <- REC$stages; rec$runs <- REC$runs; rec$hessian_calls <- REC$hessians
  rec$post_secs <- REC$post
  rec$continue_before <- REC$continue
  if (inherits(fit, "error")) {
    rec$status <- if (grepl("elapsed time limit", conditionMessage(fit))) "timeout" else "error"
    rec$error <- paste0("fit: ", conditionMessage(fit))
    .stamp("FIT", rec$status, rec$error)
    fit <- NULL
    # A timeout can leave this cell's Julia mid-computation with nobody to
    # answer. It is provably ours (its pid was recorded at spawn), so stop it
    # rather than let it load the machine under the next cell.
    if (identical(rec$status, "timeout") && is.finite(rec$julia_pid %||% NA)) {
      save_rec()
      try(tools::pskill(rec$julia_pid), silent = TRUE)
      .stamp("DONE", CELL$id, rec$status)
      quit(save = "no", status = 0)
    }
  } else {
    rec$status <- "ok"
    .stamp(sprintf("FIT %.1fs, %d warnings; iterations %s, loglik %.4f, certification %s",
      rec$secs_fit, length(warns), format(fit$optim$iterations),
      .num1(fit$estimate$loglik), format(fit$uncertainty$certification$status %||% NA)))
  }
}

# ---- what the fit says about itself -----------------------------------------------
# Where the references are evaluated: the fit's estimate, or the start itself
# for a cell that does not fit. A failed fit has no end point to score.
est <- if (!is.null(fit)) as.numeric(fit$estimate$raw) else
  if (!isTRUE(V$fit)) S$inits else NULL
if (!isTRUE(V$fit)) rec$status <- "evalonly"
if (!is.null(fit)) {
  o <- fit$optim
  cert <- fit$uncertainty$certification
  corr <- fit$laplace$correction
  rec$optim <- o[setdiff(names(o), c("trace", "gradient", "hessian_profile"))]
  rec$trace <- o$trace
  rec$estimate <- list(raw = est, se = as.numeric(fit$estimate$se %||% NA_real_),
    loglik = .num1(fit$estimate$loglik), logposterior = .num1(fit$estimate$logposterior),
    loglik_method = .chr1(fit$estimate$loglik_method),
    loglik_laplace = .num1(fit$estimate$loglik_laplace))
  rec$certification <- if (is.list(cert)) cert[intersect(names(cert), c("status",
    "verdict", "gap", "lambda_min", "ntrusted", "nflat", "nnegative",
    "residual_gain", "tolerance", "certified"))] else NULL
  rec$laplace_block <- if (!is.null(fit$laplace)) list(
    correction = corr[setdiff(names(corr), c("draws"))],
    conditioning = fit$laplace$conditioning, floor = fit$laplace$floor,
    gated_units = fit$laplace$gated_units) else NULL
  rec$identifiability <- tryCatch(list(nweak = .int1(fit$identifiability$nweak),
    weak = fit$identifiability$weak), error = function(e) NULL)
}

# ---- 5. re-score and references ----------------------------------------------------
refs <- list()
if (DOREF && !is.null(est) && length(est) == npar && all(is.finite(est))) {
  t0 <- .now()
  refs <- tryCatch({
    setTimeLimit(elapsed = REFCAP, transient = TRUE)
    out <- list()
    fresh_objectives()
    obj <- ctsem:::.ctJuliaObjective(if (!is.null(fit)) fit else spec)
    if (laplace) {
      ev <- jget(jcall("ctsem_laplace_evaluate", obj, jvec(est), gradient = FALSE))
      out$rescored <- .num1(ev$value)
      out$prior_term <- .num1(ev$value) - sum(as.numeric(ev$unit_loglik))
      qc <- tryCatch({
        fresh_objectives()
        x <- if (!is.null(fit)) fit else NULL
        if (is.null(x)) stop("no fit")
        ch <- ctLaplaceCheck(x, nodes = 5L, correction = FALSE)
        list(quadrature = .num1(ch$quadrature), laplace = .num1(ch$laplace),
          gap = .num1(ch$gap))
      }, error = function(e) {
        # evalonly cells have no fit: the same two calls ctLaplaceCheck makes.
        q <- tryCatch(.num1(jget(jcall("ctsem_laplace_quadrature", obj, jvec(est),
          nodes = 5L))$value), error = function(e2) NA_real_)
        list(quadrature = q, laplace = out$rescored, gap = q - out$rescored,
          note = conditionMessage(e))
      })
      out$quad5 <- qc
      # Units, their dimensions and curvature; then the reference.
      if (!isTRUE(REFJL)) stop("reference.jl did not load: ", REFJL_ERR)
      fresh_objectives()
      obj <- ctsem:::.ctJuliaObjective(if (!is.null(fit)) fit else spec)
      PU <- as.matrix(jget(jcall("bench_probe_units", obj, jvec(est))))
      nunit <- nrow(PU); dmax <- ncol(PU) - 6L
      eig <- PU[, 7:ncol(PU), drop = FALSE]
      out$units <- list(n = nunit, dmax = dmax,
        lambda_min = suppressWarnings(min(eig, na.rm = TRUE)),
        below_one = sum(apply(eig, 1, function(e) any(e < 1, na.rm = TRUE))))
      # Seeds that reach the same optimum share its reference: the integral is
      # a property of the build, the data and the point, and at an optimum a
      # raw difference of 1e-4 moves it by far less than anything reported.
      # Cached per results directory, so never across builds.
      cachedir <- file.path(OUT, "refcache")
      dir.create(cachedir, showWarnings = FALSE)
      hit <- NULL
      for (f in list.files(cachedir, "\\.rds$", full.names = TRUE)) {
        cc <- tryCatch(readRDS(f), error = function(e) NULL)
        if (!is.null(cc) && identical(cc$model, CELL$model) &&
            identical(cc$data_md5, D$md5) && identical(cc$route, CELL$route) &&
            length(cc$est) == length(est) && max(abs(cc$est - est)) < 1e-4) {
          hit <- cc; break
        }
      }
      if (!is.null(hit)) {
        out$exact <- hit$exact; out$exact_alt <- hit$exact_alt
        out$exact_kind <- hit$exact_kind; out$unit_ref <- hit$unit_ref
        out$is_check <- hit$is_check
        out$exact_from <- hit$from
        out$exact <- out$exact - hit$prior_term + out$prior_term
      } else if (P$reference && dmax <= 3L) {
        R <- vapply(seq_len(nunit), function(U) .num1(jget(jcall("bench_probe_reference",
          obj, jvec(est), as.integer(U), softcut = 3.5, nstiff = 9L, h = 0.1,
          halfwidth = 6.0))$value), 0)
        out$exact <- sum(R) + out$prior_term
        out$exact_kind <- "exact-softcut3.5"
        out$unit_ref <- R
        # Agreement with importance sampling on the three units where the
        # reference departs most from the unit's own Laplace term.
        pick <- utils::head(order(-abs(R - PU[, 6])), min(3L, nunit))
        out$is_check <- do.call(rbind, lapply(pick, function(U) {
          r <- jget(jcall("bench_probe_is", obj, jvec(est), as.integer(U),
            n = 100000L, seed = 1L, inflate = 1.5))
          data.frame(unit = U, reference = R[U], is = .num1(r$value),
            is_se = .num1(r$se), ess = .num1(r$ess))
        }))
      } else if (P$reference) {
        R1 <- vapply(seq_len(nunit), function(U) .num1(jget(jcall("bench_probe_is",
          obj, jvec(est), as.integer(U), n = 100000L, seed = 1L, inflate = 1.5))$value), 0)
        R2 <- vapply(seq_len(nunit), function(U) .num1(jget(jcall("bench_probe_is",
          obj, jvec(est), as.integer(U), n = 100000L, seed = 2L, inflate = 2.5))$value), 0)
        out$exact <- sum(R1) + out$prior_term
        out$exact_alt <- sum(R2) + out$prior_term
        out$exact_kind <- "is-t4-1.5"
        out$unit_ref <- R1; out$unit_ref_alt <- R2
      }
      if (is.null(hit) && !is.null(out$exact)) {
        f <- file.path(cachedir, paste0(CELL$id, ".rds"))
        saveRDS(list(model = CELL$model, data_md5 = D$md5, route = CELL$route,
          est = est, exact = out$exact, exact_alt = out$exact_alt,
          exact_kind = out$exact_kind, unit_ref = out$unit_ref,
          is_check = out$is_check, prior_term = out$prior_term, from = CELL$id),
          paste0(f, ".tmp"))
        file.rename(paste0(f, ".tmp"), f)
      }
    } else {
      ev <- jget(jcall("ctsem_evaluate", obj, jvec(est), gradient = FALSE,
        contributions = FALSE))
      out$rescored <- .num1(ev$value)
    }
    out
  }, error = function(e) list(error = conditionMessage(e)))
  setTimeLimit()
  refs$secs <- .now() - t0
  .stamp(sprintf("REFERENCE %.1fs rescored %s quad %s exact %s %s", refs$secs,
    format(refs$rescored), format(refs$quad5$quadrature), format(refs$exact),
    if (!is.null(refs$error)) paste("ERROR", refs$error) else ""))
}
rec$refs <- refs

# ---- 5b. the continuation's Hessian, by each scheme --------------------------------
# Only where the quadrature continuation reached its fixed point and its own
# Hessian became the fit's. The continuation's Hessian at the estimate by each
# scheme the build has -- `exact` (the flagged units by the Louis identity,
# the default since the quadhess job; NA on a build without it), `forward`
# (npar + 1 hybrid gradients) and `central` (2 npar) -- then the Laplace
# Hessian at the continuation's estimate (2 npar Laplace gradients) and the
# one the fit already held at the Laplace optimum (free). Standard errors from
# each by one rule -- the inverse on the directions whose curvature is above
# 1e-8 of the largest -- and each one's largest relative difference from the
# central scheme's; `reported` checks the fit's standard errors against its
# own Hessian. Timed here, outside the fit, one after another, cheapest
# first, each on its own so that one running out of time keeps the others.
bench_se <- function(H) {
  if (!is.matrix(H) || any(!is.finite(H))) return(NULL)
  I <- -(H + t(H)) / 2
  e <- eigen(I, symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values)
  if (!any(keep)) return(NULL)
  V <- e$vectors[, keep, drop = FALSE]
  list(se = sqrt(pmax(0, rowSums(sweep(V^2, 2, e$values[keep], "/")))),
    negative = sum(e$values < -1e-8 * max(abs(e$values))), kept = sum(keep))
}
se_rel <- function(a, b) {
  if (is.null(a) || is.null(b) || length(a$se) != length(b$se)) return(NA_real_)
  ok <- b$se > 0 & is.finite(a$se) & is.finite(b$se)
  if (!any(ok)) NA_real_ else max(abs(a$se[ok] / b$se[ok] - 1))
}
cse <- list()
Hq <- if (!is.null(fit)) fit$laplace$correction$hessian else NULL
if (!is.null(REC$cont) && is.matrix(Hq) && !is.null(REC$continue)) {
  cse <- tryCatch({
    setTimeLimit(elapsed = REFCAP, transient = TRUE)
    x <- as.numeric(fit$estimate$raw)
    # Not `errors`: `cse$error` below would partial-match it.
    out <- list(delta_se = suppressWarnings(as.numeric(fit$laplace$correction$delta_se)),
      scheme_failures = list())
    scheme <- function(name) {
      t0 <- .now()
      H <- tryCatch(matrix(as.numeric(jget(jcall("ctsem_laplace_continuation_hessian",
        REC$cont, jvec(x), scheme = name))), npar, npar),
        error = function(e) { out$scheme_failures[[name]] <<- conditionMessage(e); NULL })
      # A build without `scheme = "forward"` has the bench's own copy of it.
      if (is.null(H) && name == "forward") H <- tryCatch(matrix(as.numeric(jget(jcall(
        "bench_continuation_hessian_forward", REC$cont, jvec(x)))), npar, npar),
        error = function(e) NULL)
      out[[paste0("secs_", name)]] <<- if (is.null(H)) NA_real_ else .now() - t0
      H
    }
    He <- scheme("exact")
    out$exact_repro <- if (is.matrix(He)) max(abs(He - Hq)) else NA_real_
    Hf <- scheme("forward")
    Hc <- scheme("central")
    t0 <- .now()
    Hlx <- tryCatch(ctsem:::.ctBackendHessianAt(fit$model_spec, x), error = function(e) NULL)
    out$secs_laplace_x <- .now() - t0
    Hl0 <- REC$continue$hessian
    sc <- bench_se(Hc)
    sq <- bench_se(Hq)
    out$se <- list(central = sc$se, exact = bench_se(He)$se, forward = bench_se(Hf)$se,
      quadrature = sq$se, laplace_x = bench_se(Hlx)$se, laplace_est = bench_se(Hl0)$se,
      reported = as.numeric(fit$estimate$se %||% NA_real_),
      reported_before = REC$continue$se)
    out$negative <- c(central = sc$negative %||% NA, exact = bench_se(He)$negative %||% NA,
      forward = bench_se(Hf)$negative %||% NA, quadrature = sq$negative %||% NA,
      laplace_x = bench_se(Hlx)$negative %||% NA, laplace_est = bench_se(Hl0)$negative %||% NA)
    out$rel <- c(exact = se_rel(bench_se(He), sc), forward = se_rel(bench_se(Hf), sc),
      laplace_x = se_rel(bench_se(Hlx), sc), laplace_est = se_rel(bench_se(Hl0), sc),
      reported = se_rel(list(se = out$se$reported), sq))
    out$hessians <- list(quadrature = Hq, exact = He, forward = Hf, central = Hc,
      laplace_x = Hlx, laplace_est = Hl0)
    out
  }, error = function(e) list(error = conditionMessage(e)))
  setTimeLimit()
  .stamp(sprintf(paste0("CONTINUATION HESSIAN exact %.1fs forward %.1fs central %.1fs ",
    "laplace-at-x %.1fs; max rel se diff from central: exact %s forward %s laplace_x %s ",
    "laplace_est %s %s"),
    cse$secs_exact %||% NA, cse$secs_forward %||% NA, cse$secs_central %||% NA,
    cse$secs_laplace_x %||% NA,
    format(signif(cse$rel[["exact"]] %||% NA, 3)), format(signif(cse$rel[["forward"]] %||% NA, 3)),
    format(signif(cse$rel[["laplace_x"]] %||% NA, 3)),
    format(signif(cse$rel[["laplace_est"]] %||% NA, 3)),
    if (!is.null(cse$error)) paste("ERROR", cse$error) else ""))
}
rec$continuation_se <- cse
REC$cont <- NULL

# ---- 5c. the overshoot probe's cost at the end point ------------------------------
# What one run's verdict probe costs here: `ctsem_pullback` runs the same probe
# (`_ctsem_overshot`, the default mode) at the estimate, plus one evaluation
# there. Once per cell, outside the fit, on a fresh objective. Times the number
# of engine runs it estimates the probe's share of the fit; the `noprobe`
# variant measures that directly on a few cells.
probe <- list()
if (!is.null(fit) && length(est) == npar && all(is.finite(est))) {
  probe <- tryCatch({
    setTimeLimit(elapsed = REFCAP, transient = TRUE)
    fresh_objectives()
    obj <- ctsem:::.ctJuliaObjective(fit)
    t0 <- .now()
    pb <- jget(jcall("ctsem_pullback", obj, jvec(est), tolerance = 1e-6))
    list(secs = .now() - t0, found = .lgl1(pb$found), gain = .num1(pb$gain))
  }, error = function(e) list(error = conditionMessage(e)))
  setTimeLimit()
  .stamp(sprintf("PROBE at the estimate %.2fs, found %s", probe$secs %||% NA,
    format(probe$found %||% NA)))
}
rec$probe <- probe

# ---- the row ---------------------------------------------------------------------------
st <- rec$stages %||% list()
kinds <- vapply(st, function(s) s$kind, "")
stagecount <- function(kind, key) sum(vapply(st[kinds == kind], function(s)
  .num1(s$result[[key]] %||% NA_real_), 0), na.rm = TRUE)
runs <- rec$runs %||% list()
runsum <- function(key) sum(vapply(runs, function(r) .num1(r[[key]] %||% NA_real_), 0), na.rm = TRUE)
o <- rec$optim %||% list()
cert <- rec$certification %||% list()
corr <- rec$laplace_block$correction %||% list()
compact <- paste(vapply(st, function(s) sprintf("%s:%s/%s/%s", s$kind,
  format(s$result$iterations %||% NA), format(s$result$f_calls %||% NA),
  format(s$result$g_calls %||% NA)), ""), collapse = ";")
rec$row <- data.frame(
  id = CELL$id, model = CELL$model, data = CELL$data, route = CELL$route,
  start = CELL$start, variant = CELL$variant, objective = V$objective,
  label = LABEL, build = paste(rec$build %||% NA, collapse = " "),
  host = rec$host, status = rec$status, error = rec$error,
  npar = npar, nsubjects = rec$nsubjects, nrows = rec$nrows,
  data_md5 = rec$data_md5,
  secs_fit = rec$secs_fit %||% NA_real_, secs_warm = rec$warm$secs %||% NA_real_,
  secs_ref = refs$secs %||% NA_real_,
  secs_optimise = sum(vapply(st, function(s) .num1(s$secs), 0), na.rm = TRUE),
  secs_certify = .num1(rec$post_secs$certification %||% NA),
  secs_uncertainty = .num1(rec$post_secs$uncertainty %||% NA),
  secs_lapcorrect = .num1(rec$post_secs$laplace_correct %||% NA),
  load1_start = rec$load_start[1], load1_fit_start = (rec$load_fit_start %||% NA)[1],
  load1_fit_end = (rec$load_fit_end %||% NA)[1],
  ctrl1 = ctrl$reps[1], ctrl2 = ctrl$reps[2], ctrl3 = ctrl$reps[3],
  converged = .lgl1(o$converged),
  cert_status = .chr1(cert$status), cert_gap = .num1(cert$gap),
  iterations = .int1(o$iterations), stage_iterations = .int1(o$stage_iterations),
  f_calls = .int1(o$f_calls), g_calls = .int1(o$g_calls),
  engine_runs = length(runs), engine_iterations = runsum("iterations"),
  engine_f_calls = runsum("f_calls"), engine_g_calls = runsum("g_calls"),
  warm_iterations = stagecount("warmup", "iterations"),
  warm_f_calls = stagecount("warmup", "f_calls"), warm_g_calls = stagecount("warmup", "g_calls"),
  main_iterations = stagecount("main", "iterations"),
  resume_iterations = stagecount("resume", "iterations"),
  restart_iterations = stagecount("restart", "iterations"),
  n_resumes = sum(kinds == "resume"), n_restarts = sum(kinds == "restart"),
  stages = compact,
  warmup_claimed = .lgl1(o$carefulfit), warmup_ran = any(kinds == "warmup"),
  hessians_cert = .int1(o$hessian_evaluations),
  hessians_computed = sum(vapply(rec$hessian_calls %||% list(), function(h)
    isTRUE(h$computed), logical(1))),
  newton_steps = runsum("newton_steps"), newton_hessians = runsum("newton_hessians"),
  newton_subset_hessians = runsum("newton_subset_hessians"),
  batch_sizes = paste(o$batch_sizes %||% NA, collapse = ","),
  stopped_by_gap = .lgl1(o$stopped_by_gap), stopped_by_stall = .lgl1(o$stopped_by_stall),
  stall_escapes = .int1(o$stall_escapes), overshot = .lgl1(o$overshot),
  saturated = .lgl1(o$saturated), corrections = length(o$corrections %||% list()),
  loglik = .num1(rec$estimate$loglik), logposterior = .num1(rec$estimate$logposterior),
  loglik_method = .chr1(rec$estimate$loglik_method),
  loglik_laplace = .num1(rec$estimate$loglik_laplace),
  rescored = .num1(refs$rescored), prior_term = .num1(refs$prior_term),
  quad5 = .num1(refs$quad5$quadrature), quad5_gap = .num1(refs$quad5$gap),
  exact = .num1(refs$exact), exact_alt = .num1(refs$exact_alt),
  exact_kind = .chr1(refs$exact_kind),
  is_check_maxdiff = if (is.data.frame(refs$is_check)) max(abs(refs$is_check$is -
    refs$is_check$reference)) else NA_real_,
  units = .int1(refs$units$n), unit_dmax = .int1(refs$units$dmax),
  unit_lambda_min = .num1(refs$units$lambda_min), units_below_one = .int1(refs$units$below_one),
  lapcorr_status = .chr1(corr$status), lapcorr_applied = .lgl1(corr$applied),
  lapcorr_steps = .int1(corr$steps),
  nweak = .int1(rec$identifiability$nweak),
  warnings = length(warns), first_warning = if (length(warns)) substr(warns[1], 1, 200) else NA_character_,
  opcount_exp = .num1(rec$opcounts[["exp"]] %||% NA), opcount_frechet = .num1(rec$opcounts[["frechet"]] %||% NA),
  instrumented = all(rec$instrumented, na.rm = TRUE),
  warm_dx = if (length(rec$warm$raw) && length(rec$warm$raw) == length(est)) .fin(max(abs(rec$warm$raw - est))) else NA_real_,
  warm_dll = .num1(rec$estimate$logposterior) - .num1(rec$warm$logposterior %||% NA),
  runs_overshot = sum(vapply(runs, function(r) isTRUE(r$overshot), logical(1))),
  overshot_maxabs_min = .fin(suppressWarnings(min(vapply(runs, function(r)
    if (isTRUE(r$overshot)) .num1(r$maxabs_raw) else NA_real_, 0), na.rm = TRUE))),
  runs_maxabs_lt2 = sum(vapply(runs, function(r) isTRUE(.num1(r$maxabs_raw) < 2), logical(1))),
  maxabs_raw_max = .fin(suppressWarnings(max(vapply(runs, function(r) .num1(r$maxabs_raw), 0),
    na.rm = TRUE))),
  # Runs whose verdict probed: not the prior warm-up's (it runs with the probe
  # off), not one outside an optimiser stage, and none under `noprobe`.
  probe_runs = if (identical(V$optimcontrol$overshoot, "off")) 0L else
    sum(vapply(runs, function(r) { k <- .int1(r$stage)
      is.finite(k) && k >= 1L && k <= length(st) && !identical(st[[k]]$kind, "warmup")
    }, logical(1))),
  probe_secs = .num1(probe$secs %||% NA), probe_found = .lgl1(probe$found),
  lapcorr_continuation = .chr1(corr$continuation),
  lapcorr_secs = .num1(corr$seconds["total"]), lapcorr_gradients = .num1(corr$evaluations["gradients"]),
  cont_hessian = is.matrix(corr$hessian),
  cont_delta_se_max = .fin(suppressWarnings(max(abs(as.numeric(cse$delta_se)), na.rm = TRUE))),
  hess_secs_central = .num1(cse$secs_central %||% NA), hess_secs_forward = .num1(cse$secs_forward %||% NA),
  hess_secs_exact = .num1(cse$secs_exact %||% NA),
  hess_secs_laplace_x = .num1(cse$secs_laplace_x %||% NA),
  se_rel_exact = .num1(cse$rel["exact"]),
  se_rel_forward = .num1(cse$rel["forward"]), se_rel_laplace_x = .num1(cse$rel["laplace_x"]),
  se_rel_laplace_est = .num1(cse$rel["laplace_est"]), se_rel_reported = .num1(cse$rel["reported"]),
  stringsAsFactors = FALSE)
save_rec()
.stamp("DONE", CELL$id, rec$status)
