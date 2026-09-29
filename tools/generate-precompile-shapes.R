# Regenerate what the Julia engine precompiles: for each model below, every
# engine call a default ctFit makes, captured from a real fit through the bridge.
#
#     Rscript tools/generate-precompile-shapes.R
#
# It rewrites inst/julia/ContinuousTimeSEM/src/precompile_shapes.jl, which is
# checked in -- the engine must be precompilable without R present.
#
# Why captured and not written. The engine compiles code specialised to a
# model's *type* -- its matrix dimensions, its transform templates -- and to
# the exact argument types of each call, keywords included. A workload that
# builds a plausible model and calls a plausible subset of the engine compiles
# code no fit asks for: the first version of this file hand-wrote a model that
# was wrong in four ways at once, and the next captured the model but not the
# fit, so it covered a value, a gradient and two optimiser iterations, while
# the Newton finish, the certification Hessian, the smoother and the Laplace
# continuation all compiled inside the user's first fit. Replaying what R
# actually sent covers the whole default pipeline, and stays right when R
# changes what it sends -- provided this script is rerun. The test in
# test-julia-precompile.R is what notices when it has not been.
#
# Each model costs build time (roughly its own first-fit compile, once per
# engine version) and, less obviously, load time in *every* session: the image
# holds its specialisations whether or not the session fits that shape. So the
# list is short, and a model earns a place by being a likely first fit.

root <- file.path(dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(FALSE), value = TRUE)[1])), "..")
suppressMessages(devtools::load_all(root, compile = FALSE, quiet = TRUE))
suppressMessages(ctJuliaSetup(threads = 1L))
ns <- asNamespace("ctsem")

# ---- the models ---------------------------------------------------------------
# Small data on purpose: the replay runs inside the package build, and it is the
# types that are captured, not the estimates. Enough subjects that a random
# effect has a variance to estimate.

sim_latent <- function(nsub, nobs, drift = -0.3, cints = stats::rnorm(nsub, 0, 0.5)) {
  a <- exp(drift); b <- (a - 1) / drift; q <- 0.64 * (a^2 - 1) / (2 * drift)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    eta <- numeric(nobs); eta[1] <- stats::rnorm(1)
    for (t in 2:nobs) eta[t] <- a * eta[t - 1] + b * cints[i] + stats::rnorm(1, 0, sqrt(q))
    data.frame(id = i, time = seq_len(nobs) - 1, eta = eta)
  }))
}
indicators <- function(d, names, binary) {
  for (nm in names) d[[nm]] <- if (binary)
    stats::rbinom(nrow(d), 1, stats::plogis(d$eta)) else d$eta + stats::rnorm(nrow(d), 0, 0.5)
  d$eta <- NULL
  d
}
# One latent, several indicators loading 1, a random continuous intercept.
cint_model <- function(names, type) {
  n <- length(names)
  mvar <- diag(0, n)
  if (type == 0L) for (i in seq_len(n)) mvar[i, i] <- "mvar"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = n, manifestNames = names, latentNames = "eta1",
    LAMBDA = matrix(1, n, 1), MANIFESTMEANS = matrix(0, n, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
    manifesttype = rep(type, n))))
  m$pars$indvarying <- m$pars$param %in% "cint"
  m
}
# A random drift `-log1p_exp(-raw)`, raw ~ N(1, 1), 40 subjects: the bench's
# acnonlin data (dev/optimbench/cells.R, ac_data, seed 3) and gated-gaps config
# A1 (gg_genA, seed 1), copied with their seeds.
drift_data <- function() {
  set.seed(3L)
  drift <- -log1p(exp(-stats::rnorm(40L, 1, 1)))
  do.call(rbind, lapply(seq_len(40L), function(i) {
    a <- drift[i]; decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(8L); latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(7L)) latent[t + 1L] <- decay * latent[t] + stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(8L) - 1L, Y1 = latent + stats::rnorm(8L, 0, 0.3))
  }))
}
drift_cint_data <- function(nsub = 40L, ntimes = 6L) {
  set.seed(1L)
  baseline <- stats::rnorm(nsub, 2, 2)
  start <- stats::rnorm(nsub, baseline / 2, 1)
  raw <- stats::rnorm(nsub, 1 + (baseline - 2) / 2, 1)
  drift <- -log1p(exp(-raw))
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    a <- drift[i]; decay <- exp(a)
    intercept <- (baseline[i] / a) * (decay - 1)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- start[i]
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] + intercept +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L, Y1 = latent + stats::rnorm(ntimes, 0, 0.5))
  }))
}
# One latent, one indicator, a random drift kept negative by its transform:
# with a fixed zero intercept and initial mean, or a free intercept and the
# default (random) initial mean.
drift_model <- function(cint) {
  if (identical(cint, 0)) {
    m <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
      MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
      DRIFT = "drift|-log1p_exp(-param)|TRUE"))
    m$pars$indvarying <- m$pars$param %in% "drift"
    return(m)
  }
  suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = cint,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), DRIFT = "drift|-log1p_exp(-param)|TRUE"))
}

one_latent_model <- function() suppressMessages(ctModel(silent = TRUE, type = "ct",
  LAMBDA = matrix(1), manifestNames = "Y1", latentNames = "eta1"))
one_latent_data <- function() {
  d <- sim_latent(16, 8, cints = rep(0, 16))
  data.frame(id = d$id, time = d$time, Y1 = d$eta + stats::rnorm(nrow(d), 1, 0.5))
}

shapes <- list(
  # One latent, two Gaussian indicators, a random intercept.
  gaussian_augmented = list(route = "augmented", model = function() cint_model(c("y1", "y2"), 0L),
    data = function() indicators(sim_latent(8, 6), c("y1", "y2"), FALSE)),
  gaussian_laplace = list(route = "laplace", model = function() cint_model(c("y1", "y2"), 0L),
    data = function() indicators(sim_latent(8, 6), c("y1", "y2"), FALSE)),
  # The same with three binary indicators.
  binary_augmented = list(route = "augmented", model = function() cint_model(c("b1", "b2", "b3"), 1L),
    data = function() indicators(sim_latent(8, 6), c("b1", "b2", "b3"), TRUE)),
  binary_laplace = list(route = "laplace", model = function() cint_model(c("b1", "b2", "b3"), 1L),
    data = function() indicators(sim_latent(8, 6), c("b1", "b2", "b3"), TRUE)),
  # Two latents measured one each, random intercepts on both.
  two_latent_augmented = list(route = "augmented", model = function() {
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
      manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"), LAMBDA = diag(2),
      MANIFESTMEANS = matrix(0, 2, 1), CINT = matrix(c("cint1", "cint2")),
      MANIFESTVAR = matrix(c("mv1", 0, 0, "mv2"), 2))))
    m$pars$indvarying <- m$pars$param %in% c("cint1", "cint2")
    m
  }, data = function() {
    a <- sim_latent(8, 6); b <- sim_latent(8, 6)
    data.frame(id = a$id, time = a$time, Y1 = a$eta + stats::rnorm(nrow(a), 0, 0.5),
      Y2 = b$eta + stats::rnorm(nrow(b), 0, 0.5))
  }),
  # One latent, one indicator, individual differences in the drift. Here the
  # data are not small, because a replay covers the paths its own fit took:
  # a drift weakly informed by each subject's few observations is what sends
  # the inner solve down the gated floor and the continuation through its
  # rounds and, where it converges, its Hessian. Fitted on easy data none of
  # that ran, and it compiled in the first real fit instead -- 61 s of it on
  # the bench's gA1.
  drift_laplace = list(route = "laplace", model = function() drift_model(0),
    data = drift_data),
  # The same with a free intercept and the default random initial mean.
  drift_cint_laplace = list(route = "laplace", model = function() drift_model("cint"),
    data = drift_cint_data),
  # One latent and one indicator with the writer's defaults: the filter type
  # the test suite compiled most often (44 files) and three vignettes compile,
  # which none of the models above has. Without random effects on the augmented
  # route (its Hessian), and with the default random T0MEANS and MANIFESTMEANS
  # on the Laplace route, whose filter is the same type.
  one_latent_fixed = list(route = "augmented", model = function() {
    m <- one_latent_model(); m$pars$indvarying <- FALSE; m }, data = one_latent_data),
  one_latent_laplace = list(route = "laplace", model = one_latent_model,
    data = one_latent_data)
)
only <- Sys.getenv("CTSEM_PRECOMPILE_SHAPES", "")
if (nzchar(only)) shapes <- shapes[strsplit(only, ",", fixed = TRUE)[[1]]]

# ---- the recorder ---------------------------------------------------------------
# Every call goes through `.ctJuliaCall()` (R/ctJuliaBridge.R), so wrapping it
# sees the whole conversation. An argument is recorded as what Julia received:
# a value R sent is written as Julia's own `repr` of it, which round-trips the
# element type (Int64 against Float64, a Vector against a scalar); a proxy an
# earlier recorded call returned is written as a reference to that call.

raw_call <- get(".ctJuliaCall", envir = ns)
invisible(ctsem:::.ctJuliaEval(paste(
  "function ctsem_generator_literal(x)",
  "  x isa AbstractArray && length(x) > 50000 && error(\"literal too large: \", summary(x))",
  "  repr(x)",
  "end",
  # A matrix of draws is replayed with three of them, along whichever dimension
  # holds the draws: its type is what matters, and a thousand draws would be a
  # megabyte of source. Only for the calls that take draws -- the data matrix
  # must keep every column its times and timesteps describe.
  "function ctsem_generator_draws(x)",
  "  if x isa AbstractMatrix && length(x) > 200",
  "    x = size(x, 1) >= size(x, 2) ? x[1:3, :] : x[:, 1:3]",
  "  end",
  "  ctsem_generator_literal(x)",
  "end", sep = "\n")))
draws <- c("ctsem_parameter_matrices", "ctsem_laplace_population")

rec <- new.env()
skip <- c("ctsem_set_max_chunks!", "ctsem_max_chunks", "ctsem_set_interrupt!",
  "ctsem_tune_bridge!")
engine_name <- function(name) {
  if (!grepl("^ContinuousTimeSEM\\.:?", name)) return(NULL)
  fn <- sub("^ContinuousTimeSEM\\.:?", "", name)
  if (fn %in% skip) NULL else fn
}
# Written when the value is used, because only the call that uses it says
# whether it is a matrix of draws.
encode <- function(x, fn) {
  if (is.function(x)) stop("an R callback cannot be replayed at build time")
  writer <- if (fn %in% draws) "ctsem_generator_draws" else "ctsem_generator_literal"
  if (!inherits(x, "JuliaProxy")) return(raw_call(writer, x))
  for (p in rec$proxies) if (identical(p$proxy, x)) {
    return(if (is.null(p$code)) raw_call(writer, x) else p$code)
  }
  stop("an argument is a Julia object no recorded call produced")
}
recorder <- function(name, ..., .defer = FALSE) {
  if (!isTRUE(rec$on)) return(raw_call(name, ..., .defer = .defer))
  if (identical(name, "RConnector.EnforcedProxy")) {
    out <- raw_call(name, ..., .defer = .defer)
    rec$proxies[[length(rec$proxies) + 1L]] <- list(proxy = out, code = NULL)
    return(out)
  }
  fn <- engine_name(name)
  if (is.null(fn)) return(raw_call(name, ..., .defer = .defer))
  args <- list(...)
  nm <- names(args); if (is.null(nm)) nm <- rep("", length(args))
  pos <- lapply(args[!nzchar(nm)], encode, fn = fn)
  kw <- lapply(args[nzchar(nm)], encode, fn = fn)
  out <- raw_call(name, ..., .defer = .defer)
  k <- length(rec$calls) + 1L
  rec$calls[[k]] <- list(fn = fn, pos = pos, kw = kw)
  if (inherits(out, "JuliaProxy")) {
    rec$proxies[[length(rec$proxies) + 1L]] <- list(proxy = out,
      code = sprintf("_PrecompileRef(%d)", k))
  }
  out
}
unlockBinding(".ctJuliaCall", ns)
assign(".ctJuliaCall", recorder, envir = ns)
lockBinding(".ctJuliaCall", ns)

emit_call <- function(call) {
  pos <- paste(unlist(call$pos), collapse = ", ")
  kw <- if (length(call$kw)) paste0(":", names(call$kw), " => ",
    unlist(call$kw), collapse = ", ") else ""
  sprintf("        (:%s, Any[%s], Pair{Symbol,Any}[%s]),", call$fn, pos, kw)
}

body <- character()
for (name in names(shapes)) {
  s <- shapes[[name]]
  set.seed(1)
  dat <- s$data()
  model <- s$model()
  # A later model must not find an earlier one's objective in the session cache.
  cache <- get(".ct_julia_cache", envir = ns)
  cache$objectives <- new.env(parent = emptyenv())
  rec$calls <- list(); rec$proxies <- list(); rec$on <- TRUE
  # The starting values a user's `inits = NULL` fit draws, from a fixed seed.
  set.seed(1)
  fit <- tryCatch(suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "julia", intoverpop = s$route, cores = 1L))),
    finally = rec$on <- FALSE)
  # A replay that could not build its model would compile nothing and say so
  # only in a stopwatch.
  if (!any(vapply(rec$calls, function(x) x$fn == "ekf_from_columns", TRUE))) {
    stop(name, ": the fit made no ekf_from_columns call")
  }
  cat(sprintf("%-22s %3d engine calls, loglik %.3f\n", name, length(rec$calls),
    fit$estimate$loglik))
  body <- c(body, sprintf("    %s = Any[", name),
    vapply(rec$calls, emit_call, ""), "    ],")
}

out <- c(
  "# GENERATED by tools/generate-precompile-shapes.R -- do not edit by hand.",
  "#",
  "# For each model listed in that script, every engine call a default ctFit",
  "# made, in order, with the arguments Julia received: values as their own",
  "# `repr`, and `_PrecompileRef(k)` for the object the k-th call returned. The",
  "# workload in precompile_workload.jl replays them at build time. Regenerate",
  "# after any change to what R sends; a stale file costs nothing but the",
  "# optimisation -- first fits compile again, silently. test-julia-precompile.R",
  "# is what catches that.",
  "",
  "const _PRECOMPILE_SHAPES = (",
  body,
  ")")
path <- file.path(root, "inst", "julia", "ContinuousTimeSEM", "src", "precompile_shapes.jl")
# Binary, so Windows writes the LF endings the rest of the engine uses.
con <- file(path, "wb")
writeLines(out, con, sep = "\n")
close(con)
cat("wrote", path, "--", sum(nchar(out)), "bytes\n")
