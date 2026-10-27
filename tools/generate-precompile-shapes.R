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
# One model per route and indicator kind, and the list should stay that short
# -- with one exception: on the Laplace route a Gaussian intercept is integrated
# by the subject filter, so a second Gaussian model keeps a Laplace coordinate.
# A model's type does not depend on its size, layout or expressions (see the
# header of inst/julia/ContinuousTimeSEM/src/precompile_workload.jl), so the
# code one captured fit compiles is the code every model of its route and kind
# runs: a second Gaussian model on a route compiles about a second's worth on
# dev2 (2026-10-13) where it used to compile a minute or more. What differs is
# whether an indicator is categorical -- a filter workspace type, so that
# Gaussian fits carry none of that code -- hence a binary and an ordinal model
# here, which share it. Every model costs build time and peak memory once per
# engine version and load time in every session (Julia inserting the image's
# method-table backedges), so another has to be measured against that before
# it goes in: review/PRECOMPILE-workload-2026-10-12.md and
# review/PRECOMPILE-general-types-2026-10-13.md have the method and the numbers.

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
indicators <- function(d, names) {
  for (nm in names) d[[nm]] <- d$eta + stats::rnorm(nrow(d), 0, 0.5)
  d$eta <- NULL
  d
}
# One latent, several indicators loading 1, a random continuous intercept.
cint_model <- function(names) {
  n <- length(names)
  mvar <- diag(0, n)
  for (i in seq_len(n)) mvar[i, i] <- "mvar"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = n, manifestNames = names, latentNames = "eta1",
    LAMBDA = matrix(1, n, 1), MANIFESTMEANS = matrix(0, n, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
    manifesttype = rep(0L, n))))
  m$pars$indvarying <- m$pars$param %in% "cint"
  m
}
# The same with categorical indicators: binary (1) or four-category ordinal (2).
categorical <- function(d, names, type) {
  for (nm in names) d[[nm]] <- if (type == 1L) {
    stats::rbinom(nrow(d), 1, stats::plogis(d$eta))
  } else as.integer(cut(d$eta + stats::rlogis(nrow(d)), c(-Inf, -1, 0.4, 1.9, Inf)))
  d$eta <- NULL
  d
}
categorical_model <- function(names, type, asymptotes = NULL) {
  n <- length(names)
  args <- list(type = "ct", n.latent = 1, n.manifest = n, manifestNames = names,
    latentNames = "eta1", LAMBDA = matrix(1, n, 1), CINT = matrix("cint"),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, n), manifesttype = rep(type, n))
  if (type == 1L) args$MANIFESTMEANS <- matrix(0, n, 1)
  if (type == 2L) args$ncategories <- rep(4L, n)
  if (!is.null(asymptotes)) args$ASYMPTOTES <- asymptotes
  m <- suppressWarnings(suppressMessages(do.call(ctModel, args)))
  m$pars$indvarying <- m$pars$param %in% "cint"
  m
}
# What a user does with a fit after making it. Each of these compiled 1-7 s on
# its first call in a fresh session (scripts/first-call-compile, 2026-10-16:
# sampling 3.6 s and a sampled ctFit 7.2 s, ctGenerateFromFit 3.3, a profile
# 2.5, ctLOO 1.9, the particle filter 1.3), so they are recorded after the
# Gaussian fits and replayed at build time too: a longer install rather than a
# wait at first use. Settings as small as still reach each route's code, since
# the replay runs them; the types are what is captured, not the results. The
# replay runs what was recorded, on one thread, so a step that is slow here is
# slow in every build: on these 8 subjects the cint sd sits near zero, and a
# profile of every parameter walked refits to the iteration cap (25 refits,
# 4.6 min) and the Laplace marginal sampler ran at full tree depth (2.2 min),
# dev1 2026-10-16. Hence one parameter, few iterations and a shallow tree:
# values, not types, so the same code is compiled.
followups <- list(
  kalman = function(fit, dat, model) ctKalman(fit),
  loo = function(fit, dat, model) ctLOO(fit, folds = 2, cores = 1L),
  psis = function(fit, dat, model) ctLOO(fit, method = "psis", nsamples = 50),
  profile = function(fit, dat, model) ctFitProfile(fit, parameters = 1L,
    points = 2L, maxiter = 20L, cores = 1L),
  generate = function(fit, dat, model) ctGenerateFromFit(fit, nsamples = 2, cores = 1L),
  particle = function(fit, dat, model) ctParticleLik(fit, particles = 20, cores = 1L),
  sample = function(fit, dat, model) ctFitUncertainty(fit, "sample",
    control = list(chains = 1L, warmup = 5L, draws = 5L, maxdepth = 3L,
      processes = FALSE), cores = 1L))
# A sampled fit is a route of its own (`intoverpop = 'auto'` resolves to the
# joint target when `optimize = FALSE`), the same whatever the shape's route,
# so it follows one shape only. That shape's Laplace term is exact, so the
# sampled fit and the `sample` follow-up both take its exact marginal; the
# joint is asked for by name (`intoverpop = 'none'`), where it goes through the
# map, the filter integrating the intercept (`ctsem_sample_map`). SAEM's
# kernel, the joint's sampler where TI predictor values are sampled, is asked
# for by name too, and so is 'checked' (its normal check, DAIS and importance
# sampler), which a default fit no longer runs.
sampled_fit <- list(sampled = function(fit, dat, model) ctFit(dat, model,
  backend = "julia", optimize = FALSE, priors = TRUE, cores = 1L, verbose = 0,
  sampleControl = list(iter = 10L, chains = 1L, warmup = 5L, maxdepth = 3L,
    processes = FALSE)),
  map = function(fit, dat, model) ctFitUncertainty(fit, "sample",
    control = list(chains = 1L, warmup = 5L, draws = 5L, maxdepth = 3L,
      processes = FALSE, intoverpop = "none"), cores = 1L),
  saem = function(fit, dat, model) ctFitUncertainty(fit, "sample",
    control = list(chains = 1L, warmup = 5L, draws = 5L, maxdepth = 3L,
      processes = FALSE, intoverpop = "none", sampler = "saem"), cores = 1L),
  checked = function(fit, dat, model) ctFitUncertainty(fit, "checked", cores = 1L))

# What the quadrature correction runs once its shift check flags units: the
# continuation, its screen, its rounds and the continuation's Hessian. A
# shape whose own fit's check left it alone records them all the same, by
# running the continuation again without the check, whose screen then sets
# the rounds off as it did before the check came first. Where the fit's
# correction moved the estimate, the fit has recorded them already.
correction_rounds <- function(fit, dat, model) {
  if (isTRUE(fit$laplace$correction$applied)) return(invisible(NULL))
  ns$.ctLaplaceContinue(fit, control = utils::modifyList(
    ns$.ctLaplaceContinueDefaults, list(check = FALSE)))
}

shapes <- list(
  # One latent, two Gaussian indicators, a random intercept, on each route:
  # `image_model()` in tests/testthat/helper-julia.R, which tests fit so that a
  # fit costs seconds.
  gaussian_augmented = list(route = "augmented", model = function() cint_model(c("y1", "y2")),
    data = function() indicators(sim_latent(8, 6), c("y1", "y2")),
    followups = followups),
  gaussian_laplace = list(route = "laplace", model = function() cint_model(c("y1", "y2")),
    data = function() indicators(sim_latent(8, 6), c("y1", "y2")),
    followups = c(followups, sampled_fit)),
  # The same without the random effect: the plain filter, whose Hessians (the
  # Newton finish's and the certification's) take the hand-coded second order
  # (`tangent_hessian.jl`), which the shapes with random effects do not reach.
  # Its follow-ups are the augmented shape's code, so it has none.
  gaussian_fixed = list(route = "augmented", model = function() {
      m <- cint_model(c("y1", "y2"))
      m$pars$indvarying <- FALSE
      m
    },
    data = function() indicators(sim_latent(8, 6), c("y1", "y2"))),
  # The intercept above is integrated by each subject's filter, leaving the
  # Laplace units no coordinate. A random drift beside it and y2's intercept
  # stays a coordinate, and the filter integrates the two intercepts given it
  # (a full-rank block, `CTSEMHybrid`; one alone beside the drift would stay a
  # coordinate, by the cost rule): the inner solve, its probes and the
  # logdet's derivatives
  # at the Gaussian workspace type, which every Gaussian Laplace fit with a
  # nonlinear effect runs, and the quadrature correction's rounds, which the
  # drift effect sets off on these ten subjects. Its fit stops at the Hessian's
  # normal, without importance sampling. The correction was off here while
  # every call was recorded (2,100 of 2,500 were the continuation's
  # evaluations); at two calls a signature it is about 40, and without it a
  # session's first such fit compiled the rounds. The shift check comes first
  # in the correction, and the rounds follow it whatever it finds
  # (`correction_rounds`).
  gaussian_laplace_block = list(route = "laplace", model = function() {
      mvar <- diag(0, 2)
      diag(mvar) <- "mvar"
      m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
        n.manifest = 2, manifestNames = c("y1", "y2"), latentNames = "eta1",
        LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(c(0, "mm2"), 2, 1),
        CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
        manifesttype = c(0L, 0L))))
      m$pars$indvarying <- m$pars$param %in% c("cint", "mm2", "drift_eta1")
      m
    },
    data = function() indicators(sim_latent(10, 6), c("y1", "y2")),
    optimcontrol = list(uncertainty = "hessian"),
    followups = list(rounds = correction_rounds)),
  # A model with a categorical indicator is a different filter workspace type
  # -- the categorical update is gated on it, so Gaussian fits do not carry
  # that code -- so each route needs one. Binary and ordinal share it.
  binary_augmented = list(route = "augmented",
    model = function() categorical_model(c("b1", "b2"), 1L),
    data = function() categorical(sim_latent(8, 6), c("b1", "b2"), 1L)),
  ordinal_laplace = list(route = "laplace",
    model = function() categorical_model(c("o1", "o2"), 2L),
    data = function() categorical(sim_latent(10, 6), c("o1", "o2"), 2L)),
  # A guessing asymptote on one binary item, on enough subjects that the
  # quadrature correction runs its rounds and forms the continuation's Hessian
  # at the categorical workspace type, which without it a session's first such
  # fit compiled; and the asymptote's own mixture moments
  # (`_asymptote_partials!`), which the binary and ordinal shapes do not reach.
  # Recorded while a failed unit sweep still fell back to forward mode over
  # every unit, it set that off once and took the build from 11.8 to 17.0 GB
  # (dev2); without that fallback it costs about 50 s of build and no memory.
  asymptote_laplace = list(route = "laplace",
    model = function() categorical_model(c("b1", "b2"), 1L,
      asymptotes = c("guess", NA)),
    data = function() categorical(sim_latent(40, 12, cints = stats::rnorm(40, 0, 1)),
      c("b1", "b2"), 1L),
    followups = list(rounds = correction_rounds))
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
  "end",
  # The batch evaluations take their draws as columns, whichever dimension
  # is the larger: the uncertainty stage sends hundreds of them, each with its
  # draws, and written whole they made the file thirty times its size.
  "function ctsem_generator_draw_columns(x)",
  "  x isa AbstractMatrix && size(x, 2) > 3 && (x = x[:, 1:3])",
  "  ctsem_generator_literal(x)",
  "end",
  # An argument's type, for a call's signature (`signature()` below).
  "ctsem_generator_type(x) = string(typeof(x))", sep = "\n")))
draws <- c("ctsem_parameter_matrices")
drawcols <- c("ctsem_evaluate_batch", "ctsem_evaluate_batch_gradient",
  "ctsem_coreset_batch", "ctsem_coreset_batch_gradient", "ctsem_coreset_terms")

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
  writer <- if (fn %in% draws) "ctsem_generator_draws" else
    if (fn %in% drawcols) "ctsem_generator_draw_columns" else "ctsem_generator_literal"
  if (!inherits(x, "JuliaProxy")) return(raw_call(writer, x))
  for (p in rec$proxies) if (identical(p$proxy, x)) {
    return(if (is.null(p$code)) raw_call(writer, x) else p$code)
  }
  stop("an argument is a Julia object no recorded call produced")
}
# A call's signature: the engine function, its keywords' names, and every
# argument's Julia type -- what decides which code a replay compiles. A
# second call of one signature compiles nothing the first did not, so at most
# `repeats` of each are kept: the uncertainty stage makes hundreds of batch
# evaluations of one signature a fit (DAIS's iterations), and every one of them
# was replayed in every build. Calls that return an object a later call may
# refer to (a proxy) are always kept, and so are calls that change the
# engine's state (a name ending in `!`), whose effect later calls may need.
repeats <- 2L
signature <- function(fn, args) {
  nm <- names(args); if (is.null(nm)) nm <- rep("", length(args))
  types <- vapply(args, function(x) as.character(raw_call("ctsem_generator_type", x)), "")
  paste(fn, paste(nm, types, sep = ":", collapse = ","), sep = "|")
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
  sig <- signature(fn, args)
  pos <- lapply(args[!nzchar(nm)], encode, fn = fn)
  kw <- lapply(args[nzchar(nm)], encode, fn = fn)
  out <- raw_call(name, ..., .defer = .defer)
  seen <- if (is.null(rec$seen[[sig]])) 0L else rec$seen[[sig]]
  rec$seen[[sig]] <- seen + 1L
  if (seen >= repeats && !inherits(out, "JuliaProxy") && !endsWith(fn, "!")) {
    rec$skipped <- rec$skipped + 1L
    return(out)
  }
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
  rec$calls <- list(); rec$proxies <- list(); rec$seen <- list(); rec$skipped <- 0L
  rec$on <- TRUE
  # The starting values a user's `inits = NULL` fit draws, from a fixed seed.
  set.seed(1)
  args <- list(dat, model, backend = "julia", intoverpop = s$route, cores = 1L)
  if (!is.null(s$optimcontrol)) args$optimcontrol <- s$optimcontrol
  fit <- tryCatch(suppressWarnings(suppressMessages(do.call(ctFit, args))),
    finally = rec$on <- FALSE)
  # A replay that could not build its model would compile nothing and say so
  # only in a stopwatch.
  if (!any(vapply(rec$calls, function(x) x$fn == "ekf_from_columns", TRUE))) {
    stop(name, ": the fit made no ekf_from_columns call")
  }
  cat(sprintf("%-22s %3d engine calls (%d repeats left out), loglik %.3f\n", name,
    length(rec$calls), rec$skipped, fit$estimate$loglik))
  # Recorded into the same replay, after the fit whose objects they use. One
  # that cannot be recorded -- an R callback handed to Julia cannot be replayed
  # at build time -- is said so and left out; the calls it made before that
  # still replay.
  for (f in names(s$followups)) {
    before <- length(rec$calls)
    rec$on <- TRUE
    set.seed(1)
    tryCatch(suppressWarnings(suppressMessages(s$followups[[f]](fit, dat, model))),
      error = function(e) cat(sprintf("  %s: %s stopped recording: %s\n", name, f,
        conditionMessage(e))),
      finally = rec$on <- FALSE)
    cat(sprintf("  %-20s %3d engine calls\n", f, length(rec$calls) - before))
  }
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
