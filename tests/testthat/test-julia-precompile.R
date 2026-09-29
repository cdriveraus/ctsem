# The engine compiles a whole default fit for a few model shapes at build time,
# by replaying the engine calls a real fit made (tools/generate-precompile-
# shapes.R writes them). That only helps while two things stay true, and
# nothing about either failing is visible -- the package loads, the fit runs,
# and it is merely slow again:
#
#   * the models it captured still have the types ctsem's model writer emits,
#     and the calls it replays are still the calls a fit makes;
#   * the package image is still usable in a session started by the R bridge.
#     For months it was not: JuliaConnectoR loads Pkg before the engine, and
#     that alone failed the verification of ~4800 of the image's methods, so
#     every R session compiled the precompiled code again.
#
# The first test asks the cheap question about the first point. The second
# asks the real one about both: a fresh session's first fit of a captured model
# must compile almost nothing.
#
# When either fails, regenerate:  Rscript tools/generate-precompile-shapes.R

# The generator's `gaussian_augmented` / `gaussian_laplace` model: one latent,
# two Gaussian indicators, a random intercept. Must stay identical to it.
.precompile_model <- function() {
  mvar <- diag(0, 2); diag(mvar) <- "mvar"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("y1", "y2"), latentNames = "eta1",
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
    manifesttype = c(0L, 0L))))
  m$pars$indvarying <- m$pars$param %in% "cint"
  m
}
.precompile_data <- function(nsub = 20L, nobs = 8L) {
  set.seed(4)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    eta <- cumsum(stats::rnorm(nobs, 0, 0.5)) + stats::rnorm(1)
    data.frame(id = i, time = seq_len(nobs) - 1,
      y1 = eta + stats::rnorm(nobs, 0, 0.5), y2 = eta + stats::rnorm(nobs, 0, 0.5))
  }))
}

test_that("the precompiled models still match what ctsem's model writer emits", {
  skip_without_julia()
  module <- ctsem:::.ctJuliaModule(NULL)
  skip_if(is.null(module$ctsem_shape_is_precompiled),
    "engine predates ctsem_shape_is_precompiled")
  nz <- function(x, empty) { x[is.na(x)] <- empty; x }
  v <- ctsem:::.ctJuliaVector

  # The captured model, and the same with its second loading free: a default
  # template the captured model does not use, which must not make it another
  # type, because every spec carries every default template's group.
  free_loading <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("y1", "y2"), latentNames = "eta1",
    LAMBDA = matrix(c(1, "lam2"), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0),
    MANIFESTVAR = matrix(c("mvar", 0, 0, "mvar"), 2), manifesttype = c(0L, 0L))))
  free_loading$pars$indvarying <- free_loading$pars$param %in% "cint"
  models <- list(captured = .precompile_model(), free_loading = free_loading)
  for (route in c("augmented", "laplace")) for (nm in names(models)) {
    spec <- suppressWarnings(suppressMessages(ctFit(.precompile_data(3L, 4L),
      models[[nm]], backend = "julia", fit = FALSE, intoverpop = route)))
    table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
    matched <- ctsem:::.ctBackendJuliaValue(module$ctsem_shape_is_precompiled(
      v(as.character(table$matrix)), v(as.integer(table$row)), v(as.integer(table$col)),
      v(nz(as.integer(table$parnumber), 0L)), v(nz(as.numeric(table$value), NaN)),
      v(nz(as.character(table$transform), "")),
      v(nz(as.character(table$predicttransform), "")),
      v(nz(as.character(table$updatetransform), "")),
      v(nz(as.character(table$tdtransform), ""))))
    expect_true(matched, label = sprintf("the %s gaussian model (%s) is a precompiled shape",
      route, nm))
  }
})

test_that("every template a default model writes is one of the engine's default templates", {
  skip_without_julia()
  # Every spec carries a group for each default template, so that which of them
  # a model uses does not change its type. A default the writer starts emitting
  # and the engine's list lacks would quietly bring that back: a model freeing
  # one more cell compiling the whole pipeline again. Fit-free: parameter
  # tables only, and one engine call.
  module <- ctsem:::.ctJuliaModule(NULL)
  skip_if(is.null(module$ctsem_nondefault_templates),
    "engine predates the default template groups")
  q <- function(e) suppressWarnings(suppressMessages(e))
  sim <- function(man, type = rep(0L, length(man)), td = FALSE, ti = FALSE) {
    set.seed(2)
    d <- data.frame(id = rep(1:6, each = 5), time = rep(0:4, 6))
    for (j in seq_along(man)) d[[man[j]]] <- switch(as.character(type[j]),
      "0" = stats::rnorm(30), "1" = stats::rbinom(30, 1, 0.5), "2" = sample(1:3, 30, TRUE))
    if (td) d$td1 <- stats::rbinom(30, 1, 0.3)
    if (ti) d$ti1 <- rep(stats::rnorm(6), each = 5)
    d
  }
  cases <- list(
    list(q(ctModel(type = "ct", n.latent = 2, n.manifest = 3, manifestNames = c("y1", "y2", "y3"),
      LAMBDA = matrix(c(1, "auto", "auto", 0, 0, 1), 3, 2))), sim(c("y1", "y2", "y3")), "augmented"),
    list(q(ctModel(type = "dt", n.latent = 1, n.manifest = 2, manifestNames = c("y1", "y2"),
      LAMBDA = matrix(c(1, "auto"), 2, 1))), sim(c("y1", "y2")), "laplace"),
    list(q(ctModel(type = "ct", n.latent = 1, n.manifest = 2, manifestNames = c("o1", "o2"),
      LAMBDA = matrix(c(1, "auto"), 2, 1), manifesttype = c(2L, 2L), ncategories = c(3L, 3L))),
      sim(c("o1", "o2"), c(2L, 2L)), "laplace"),
    list(q(ctModel(type = "ct", n.latent = 1, n.manifest = 1, manifestNames = "y1",
      LAMBDA = matrix(1), TDpredNames = "td1", TIpredNames = "ti1")),
      sim("y1", td = TRUE, ti = TRUE), "augmented"))
  transforms <- unlist(lapply(cases, function(cs) {
    spec <- q(ctFit(cs[[2]], cs[[1]], backend = "julia", fit = FALSE, intoverpop = cs[[3]]))
    tb <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
    t <- tb$transform[!is.na(tb$parnumber) & tb$parnumber > 0]
    t[!is.na(t) & nzchar(t)]
  }))
  expect_gt(length(transforms), 20)
  missing <- as.character(ctsem:::.ctBackendJuliaValue(
    module$ctsem_nondefault_templates(ctsem:::.ctJuliaVector(transforms))))
  expect_identical(missing, "", label = "templates default models write outside the default groups")
})

test_that("every captured model replayed when the engine image was built", {
  skip_without_julia()
  # A replay that throws is caught, so that a broken one cannot stop the
  # package loading: the build prints a warning and the image simply lacks that
  # model. This is where that becomes a failure.
  enabled <- tryCatch(isTRUE(ctsem:::.ctJuliaEval(
    "ContinuousTimeSEM._PRECOMPILE_WORKLOAD_ENABLED")), error = function(e) NA)
  skip_if(is.na(enabled), "engine predates the replayed workload")
  skip_if(!enabled, "engine image built with CTSEM_PRECOMPILE_WORKLOAD=false")
  missing <- as.character(ctsem:::.ctJuliaEval(paste0("join(string.(setdiff(",
    "collect(keys(ContinuousTimeSEM._PRECOMPILE_SHAPES)), ",
    "ContinuousTimeSEM._PRECOMPILE_REPLAYED)), \", \")")))
  expect_identical(missing, "", label = "captured models whose replay failed at build time")
})

test_that("a fresh session's first fit of a precompiled model compiles almost nothing", {
  skip_without_julia()
  # Its own R process, and so its own Julia: this session's Julia has compiled
  # whatever earlier tests ran, and a first fit is a property of a fresh one.
  # Most of its time is starting R, ctsem and Julia -- nothing cheaper sees the
  # bridge's load order, which is what broke this before.
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  from_source <- file.exists(file.path(root, "DESCRIPTION")) &&
    dir.exists(file.path(root, "R"))
  loader <- if (from_source) {
    sprintf("suppressMessages(pkgload::load_all(%s, compile = FALSE, quiet = TRUE))",
      deparse(root))
  } else "suppressMessages(library(ctsem))"
  trace <- normalizePath(tempfile("firstfit", fileext = ".jl"), winslash = "/",
    mustWork = FALSE)
  out <- tempfile("firstfit", fileext = ".rds")
  script <- tempfile("firstfit", fileext = ".R")
  # JuliaConnectoR pastes these into a command line, hence the quoting.
  opts <- trimws(paste(Sys.getenv("JULIACONNECTOR_JULIAOPTS"),
    paste0("--trace-compile=", shQuote(trace)), "--trace-compile-timing"))
  writeLines(c(
    sprintf("Sys.setenv(JULIACONNECTOR_JULIAOPTS = %s)", deparse(opts)),
    loader,
    paste(".precompile_model <-", paste(deparse(.precompile_model), collapse = "\n")),
    paste(".precompile_data <-", paste(deparse(.precompile_data), collapse = "\n")),
    "suppressMessages(ctJuliaSetup(threads = 1L))",
    sprintf("before <- length(readLines(%s, warn = FALSE))", deparse(trace)),
    "fit <- suppressWarnings(suppressMessages(ctFit(.precompile_data(),",
    "  .precompile_model(), backend = 'julia', intoverpop = 'augmented', cores = 1L)))",
    sprintf("saveRDS(list(before = before, lines = readLines(%s, warn = FALSE),",
      deparse(trace)),
    sprintf("  loglik = fit$estimate$loglik), %s)", deparse(out))), script)
  # The child inherits the environment, so it finds this session's libraries
  # (an installed ctsem under R CMD check) and its Julia.
  log <- withr::with_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
      stdout = TRUE, stderr = TRUE)))
  finished <- file.exists(out)
  expect_true(finished, label = paste(c("the fresh session finished", tail(log, 15)),
    collapse = "\n"))
  skip_if(!finished)
  res <- readRDS(out)
  expect_true(is.finite(res$loglik))
  compiled <- res$lines[-seq_len(res$before)]
  # What the engine itself compiled: its own methods and closures, and keyword
  # calls into them. The bridge compiles a little for new result types, which
  # is not what this is about.
  body <- sub("^#=.*=# ", "", compiled)
  engine <- grepl("^precompile\\(Tuple\\{(typeof\\()?ContinuousTimeSEM\\.", body) |
    (startsWith(body, "precompile(Tuple{typeof(Core.kwcall)") &
      grepl(", typeof(ContinuousTimeSEM.", body, fixed = TRUE))
  ms <- suppressWarnings(as.numeric(sub("^#=\\s*([0-9.]+) ms.*", "\\1", compiled)))
  # A healthy replay leaves a handful; a rotted one, or an image the session
  # cannot use, leaves hundreds and most of a minute. The bound tells those
  # apart and is not for tuning.
  expect_lt(sum(engine), 60, label = sprintf(
    "engine methods compiled in the first fit (%d, %.1f s of compile in all)",
    sum(engine), sum(ms, na.rm = TRUE) / 1000))
})

# The engine's numbers from its package image, and from the same source built
# without the workload -- so compiled in the session -- at the points the
# workload was captured at: value, gradient and exact Hessian of every captured
# model. The image carries compiled code and module state from the process
# that built it, and nothing else checks that they compute what the source
# says (review/ENGINE-image-state-2026-09-29.md, F6). `mutate` edits the
# second copy's source, which is how the check was shown to fail when the two
# differ; the test itself never passes one.
.image_vs_session <- function(mutate = NULL) {
  julia <- file.path(ctsem:::.ctJuliaBin(), ctsem:::.ctJuliaExeName())
  script <- normalizePath(testthat::test_path("image-vs-session.jl"), winslash = "/")
  image_env <- ctsem:::.ctJuliaEnvDir()
  # R's LD_LIBRARY_PATH (R's own and the system's libraries) is unset, as
  # JuliaConnectoR does for its Julia: with it, Julia's precompile workers load
  # the wrong libraries and die with a segmentation fault (Linux).
  run <- function(project, out, env = character()) {
    withr::with_envvar(c(LD_LIBRARY_PATH = NA, env), suppressWarnings(system2(julia, c("--startup-file=no",
      shQuote(script), shQuote(normalizePath(project, winslash = "/")), shQuote(out)),
      stdout = TRUE, stderr = TRUE)))
    if (file.exists(out)) readLines(out, warn = FALSE) else character()
  }
  from_image <- run(image_env, tempfile("image", fileext = ".txt"))
  # A copy at another path is another cache entry, and with its own depot
  # first it neither takes one of the user's image slots nor is found by a
  # later session of this source. The session's own depots follow, where the
  # packages and the dependencies' images are: an empty entry would add only
  # Julia's bundled depots, not the user's, and nothing would load.
  copy <- file.path(tempfile("ctsem-session-engine"), "engine")
  dir.create(copy, recursive = TRUE)
  file.copy(list.files(ctsem:::.ctJuliaEnginePath(), full.names = TRUE), copy,
    recursive = TRUE)
  if (!is.null(mutate)) mutate(copy)
  depot <- tempfile("ctsem-depot")
  dir.create(depot)
  depots <- as.character(ctsem:::.ctJuliaEval(sprintf("join(DEPOT_PATH, %s)",
    deparse(.Platform$path.sep))))
  from_session <- run(copy, tempfile("session", fileext = ".txt"),
    c(CTSEM_PRECOMPILE_WORKLOAD = "false", JULIA_DEPOT_PATH = paste(
      normalizePath(depot, winslash = "/"), depots, sep = .Platform$path.sep)))
  list(image = from_image, session = from_session)
}
.image_vs_session_diff <- function(lines) {
  parse <- function(x) {
    x <- x[!startsWith(x, "workload") & x != "done"]
    keys <- sub("^(\\S+ \\S+ \\S+) .*$", "\\1", x)
    stats::setNames(lapply(sub("^\\S+ \\S+ \\S+ ", "", x),
      function(v) as.numeric(strsplit(v, " ", fixed = TRUE)[[1]])), keys)
  }
  a <- parse(lines$image); b <- parse(lines$session)
  keys <- union(names(a), names(b))
  vapply(keys, function(k) {
    if (is.null(a[[k]]) || is.null(b[[k]]) || length(a[[k]]) != length(b[[k]])) return(Inf)
    max(abs(a[[k]] - b[[k]])) / max(1, max(abs(a[[k]])))
  }, numeric(1))
}

test_that("the package image computes what the engine compiled in the session computes", {
  skip_without_julia()
  # Needs its own engine build, a workload-free one in a temporary depot (8 s
  # on dev1), and then compiles every captured model's value, gradient and
  # Hessian in that session (about 8.5 minutes). Nothing cheaper sees the
  # image's code at all: every other test of a fit runs whatever the image
  # holds, and a workload-free engine has to be a second build.
  skip_unless_slow("building a workload-free engine to compare the image with")
  enabled <- tryCatch(isTRUE(ctsem:::.ctJuliaEval(
    "ContinuousTimeSEM._PRECOMPILE_WORKLOAD_ENABLED")), error = function(e) NA)
  skip_if(!isTRUE(enabled), "the engine image was built without the workload")
  lines <- .image_vs_session()
  expect_identical(tail(lines$image, 1), "done", label = paste(c("the image run finished",
    tail(lines$image, 5)), collapse = "\n"))
  expect_identical(tail(lines$session, 1), "done", label = "the workload-free run finished")
  expect_identical(lines$image[1], "workload true")
  expect_identical(lines$session[1], "workload false")
  diffs <- .image_vs_session_diff(lines)
  expect_gt(length(diffs), 20)
  # Rounding differences are allowed for (the image's code is built for the
  # CPU targets the system image lists, the session's for this CPU); anything
  # the source actually changes is orders of magnitude above this.
  expect_lt(max(diffs), 1e-8, label = paste("largest relative difference, at",
    names(diffs)[which.max(diffs)]))
})
