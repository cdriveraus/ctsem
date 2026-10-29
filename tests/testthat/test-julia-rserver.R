# Sessions start from a precompiled copy of JuliaConnectoR's server
# (R/ctJuliaServer.R, inst/julia/rserver/): the installed JuliaConnectoR's own
# server files as a package, built once, with JuliaConnectoR's start as the
# fallback wherever anything does not go to plan.

test_that("the server's edits apply to the installed JuliaConnectoR, once each", {
  skip_if_not_installed("JuliaConnectoR")
  internals <- ctsem:::.ctJuliaServerInternals()
  # Internals that have moved are the fallback's case, said here so a
  # JuliaConnectoR update that moves them is seen rather than only slower.
  expect_false(is.null(internals))
  dir <- tempfile("rserver-layout")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  expect_true(ctsem:::.ctJuliaServerLayout(dir, internals))
  src <- file.path(dir, "CtsemRServer", "src", "connector")
  for (patch in ctsem:::.ctJuliaServerPatches) {
    text <- paste(readLines(file.path(src, patch[[1L]]), warn = FALSE), collapse = "\n")
    expect_true(grepl(patch[[3L]], text, fixed = TRUE), label = patch[[1L]])
    expect_false(grepl(patch[[2L]], text, fixed = TRUE), label = patch[[1L]])
  }
  expect_true(file.exists(file.path(dir, "serve.jl")))
  expect_true(file.exists(file.path(dir, "CtsemRServer", "Project.toml")))
})

test_that("without JuliaConnectoR's internals a session starts the old way", {
  cache <- ctsem:::.ct_julia_cache
  saved <- cache$rserver_off
  on.exit(assign("rserver_off", saved, envir = cache), add = TRUE)
  assign("rserver_off", NULL, envir = cache)
  local_mocked_bindings(.ctJuliaServerInternals = function() NULL)
  expect_null(ctsem:::.ctJuliaServerPrepare())
  # And it is not asked again this R session.
  expect_true(isTRUE(cache$rserver_off))
})

test_that("a depot JuliaConnectoR has not connected from builds nothing and records nothing", {
  skip_on_cran()
  skip_without_julia()
  internals <- ctsem:::.ctJuliaServerInternals()
  skip_if(is.null(internals))
  bin <- ctsem:::.ctJuliaBin()
  # This session's depots: JuliaConnectoR has connected from them.
  expect_true(ctsem:::.ctJuliaServerBuildable(internals, bin))
  # A first install's: no registry, no Tables. The offline build could only
  # fail there, and a recorded failure held the next week's sessions to the
  # old start, though the first one installs Tables.
  empty <- tempfile("depot")
  cache <- tempfile("cache")
  dir.create(empty)
  on.exit(unlink(c(empty, cache), recursive = TRUE), add = TRUE)
  withr::local_envvar(JULIA_DEPOT_PATH = empty, R_USER_CACHE_DIR = cache)
  expect_false(ctsem:::.ctJuliaServerBuildable(internals, bin))
  expect_null(ctsem:::.ctJuliaServerReady(internals, bin))
  expect_length(list.files(cache, pattern = "^rserver-", recursive = TRUE,
    include.dirs = TRUE), 0L)
})

test_that("a fresh session runs the precompiled server, with the real Tables", {
  skip_on_cran()
  skip_without_julia()
  # In a fresh R process: this one's Julia was started by the helper's own
  # JuliaConnectoR check, which ctsem adopts as it finds it.
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  from_source <- file.exists(file.path(root, "DESCRIPTION")) &&
    dir.exists(file.path(root, "R"))
  loader <- if (from_source) {
    sprintf("suppressMessages(pkgload::load_all(%s, compile = FALSE, quiet = TRUE))",
      deparse(root))
  } else "suppressMessages(library(ctsem))"
  out <- tempfile("rserver", fileext = ".rds")
  script <- tempfile("rserver", fileext = ".R")
  writeLines(c(loader,
    "suppressMessages(ctJuliaSetup(threads = 1L))",
    "res <- list(",
    "  server = as.character(ctsem:::.ctJuliaEval('string(parentmodule(RConnector))')),",
    "  dummy = isTRUE(ctsem:::.ctJuliaEval('isdefined(Tables, :JuliaConnectoR_DummyTables)')),",
    "  sum = as.numeric(ctsem:::.ctJuliaEval('sum([1.5, 2.5])')),",
    "  pid = ctsem:::.ct_julia_cache$pid)",
    sprintf("saveRDS(res, %s)", deparse(out))), script)
  log <- withr::with_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
      stdout = TRUE, stderr = TRUE)))
  expect_true(file.exists(out), info = paste(log, collapse = "\n"))
  skip_if(!file.exists(out))
  res <- readRDS(out)
  expect_identical(res$server, "CtsemRServer")
  expect_false(res$dummy)
  expect_equal(res$sum, 4)
  # Recorded as ctsem's own session, not as one someone else started.
  expect_false(is.null(res$pid))
})
