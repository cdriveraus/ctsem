# Escape during a julia call, and the session start (R/ctJuliaBridge.R).
#
# A real Escape cannot be sent from inside a test, so the state it leaves is
# built directly: the request is written and the reply walked away from, which
# is exactly what the interrupt handler does when Escape lands while R waits.
# What is asserted is everything after that -- the reply collected rather than
# read by the wrong call, the engine stopping at a checkpoint, the session and
# its process kept -- which is where the old behaviour went wrong.

send_and_walk_away <- function(expr) {
  ns <- asNamespace("ctsem")
  wire <- ns$.ctJuliaWire()
  ns$.ctJuliaEnsureStarted(wire)
  con <- wire$pkgLocal$con
  writeBin(wire$CALL, con)
  wire$writeString("RConnector.mainevalcmd")
  wire$writeList(list(expr))
  ns$.ctJuliaAbandon("RConnector.mainevalcmd")
}

process_exists <- function(pid) {
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system2("tasklist",
      c("/FI", shQuote(paste("PID eq", pid)), "/NH"), stdout = TRUE))
    any(grepl(paste0("\\b", pid, "\\b"), out))
  } else {
    isTRUE(tools::pskill(pid, 0L))
  }
}

test_that("a reply left owed by Escape is collected before the next call", {
  skip_without_julia()
  ctJuliaSetup()
  cache <- ctsem:::.ct_julia_cache
  pid <- cache$pid
  expect_true(is.numeric(pid) && !is.na(pid))
  send_and_walk_away("sleep(2); :abandoned")
  expect_false(is.null(cache$inflight))
  expect_true(file.exists(cache$interrupt_file))
  # Read by the next call and thrown away; before this existed, the next call
  # read it as its own answer.
  expect_message(x <- ctsem:::.ctJuliaEval("1 + 1"), "Waiting for Julia")
  expect_equal(x, 2)
  expect_null(cache$inflight)
  expect_false(file.exists(cache$interrupt_file))
  # Nothing was restarted, so nothing compiled was lost.
  expect_identical(cache$pid, pid)
})

test_that("a restore deferred behind an owed reply runs before the next call", {
  skip_without_julia()
  ctJuliaSetup()
  on.exit(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!", 0L), add = TRUE)
  send_and_walk_away("sleep(1); 0")
  # As `.ctBackendRestoreMaxChunks()` does from `on.exit` while Escape unwinds.
  expect_null(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!", 3L,
    .defer = TRUE))
  expect_length(ctsem:::.ct_julia_cache$deferred, 1L)
  now <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_max_chunks()"))
  expect_equal(now$max_chunks, 3)
  expect_null(ctsem:::.ct_julia_cache$deferred)
})

test_that("the engine stops at its next checkpoint once asked to", {
  skip_without_julia()
  ctJuliaSetup()
  # A minute of iterations through the same hook the optimiser records its
  # trace with. Without the checkpoint the next call would wait all of it.
  send_and_walk_away(paste0("let t = ContinuousTimeSEM.CTSEMTrace(:a); ",
    "for i in 1:6000; ContinuousTimeSEM._record!(t, i, 1.0); sleep(0.01); end; ",
    ":finished end"))
  elapsed <- system.time(x <- suppressMessages(ctsem:::.ctJuliaEval("2 + 2")))[["elapsed"]]
  expect_equal(x, 4)
  expect_lt(elapsed, 30)
})

test_that("an argument that is itself a Julia call is sent before the call it is for", {
  skip_without_julia()
  ctJuliaSetup()
  # Left as a promise, the inner call used to run while the outer request was
  # half written, and both ends then waited on each other indefinitely.
  expect_equal(ctsem:::.ctJuliaCall("sum", ctsem:::.ctJuliaPut(c(1, 2, 3))), 6)
  # And through the engine module, as every fit calls it.
  module <- ctsem:::.ctJuliaModule()
  expect_equal(module$scalar_square(ctsem:::.ctJuliaEval("3.0")), 9)
})

test_that("an interrupt raised in Julia arrives as an interrupt, not an error", {
  skip_without_julia()
  ctJuliaSetup()
  got <- tryCatch(ctsem:::.ctJuliaEval("throw(InterruptException())"),
    interrupt = function(e) "interrupt", error = function(e) conditionMessage(e))
  expect_identical(got, "interrupt")
  # The reply was read whole, so the session carries on.
  expect_equal(ctsem:::.ctJuliaEval("5"), 5)
})

test_that("stopping Julia outright ends its process and the next call starts anew", {
  skip_without_julia()
  ctJuliaSetup()
  cache <- ctsem:::.ct_julia_cache
  old <- cache$pid
  ctsem:::.ctJuliaKill()
  Sys.sleep(1)
  expect_false(process_exists(old))
  msgs <- testthat::capture_messages(x <- ctsem:::.ctJuliaEval("3 + 3"))
  expect_equal(x, 6)
  expect_true(any(grepl("^Starting Julia \\d+\\.\\d+\\.\\d+ with \\d+ threads? \\.\\.\\.",
    msgs)), info = paste(msgs, collapse = " | "))
  expect_false(identical(cache$pid, old))
  ctJuliaSetup()
})

test_that("the engine tells a running R process from one that has exited", {
  skip_without_julia()
  ctJuliaSetup()
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  gone <- as.integer(system2(rscript, c("-e", shQuote("cat(Sys.getpid())")), stdout = TRUE))
  alive <- function(pid) ctsem:::.ctJuliaEval(sprintf(paste0(
    "let p = ContinuousTimeSEM._CTSEM_PARENT_PID, old = p[]; p[] = %d; ",
    "r = ContinuousTimeSEM._ctsem_parent_alive(); p[] = old; r end"), pid))
  expect_true(alive(Sys.getpid()))
  expect_false(alive(gone))
})

test_that("a newer Julia is offered only within the running or tested series", {
  newer <- function(version, managed = FALSE, latest = function(series) NA_character_) {
    got <- ctsem:::.ctJuliaNewerRelease(version, bin = "x", managed = managed,
      pinned = "1.12.7", latest = latest)
    if (is.null(got)) NA_character_ else paste0(got$version, got$how)
  }
  lookup <- function(table) function(series) {
    if (series %in% names(table)) table[[series]] else NA_character_
  }
  # A patch of the running series, from the pin or from juliaup's list.
  expect_identical(newer("1.12.5"), "1.12.7")
  expect_identical(newer("1.12.5", latest = lookup(c(`1.12` = "1.12.8"))), "1.12.8")
  # An older series is pointed at the tested one, and says so.
  expect_identical(newer("1.11.9"), "1.12.7 (ctsem is tested on 1.12)")
  # A newer series than the tested one hears only about its own patches.
  expect_identical(newer("1.13.1", latest = lookup(c(`1.13` = "1.13.1"))), NA_character_)
  expect_identical(newer("1.13.1", latest = lookup(c(`1.13` = "1.13.2"))), "1.13.2")
  # Nothing newer than what is running.
  expect_identical(newer("1.12.7"), NA_character_)
  # A Julia ctsem installed is offered the pin, with the call that installs it.
  expect_identical(newer("1.12.5", managed = TRUE,
    latest = lookup(c(`1.12` = "1.12.9"))),
    "1.12.7; ctJuliaInstall(force = TRUE) installs it")
  expect_identical(newer(NA_character_), NA_character_)
})

test_that("the thread count a session will start with is read as Julia reads it", {
  threads <- function(value) withr::with_envvar(c(JULIA_NUM_THREADS = value),
    ctsem:::.ctJuliaStartThreads())
  expect_identical(threads(NA), "1")
  expect_identical(threads("4"), "4")
  expect_identical(threads("4,1"), "4")
  expect_identical(threads("auto"), "auto")
  expect_identical(threads("x"), NA_character_)
})
