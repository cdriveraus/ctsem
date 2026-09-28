# ctJuliaProcesses() and ctJuliaKill(): telling ctsem's Julia processes from
# any other, and stopping the ones nothing should be running. Every test points
# the record at a temporary directory, so none of it touches the record of the
# machine's real sessions.

gone_pid <- function() {
  rscript <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  as.integer(system2(rscript, c("-e", shQuote("cat(Sys.getpid())")), stdout = TRUE))
}

test_that("records of processes that have gone, or are not Julia, are dropped", {
  withr::local_options(ctsem.julia.registry = withr::local_tempdir())
  ctsem:::.ctJuliaRegister(gone_pid(), 2L)
  # This R process is alive, and is not Julia: an id that has been handed on.
  ctsem:::.ctJuliaRegister(Sys.getpid(), 2L)
  expect_length(list.files(getOption("ctsem.julia.registry")), 2L)
  expect_equal(nrow(ctsem:::.ctJuliaRegistry()), 0L)
  expect_length(list.files(getOption("ctsem.julia.registry")), 0L)
})

test_that("without ps, R sessions still count as running", {
  # The fallback once listed Julia processes only, so every other session's R
  # looked ended and its engine an orphan -- which ctJuliaKill() stops by
  # default.
  info <- ctsem:::.ctProcInfo(Sys.getpid(), use_ps = FALSE)
  expect_true(info$alive)
  expect_match(info$name, "^R", ignore.case = TRUE)
  expect_false(ctsem:::.ctProcInfo(gone_pid(), use_ps = FALSE)$alive)
})

test_that("this session's engine is listed with the threads it started with", {
  skip_without_julia()
  skip_if_not_installed("ps")
  withr::local_options(ctsem.julia.registry = withr::local_tempdir())
  ctsem:::.ctJuliaClearSession()
  suppressMessages(ctJuliaSetup())
  cache <- ctsem:::.ct_julia_cache
  procs <- ctJuliaProcesses(interval = 0)
  mine <- procs[procs$role == "this session", , drop = FALSE]
  expect_true(cache$pid %in% mine$pid)
  expect_equal(mine$julia_threads[mine$pid == cache$pid],
    as.integer(ctsem:::.ctJuliaEval("Threads.nthreads()")))
})

test_that("an engine whose R session has ended is reported and stopped by default", {
  skip_without_julia()
  skip_if_not_installed("ps")
  registry <- withr::local_tempdir()
  withr::local_options(ctsem.julia.registry = registry)
  ctsem:::.ctJuliaClearSession()
  suppressMessages(ctJuliaSetup())
  cache <- ctsem:::.ct_julia_cache
  engine <- cache$pid
  # Its record now says it serves an R process that has exited.
  file <- list.files(registry, full.names = TRUE)
  expect_length(file, 1L)
  record <- read.dcf(file)
  record[, "r_pid"] <- gone_pid()
  write.dcf(record, file)
  # And this session forgets it, as the ended one would have: its own engine is
  # never an orphan, whatever a record says.
  cache$pid <- NULL
  procs <- ctJuliaProcesses(interval = 0)
  expect_identical(procs$role[procs$pid == engine], "orphaned")
  cache$orphans_checked <- NULL
  expect_message(ctsem:::.ctJuliaWarnOrphans(), "still running for R sessions that have ended")
  stopped <- suppressMessages(ctJuliaKill())
  expect_true(engine %in% stopped$pid[stopped$stopped])
  Sys.sleep(1)
  expect_false(isTRUE(ctsem:::.ctProcInfo(engine)$alive))
  # It was in fact this session's engine, so the session notices and starts
  # again: the first call says the process went, the next starts a new one.
  try(ctsem:::.ctJuliaEval("0"), silent = TRUE)
  expect_equal(ctsem:::.ctJuliaEval("1 + 1"), 2)
})

test_that("stopping this session's engine forgets it, and Escape says when work goes on", {
  skip_without_julia()
  withr::local_options(ctsem.julia.registry = withr::local_tempdir())
  suppressMessages(ctJuliaSetup())
  cache <- ctsem:::.ct_julia_cache
  wire <- ctsem:::.ctJuliaWire()
  con <- wire$pkgLocal$con
  writeBin(wire$CALL, con)
  wire$writeString("RConnector.mainevalcmd")
  wire$writeList(list("sleep(5); 1"))
  expect_message(ctsem:::.ctJuliaAbandon("RConnector.mainevalcmd"), "in the background")
  suppressMessages(ctJuliaKill("session"))
  expect_null(cache$pid)
  expect_null(cache$inflight)
  expect_length(list.files(getOption("ctsem.julia.registry")), 0L)
  expect_equal(suppressMessages(ctsem:::.ctJuliaEval("2 + 2")), 4)
})
