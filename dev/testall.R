# Run the whole test suite against the INSTALLED ctsem, in parallel processes.
#
#   source("dev/testall.R")   # from the package root
#   testall()                 # 4 workers, examples first
#   testall(cores = 1)        # one process, verbose
#   testall(examples = FALSE)
#
# This is the complement to dev/test-local.R, not a duplicate of it, and which
# one to reach for follows from the difference:
#
#   dev/test-local.R  one R session, devtools::load_all(compile = FALSE), a
#                     pattern or a changed-file heuristic. The inner loop --
#                     the julia engine precompiles once per session and tt()
#                     reuses it.
#   dev/testall.R     worker processes, library(ctsem), every file including
#                     the examples. The before-release run -- it tests what a
#                     user would install, and it is the only one of the two
#                     that runs the examples.
#
# Because the workers call library(ctsem), install the tree under test first:
# a green run here against a stale install proves nothing about the working
# tree. See dev/RELEASE-CHECKLIST.md.
#
# It lived in R/ until 2026-09, where it shipped to CRAN and was byte-compiled
# at install for no reason: unexported, and its only reference anywhere was
# its own name inside its own error message.

testall <- function(cores = 4, folder = '/tests/testthat', examples = TRUE) {
  if (!requireNamespace('testthat', quietly = TRUE)) {
    stop("Package 'testthat' is required to run testall().", call. = FALSE)
  }
  .testall_setup <- function(testfolder) {
    Sys.setenv(NOT_CRAN = 'true')
    suppressPackageStartupMessages(library(ctsem))
    supportfiles <- list.files(testfolder, pattern = '^(helper|setup).*\\.[rR]$',
      full.names = TRUE)
    for (supportfile in supportfiles) sys.source(supportfile, envir = globalenv())
    pdf(NULL)
    invisible(TRUE)
  }
  testfolder <- normalizePath(paste0('.', folder))
  Sys.setenv(NOT_CRAN = 'true')
  tests <- dir(paste0('.', folder))
  tests <- tests[grepl('^test', tests)]
  # `x[-integer(0)]` is empty, not everything, so an unmatched grep here
  # silently reduced the whole list to nothing. That is what happened when
  # test-runExamples.R moved out of tests/testthat: this helper ran zero
  # tests in every mode and said so only by finishing instantly.
  runex <- grep('runExamples', tests)
  if (length(runex)) {
    tests <- if (examples) c(tests[runex], tests[-runex]) else tests[-runex]
  }
  a <- Sys.time()

  if (cores > 1) {
    # rscript_libs: without it a worker searches its own default .libPaths(),
    # not the caller's, so library(ctsem) below can silently load a different
    # install than the one running this code -- exactly wrong for a test
    # runner, whose whole point is to test the development tree under test
    # rather than whatever ctsem happens to be installed globally. Same
    # failure mode makeClusterID() (stanoptimis.R) was fixed for.
    cl <- parallelly::makeClusterPSOCK(cores, rscript_libs = .libPaths())
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    out <- parallel::parLapplyLB(cl, paste0(getwd(), folder, '/', tests), function(x, testfolder) {
      Sys.setenv(NOT_CRAN = 'true')
      suppressPackageStartupMessages(library(ctsem))
      supportfiles <- list.files(testfolder, pattern = '^(helper|setup).*\\.[rR]$',
        full.names = TRUE)
      for (supportfile in supportfiles) sys.source(supportfile, envir = globalenv())
      pdf(NULL)
      on.exit(dev.off(), add = TRUE)
      out <- testthat::test_file(x, reporter = "minimal")
      return(out)
    }, testfolder = testfolder)
  }
  if (cores == 1) {
    .testall_setup(testfolder)
    out <- lapply(paste0(getwd(), folder, '/', tests), function(x) {
      cat(x)
      out <- testthat::test_file(x, reporter = "minimal")
      print(out)
      return(out)
    })
  }
  out2 <- do.call(what = rbind, lapply(out, utils::getS3method('as.data.frame', 'testthat_results')))
  if (dev.cur() > 1) dev.off()
  print(out2[, colnames(out2) != 'result'])
  print(Sys.time() - a)
  if (cores > 1) parallel::stopCluster(cl)
  return(invisible(out2))
}
