# Fast local test loop.
#
#   source("dev/test-local.R")   # once per session
#   tt("julia-backend")          # files matching a pattern
#   tt_changed()                 # files for the R code you have edited
#   tt_fast()                    # the quick tier
#
# Why a session rather than a script per file. The julia engine precompiles on
# first use, about two minutes, and that cost is per R SESSION, not per file.
# A loop of `Rscript -e 'test_file(...)'` calls pays it every time and spends
# most of its wall clock there. Keeping one session and calling tt() repeatedly
# pays it once. This is the single largest speedup available locally and it
# needs no changes to any test.
#
# NOT_CRAN is set here because every julia test is skip_on_cran: without it a
# run reports zero assertions and looks identical to a clean pass.

Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("CTSEM_JULIA_AGREE"))) Sys.setenv(CTSEM_JULIA_AGREE = "yes")

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is needed for the local loop")
}

# compile = FALSE is deliberate: a debug rebuild of the stan exports cannot
# succeed in this environment and it DELETES the objects on the way out. Use
# R CMD INSTALL when C++ genuinely changed.
message("loading ctsem (compile = FALSE) ...")
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
message("ready. tt('pattern'), tt_changed(), tt_fast()")

# Both separators. testthat collects anything matching `^test`, and one file is
# spelled `test_behavGenNLcor.R`, so a glob of `test-*.R` alone left one of the
# most expensive files in the suite out of tt(), tt_changed(), tt_fast() and
# tt_time() -- silently, since a file that is never selected looks exactly like
# a file that passed.
.tt_files <- function() sort(basename(c(Sys.glob("tests/testthat/test-*.R"),
  Sys.glob("tests/testthat/test_*.R"))))

.tt_run <- function(files) {
  if (!length(files)) { message("no matching test files"); return(invisible(NULL)) }
  res <- data.frame(file = files, secs = NA_real_, pass = NA_integer_,
    fail = NA_integer_, skip = NA_integer_, stringsAsFactors = FALSE)
  for (i in seq_along(files)) {
    t0 <- Sys.time()
    r <- try(testthat::test_file(file.path("tests/testthat", files[i]),
      reporter = "silent"), silent = TRUE)
    res$secs[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(r, "try-error")) {
      cat(sprintf("%-44s %6.1fs  ABORTED\n", files[i], res$secs[i]))
      next
    }
    d <- as.data.frame(r)
    res$pass[i] <- sum(d$passed); res$fail[i] <- sum(d$failed)
    res$skip[i] <- sum(d$skipped)
    # A file reporting zero assertions is a failure, not a pass: it means every
    # test skipped, which is what an unset NOT_CRAN or a missing julia looks like.
    flag <- if (sum(d$failed) > 0) "FAIL" else if (sum(d$passed) == 0) "ZERO ASSERTIONS" else "ok"
    cat(sprintf("%-44s %6.1fs  pass=%-5d fail=%-3d skip=%-3d %s\n",
      files[i], res$secs[i], sum(d$passed), sum(d$failed), sum(d$skipped), flag))
    bad <- d[d$failed > 0 | d$error > 0, "test"]
    if (length(bad)) cat("    failing:", paste(bad, collapse = " | "), "\n")
    flush.console()
  }
  cat(sprintf("-- %d file(s), %.0fs, %d pass, %d fail\n", nrow(res),
    sum(res$secs, na.rm = TRUE), sum(res$pass, na.rm = TRUE),
    sum(res$fail, na.rm = TRUE)))
  invisible(res)
}

#' Run test files whose name matches a pattern.
tt <- function(pattern = NULL) {
  f <- .tt_files()
  if (!is.null(pattern)) f <- grep(pattern, f, value = TRUE)
  .tt_run(f)
}

#' Run the test files most likely to cover the R files you have edited.
#'
#' Maps a changed `R/ctFoo.R` to any test file whose name shares a stem with it,
#' which is a heuristic and not a coverage guarantee: it will miss a test that
#' exercises your change under an unrelated name. It is meant for the inner loop,
#' not as the check before you commit.
tt_changed <- function(ref = "HEAD") {
  changed <- system(paste("git diff --name-only", shQuote(ref), "-- R/"),
    intern = TRUE)
  changed <- c(changed, system("git diff --name-only --cached -- R/", intern = TRUE))
  changed <- unique(basename(changed[nzchar(changed)]))
  if (!length(changed)) { message("no changed files under R/"); return(invisible(NULL)) }
  stems <- tolower(sub("\\.R$", "", changed))
  stems <- sub("^ct", "", stems)
  f <- .tt_files()
  keep <- vapply(f, function(x) {
    any(vapply(stems, function(s) nchar(s) > 3 &&
      grepl(s, tolower(x), fixed = TRUE), logical(1)))
  }, logical(1))
  message("changed: ", paste(changed, collapse = ", "))
  .tt_run(f[keep])
}

#' The quick tier: files that do not fit models, or fit only trivial ones.
#'
#' Populated from a measured timing pass rather than by guess; see
#' dev/test-timings.csv. Regenerate that with tt_time() when fixtures change,
#' since a tier built on stale timings quietly stops being quick.
tt_fast <- function(max_secs = 20) {
  p <- "dev/test-timings.csv"
  if (!file.exists(p)) {
    stop("no dev/test-timings.csv -- run tt_time() once to measure, ",
      "then commit it")
  }
  tm <- utils::read.csv(p, stringsAsFactors = FALSE)
  f <- tm$file[!is.na(tm$secs) & tm$secs <= max_secs]
  message(length(f), " file(s) under ", max_secs, "s")
  .tt_run(intersect(f, .tt_files()))
}

#' Measure every file and write dev/test-timings.csv.
#'
#' Do this on a quiet machine. Timings taken while several jobs are running are
#' not merely noisy, they are wrong by enough to put files in the wrong tier.
tt_time <- function() {
  res <- .tt_run(.tt_files())
  utils::write.csv(res[, c("file", "secs", "pass", "fail", "skip")],
    "dev/test-timings.csv", row.names = FALSE)
  message("wrote dev/test-timings.csv")
  invisible(res)
}
