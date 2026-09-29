# Fast local test loop.
#
#   source("dev/test-local.R")   # once per session
#   tt("julia-backend")          # files matching a pattern
#   tt_changed()                 # files for the R code you have edited
#   tt_fast()                    # the quick tier
#   tt_area("laplace")           # the core tests of the area a change touches
#
# Which to run. The whole suite is about four hours, so it is for a release or
# a change that cuts across many areas, not for every change: run tt_fast()
# and the core set of each area the change touches (tt_area(), listed in
# dev/test-areas.csv), and add the area's extended set when the change is deep.
#
# Why a session rather than a script per file. Each R session loads the engine
# and compiles every model shape its fits use that the package image does not
# hold -- seconds for the structures the image captures, a minute or more for
# the rest. A loop of `Rscript -e 'test_file(...)'` calls pays that for every
# file; one session calling tt() repeatedly pays it once per shape.
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
message("ready. tt('pattern'), tt_changed(), tt_fast(), tt_area('area')")

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
    # A test that errors is a failure too. It was counted nowhere, so a file
    # whose one test died before its first expectation printed as ZERO
    # ASSERTIONS, the flag for a file that skipped, rather than as a failure --
    # test_behavGenNLcor.R, whose stan compile fails on Windows, read that way.
    errored <- sum(d$error)
    res$pass[i] <- sum(d$passed); res$fail[i] <- sum(d$failed) + errored
    res$skip[i] <- sum(d$skipped)
    # A file reporting zero assertions is a failure, not a pass: it means every
    # test skipped, which is what an unset NOT_CRAN or a missing julia looks
    # like. The reason is printed with it, so a release-tier file skipping for
    # the reason it states can be told from one that skipped for no good one.
    flag <- if (res$fail[i] > 0) "FAIL" else
      if (sum(d$passed) == 0) "ZERO ASSERTIONS" else "ok"
    cat(sprintf("%-44s %6.1fs  pass=%-5d fail=%-3d skip=%-3d %s\n",
      files[i], res$secs[i], sum(d$passed), res$fail[i], sum(d$skipped), flag))
    bad <- d[d$failed > 0 | d$error > 0, "test"]
    if (length(bad)) cat("    failing:", paste(bad, collapse = " | "), "\n")
    if (identical(flag, "ZERO ASSERTIONS")) {
      why <- unlist(lapply(r, function(t) vapply(t$results, function(e)
        if (inherits(e, "expectation_skip")) conditionMessage(e) else NA_character_,
        character(1))))
      why <- unique(why[!is.na(why)])
      if (length(why)) cat("    skipped:", paste(why, collapse = " | "), "\n")
    }
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
  have <- .tt_files()
  # A file with no row is never selected, which looks exactly like a file that
  # is slow; say so, or a stale csv shrinks the tier without anyone noticing.
  untimed <- setdiff(have, tm$file)
  gone <- setdiff(tm$file, have)
  if (length(untimed)) message(length(untimed), " file(s) have no timing and ",
    "are left out -- rerun tt_time():\n  ", paste(untimed, collapse = "\n  "))
  if (length(gone)) message(length(gone), " timed file(s) no longer exist: ",
    paste(gone, collapse = ", "))
  f <- intersect(tm$file[!is.na(tm$secs) & tm$secs <= max_secs], have)
  # A release-tier file measures as quick because it skipped; it is not quick.
  f <- setdiff(f, .tt_release_only())
  message(length(f), " file(s) under ", max_secs, "s")
  .tt_run(f)
}

# Files whose every entry in dev/test-areas.csv is the release tier.
.tt_release_only <- function(map = NULL) {
  if (is.null(map)) {
    if (!file.exists("dev/test-areas.csv")) return(character())
    map <- utils::read.csv("dev/test-areas.csv", stringsAsFactors = FALSE)
  }
  tiers <- tapply(map$tier, map$file, function(t) all(t == "release"))
  names(tiers)[tiers]
}

# The test blocks gated by skip_unless_slow() in each of `files`: release-tier
# checks living inside a core or extended file.
.tt_slow_blocks <- function(files) {
  n <- vapply(files, function(f) {
    lines <- readLines(file.path("tests/testthat", f), warn = FALSE)
    sum(grepl("^\\s*skip_unless_slow\\(", lines))
  }, integer(1))
  n[n > 0L]
}

#' The core tests for the areas a change touches.
#'
#' `dev/test-areas.csv` gives every test file one or more areas, and within
#' each a tier: `core`, the files a change in that area runs; `extended`, the
#' rest of what bears on it, for a change that goes deep; and `release`, checks
#' that cannot be made cheap -- a recovery that needs hundreds of subjects to
#' mean anything, a stan program compiled at test time -- run before a release
#' and weekly in CI. A release check skips unless CTSEM_SLOW_TESTS=true
#' (`skip_unless_slow()`, helper-julia.R), and the same gate marks single
#' release blocks inside core and extended files, so `release = TRUE` sets it
#' for the run and adds the release files to the extended ones. Without it,
#' `extended = TRUE` says what it left out. Every file is in some area, so all
#' areas at `release = TRUE` is the whole suite. Called with no area, lists the
#' areas and what each costs, from dev/test-timings.csv. The quick tier runs
#' first unless `fast = FALSE`: it is about a minute and covers the
#' specification code every area builds on.
tt_area <- function(area = NULL, extended = FALSE, fast = TRUE, release = FALSE) {
  map <- utils::read.csv("dev/test-areas.csv", stringsAsFactors = FALSE)
  # A file in no area is never selected, which looks exactly like a file that
  # passed; say so, as tt_fast() does for an untimed one.
  unmapped <- setdiff(.tt_files(), map$file)
  if (length(unmapped)) message(length(unmapped), " test file(s) are in no area ",
    "and never run here -- add them to dev/test-areas.csv:\n  ",
    paste(unmapped, collapse = "\n  "))
  tm <- if (file.exists("dev/test-timings.csv"))
    utils::read.csv("dev/test-timings.csv", stringsAsFactors = FALSE) else NULL
  minutes <- function(files) if (is.null(tm)) NA_real_ else
    round(sum(tm$secs[match(files, tm$file)], na.rm = TRUE) / 60, 1)
  if (is.null(area)) {
    out <- do.call(rbind, lapply(unique(map$area), function(a) {
      m <- map[map$area == a, ]
      data.frame(area = a, core_files = sum(m$tier == "core"),
        core_min = minutes(m$file[m$tier == "core"]),
        extended_min = minutes(m$file[m$tier != "release"]),
        release_files = sum(m$tier == "release"))
    }))
    print(out, row.names = FALSE)
    message("release_files are timed as skipped; tt_area(area, release = TRUE) runs them")
    return(invisible(out))
  }
  unknown <- setdiff(area, map$area)
  if (length(unknown)) stop("no such area: ", paste(unknown, collapse = ", "),
    ". Areas: ", paste(unique(map$area), collapse = ", "))
  tiers <- c("core", if (extended || release) "extended", if (release) "release")
  m <- map[map$area %in% area & map$tier %in% tiers, ]
  files <- intersect(.tt_files(), unique(m$file))
  if (fast && !is.null(tm)) {
    quick <- intersect(tm$file[!is.na(tm$secs) & tm$secs <= 20], .tt_files())
    if (!release) quick <- setdiff(quick, .tt_release_only(map))
    files <- union(quick, files)
  }
  if (release) {
    previous <- Sys.getenv("CTSEM_SLOW_TESTS", unset = NA)
    Sys.setenv(CTSEM_SLOW_TESTS = "true")
    on.exit(if (is.na(previous)) Sys.unsetenv("CTSEM_SLOW_TESTS") else
      Sys.setenv(CTSEM_SLOW_TESTS = previous), add = TRUE)
  } else if (extended) {
    # What `release = TRUE` would add: silence here reads as "that was all".
    left <- setdiff(intersect(.tt_files(),
      map$file[map$area %in% area & map$tier == "release"]), files)
    blocks <- .tt_slow_blocks(files)
    if (length(left) || length(blocks)) message("left out, as the release tier ",
      "(tt_area(..., release = TRUE) runs it):",
      if (length(left)) paste0("\n  files: ", paste(left, collapse = ", ")),
      if (length(blocks)) paste0("\n  skip_unless_slow() blocks: ",
        paste0(names(blocks), " (", blocks, ")", collapse = ", ")))
  }
  message(length(files), " file(s), about ", minutes(files), " min by dev/test-timings.csv")
  .tt_run(files)
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
