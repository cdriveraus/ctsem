# Ratchets, not thresholds.
#
# Two shapes of duplication in this package are countable, and the count is only
# useful as a direction of travel: it may go down, and it may not go up without
# someone deciding it should. Both numbers below were the state at the commit
# that added this file, and both are things a new function acquires by accident
# rather than by intent.
#
# If one of these fails after a change that *should* add a fork -- a genuinely
# new backend-specific report -- raise the number and say in the commit message
# which fork was added and why it could not go behind an existing seam. That
# sentence is the whole point of the test; the number is just what forces it to
# be written.
#
# Cheap and always runs: this reads source text, it does not fit anything.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

.dup_sources <- function() {
  rdir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"),
    mustWork = FALSE)
  if (!dir.exists(rdir)) return(NULL)
  files <- list.files(rdir, pattern = "[.][Rr]$", full.names = TRUE)
  stats::setNames(lapply(files, readLines, warn = FALSE), basename(files))
}

# Lines that are wholly a comment do not count: this is about code, and the
# comments in this package quote the very constructs being counted.
.dup_code_lines <- function(src) {
  lines <- unlist(src, use.names = FALSE)
  lines[!grepl("^\\s*#", lines)]
}

test_that("the backend is not asked by class literal in more places than it is", {
  src <- .dup_sources()
  skip_if(is.null(src), "package source not available")
  code <- .dup_code_lines(src)

  # `.ctFitIsJulia(fit)` (R/ctBackendKalman.R) is the named predicate for this
  # question. `inherits(fit, 'ctJuliaFit')` is the same question spelled out,
  # and it is spelled out far more often -- in both quote styles, which is how
  # a grep for one of them under-reports it.
  literal <- length(grep("inherits\\([^)]*ctJuliaFit", code))

  # The ratchet. Bring this down opportunistically when touching a file for
  # another reason; a sweep is not worth the churn, and this is what makes the
  # opportunistic version converge.
  expect_lte(literal, 56L)
})

test_that("the number of two-implementation forks does not grow", {
  rdir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"),
    mustWork = FALSE)
  skip_if(!dir.exists(rdir), "package source not available")

  # ctsem's fork idiom:
  #
  #     f <- function(fit, ...) {
  #       if (<predicate>) return(g(fit, ...))
  #       ...the other implementation, inline...
  #     }
  #
  # By construction `f`'s own body and `g` answer the same question, so each of
  # these is a place the package keeps two implementations of one thing. Some
  # are right -- the two backends genuinely differ -- and the number is here to
  # make adding one a decision rather than a habit.
  #
  # Walked over the parse tree, not grepped. The first version of this test used
  # a regex and counted 151 where the real number is 23, because `return(NULL)`,
  # `return(list(...))` and `return(invisible(x))` all look the same to a regex.
  # The callee has to be a function this package defines, which is the condition
  # that separates "delegates to the other implementation" from "returns early".
  # Same logic as dev/duplication/siblings.R, which is the authority if the two
  # ever disagree.
  files <- list.files(rdir, pattern = "[.][Rr]$", full.names = TRUE)
  defined <- character(0)
  asts <- list()
  for (f in files) {
    ex <- try(parse(f, keep.source = FALSE), silent = TRUE)
    if (inherits(ex, "try-error")) next
    asts[[f]] <- ex
    for (e in as.list(ex)) {
      if (is.call(e) && length(e) >= 3 &&
          as.character(e[[1]])[1] %in% c("<-", "=") &&
          is.call(e[[3]]) && as.character(e[[3]][[1]])[1] == "function") {
        defined <- c(defined, paste(deparse(e[[2]]), collapse = ""))
      }
    }
  }

  nforks <- 0L
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (identical(as.character(e[[1]])[1], "if") && length(e) >= 3) {
      arm <- e[[3]]
      if (is.call(arm) && identical(as.character(arm[[1]])[1], "{") &&
          length(arm) == 2L) arm <- arm[[2]]
      if (is.call(arm) && identical(as.character(arm[[1]])[1], "return") &&
          length(arm) == 2L) {
        inner <- arm[[2]]
        if (is.call(inner) && is.name(inner[[1]]) &&
            as.character(inner[[1]]) %in% defined) nforks <<- nforks + 1L
      }
    }
    for (k in seq_along(e)[-1]) if (is.call(e[[k]])) walk(e[[k]])
    invisible(NULL)
  }
  for (ex in asts) for (e in as.list(ex)) walk(e)

  expect_gt(length(defined), 500L)   # the walk found the package, not nothing
  expect_lte(nforks, 23L)
})

test_that("the duplication detectors are present and runnable", {
  # A detector nobody can find is a detector nobody runs. This fails if the
  # scripts are moved or renamed without the README and this test following.
  dev <- normalizePath(file.path(testthat::test_path(), "..", "..", "dev",
    "duplication"), mustWork = FALSE)
  skip_if(!dir.exists(dev), "dev/ not shipped in this build")
  for (f in c("siblings.R", "clones.R", "clones-jl.R", "README.md")) {
    expect_true(file.exists(file.path(dev, f)), info = f)
  }
  # And that they parse, since an R script that does not is found at the moment
  # someone reaches for it, which is the moment they least want to debug it.
  for (f in c("siblings.R", "clones.R", "clones-jl.R")) {
    expect_silent(invisible(parse(file.path(dev, f))))
  }
})
