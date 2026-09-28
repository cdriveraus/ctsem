#!/usr/bin/env Rscript

# Render every Quarto/R Markdown vignette in a temporary directory using a
# PSOCK cluster. Completed artifacts are then copied to dev/vignettes. In
# RStudio, open this file, change interactive_cores if needed, and click
# Source. From a terminal: Rscript dev/render-vignettes-parallel.R --cores 4
#
# Each vignette gets its OWN staged directory. Quarto creates a .quarto
# scratch directory beside the input it is rendering and deletes it on the
# way out, so two workers sharing one directory race: whoever finishes first
# removes it under the others, and they die with "The process cannot access
# the file because it is being used by another process (os error 32)" after
# having rendered the document perfectly well. One directory per vignette is
# the whole fix; do not consolidate them again.

# Used only when this file is sourced interactively (for example, in RStudio).
interactive_cores <- 10L
#cores=10
arguments <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("Usage: Rscript dev/render-vignettes-parallel.R [--cores N | --cores=N | N]\n")
}

if (any(arguments %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

core_argument <- if (length(arguments) == 0L) {
  NA_character_
} else if (length(arguments) == 1L) {
  sub("^--cores=", "", arguments[[1L]])
} else if (length(arguments) == 2L && identical(arguments[[1L]], "--cores")) {
  arguments[[2L]]
} else {
  NA_character_
}

if (length(arguments) > 2L ||
    (length(arguments) == 2L && !identical(arguments[[1L]], "--cores")) ||
    (length(arguments) == 1L && identical(arguments[[1L]], "--cores"))) {
  usage()
  stop("Specify at most one core count.", call. = FALSE)
}

if (!requireNamespace("quarto", quietly = TRUE)) {
  stop("The 'quarto' R package must be installed to render the vignettes.", call. = FALSE)
}
if (!requireNamespace("ctsem", quietly = TRUE)) {
  stop("Install ctsem before rendering its vignettes.", call. = FALSE)
}

find_script_path <- function() {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) == 1L) {
    return(normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE))
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    editor_path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(editor_path) && file.exists(editor_path)) {
      return(normalizePath(editor_path, mustWork = TRUE))
    }
  }

  candidate <- file.path(getwd(), "dev", "render-vignettes-parallel.R")
  if (file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }

  stop(
    "Could not determine this script's location. In RStudio, save and source it ",
    "from the ctsem project.",
    call. = FALSE
  )
}

script_path <- find_script_path()
package_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
vignettes_dir <- file.path(package_root, "vignettes")
output_dir <- file.path(package_root, "dev", "vignettes")

publish_files <- function(files, output_dir) {
  if (!length(files)) {
    return(invisible(NULL))
  }
  copied <- file.copy(files, output_dir, recursive = TRUE, overwrite = TRUE)
  if (!all(copied)) {
    warning("Could not copy to ", output_dir, ": ",
            paste(basename(files[!copied]), collapse = ", "), call. = FALSE)
  }
  invisible(files[copied])
}

# Renders one vignette in its own staged directory and returns where the
# artifacts landed. Everything it prints, including knitr's progress and any
# quarto error, goes to <name>.log in that directory.
render_one <- function(vignette) {
  work_dir <- dirname(vignette)
  name <- tools::file_path_sans_ext(basename(vignette))
  log_file <- file.path(work_dir, paste0(name, ".log"))
  log_connection <- file(log_file, open = "wt")
  sink(log_connection)
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_connection)
  }, add = TRUE)

  message("Rendering ", basename(vignette), " at ", format(Sys.time(), tz = "UTC"))
  result <- tryCatch(
    {
      quarto::quarto_render(
        input = vignette,
        output_format = "html",
        as_job = FALSE
      )
      list(vignette = vignette, success = TRUE, log_file = log_file, error = NA_character_)
    },
    error = function(error) {
      message(conditionMessage(error))
      list(vignette = vignette, success = FALSE, log_file = log_file,
           error = conditionMessage(error))
    }
  )
  message("Finished ", basename(vignette), " at ", format(Sys.time(), tz = "UTC"))

  # The document may be on disk even when quarto exited non-zero, so report
  # what exists rather than what the status code implies.
  result$html_file <- file.path(work_dir, paste0(name, ".html"))
  if (!file.exists(result$html_file)) {
    result$html_file <- NA_character_
  }
  result
}

# Wrapped in a function so that on.exit actually runs: at the top level of a
# sourced script it registers against nothing, which left temporary render
# directories and PSOCK clusters behind on every interactive run.
render_vignettes <- function(vignettes_dir, output_dir, requested_cores) {
  run_dir <- tempfile("ctsem-vignettes-", tmpdir = tempdir())
  on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)

  source_entries <- list.files(vignettes_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  names <- tools::file_path_sans_ext(basename(
    grep("\\.(qmd|rmd)$", source_entries, value = TRUE, ignore.case = TRUE)
  ))
  if (!length(names)) {
    stop("No .qmd or .Rmd vignettes were found in ", vignettes_dir, ".", call. = FALSE)
  }

  # One staged copy of the whole vignettes directory per vignette, so each
  # render owns its .quarto scratch directory and its intermediate files.
  vignettes <- vapply(names, function(name) {
    work_dir <- file.path(run_dir, name)
    dir.create(work_dir, recursive = TRUE)
    if (!all(file.copy(source_entries, work_dir, recursive = TRUE, copy.date = TRUE))) {
      stop("Could not copy the vignette sources to ", work_dir, ".", call. = FALSE)
    }
    matched <- list.files(work_dir, pattern = paste0("^", name, "\\.(qmd|rmd)$"),
                          full.names = TRUE, ignore.case = TRUE)
    normalizePath(matched[[1L]], winslash = "/", mustWork = TRUE)
  }, character(1L), USE.NAMES = FALSE)

  workers <- min(requested_cores, length(vignettes))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message(
    "Staging vignette source and Quarto output in ", run_dir, ".\n",
    "Rendering ", length(vignettes), " vignette(s) with ", workers, " worker(s)."
  )

  cluster <- parallel::makePSOCKcluster(workers, outfile = "")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  results <- parallel::parLapply(cluster, vignettes, render_one)

  # Publish before reporting, and publish regardless of failures: a vignette
  # that took twenty minutes should not be discarded because another one broke.
  publish_files(vapply(results, `[[`, character(1L), "log_file"), output_dir)
  rendered <- Filter(Negate(is.na), lapply(results, `[[`, "html_file"))
  publish_files(unlist(rendered), output_dir)

  for (result in results) {
    status <- if (result$success) "OK" else if (!is.na(result$html_file)) "FAILED (html written)" else "FAILED"
    message(status, ": ", basename(result$vignette),
            " (log: ", file.path(output_dir, basename(result$log_file)), ")")
  }

  failures <- Filter(function(result) !result$success, results)
  if (length(failures)) {
    stop(length(failures), " vignette(s) failed; see their logs in ", output_dir, ".",
         call. = FALSE)
  }
  message("Rendered vignettes are in ", output_dir)
  invisible(results)
}

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores) || available_cores < 1L) {
  available_cores <- 1L
}

requested_cores <- if (!is.na(core_argument)) {
  suppressWarnings(as.integer(core_argument))
} else if (interactive()) {
  interactive_cores
} else {
  available_cores
}
if (!is.numeric(requested_cores) || length(requested_cores) != 1L ||
    is.na(requested_cores) || requested_cores < 1L ||
    requested_cores != as.integer(requested_cores)) {
  usage()
  stop("The core count must be a positive integer.", call. = FALSE)
}

render_vignettes(vignettes_dir, output_dir, requested_cores)
