#!/usr/bin/env Rscript

# Render every Quarto/R Markdown vignette in a temporary directory using a
# PSOCK cluster. Completed artifacts are then copied to dev/vignettes. In
# RStudio, open this file, change interactive_cores if needed, and click
# Source. From a terminal: Rscript dev/render-vignettes-parallel.R --cores 4

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
run_dir <- tempfile("ctsem-vignettes-", tmpdir = tempdir())
staged_vignettes_dir <- file.path(run_dir, "vignettes")
staged_output_dir <- file.path(run_dir, "output")

dir.create(staged_vignettes_dir, recursive = TRUE)
source_entries <- list.files(vignettes_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
if (!all(file.copy(source_entries, staged_vignettes_dir, recursive = TRUE, copy.date = TRUE))) {
  stop("Could not copy the vignette sources to the temporary render directory.", call. = FALSE)
}
on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)

vignettes <- list.files(
  staged_vignettes_dir,
  pattern = "\\.(qmd|rmd)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(vignettes)) {
  stop("No .qmd or .Rmd vignettes were found in the temporary vignette copy.", call. = FALSE)
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

workers <- min(requested_cores, length(vignettes))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staged_output_dir, recursive = TRUE)

publish_files <- function(files, output_dir) {
  if (!length(files)) {
    return(invisible(NULL))
  }
  if (!all(file.copy(files, output_dir, recursive = TRUE, overwrite = TRUE))) {
    stop("Could not copy rendered vignette artifacts to ", output_dir, ".", call. = FALSE)
  }
}

render_one <- function(vignette, output_dir) {
  name <- tools::file_path_sans_ext(basename(vignette))
  log_file <- file.path(output_dir, paste0(name, ".log"))
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
        quarto_args = c("--output-dir", output_dir),
        as_job = FALSE
      )
      list(vignette = vignette, success = TRUE, log_file = log_file, error = NA_character_)
    },
    error = function(error) {
      message(conditionMessage(error))
      list(vignette = vignette, success = FALSE, log_file = log_file, error = conditionMessage(error))
    }
  )
  message("Finished ", basename(vignette), " at ", format(Sys.time(), tz = "UTC"))
  result
}

message(
  "Staging vignette source and Quarto output in ", run_dir, ".\n",
  "Rendering ", length(vignettes), " vignette(s) with ", workers, " worker(s)."
)
cluster <- parallel::makePSOCKcluster(workers, outfile = "")
on.exit(parallel::stopCluster(cluster), add = TRUE)
results <- parallel::parLapply(cluster, vignettes, render_one, output_dir = staged_output_dir)

failures <- Filter(function(result) !result$success, results)
publish_files(list.files(staged_output_dir, pattern = "\\.log$", full.names = TRUE), output_dir)
for (result in results) {
  status <- if (result$success) "OK" else "FAILED"
  message(status, ": ", basename(result$vignette), " (log: ", file.path(output_dir, basename(result$log_file)), ")")
}

if (length(failures)) {
  stop(length(failures), " vignette(s) failed; see their logs in ", output_dir, ".", call. = FALSE)
}

publish_files(list.files(staged_output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE), output_dir)
message("Rendered vignettes are in ", output_dir)
