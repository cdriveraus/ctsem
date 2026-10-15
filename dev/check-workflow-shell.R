# Parse every `run:` script in .github/workflows with `bash -n` before pushing.
#
# The R and Julia those steps run sit inside single-quoted shell strings
# (`Rscript -e '...'`, `julia -e '...'`), so one apostrophe anywhere in them --
# a comment's "test's" -- ends the string, and the step fails on GitHub after
# its whole suite has run: all four julia-tests shards did, on juliaFit
# feba1657. Parsing the YAML alone does not see it; parsing each script as the
# runner's bash does, does.
#
#   Rscript dev/check-workflow-shell.R        # from the package root
files <- list.files(".github/workflows", pattern = "\\.ya?ml$", full.names = TRUE)
bad <- 0L
for (f in files) {
  jobs <- yaml::read_yaml(f)$jobs
  k <- 0L
  for (job in names(jobs)) for (step in jobs[[job]]$steps) {
    if (is.null(step$run) || identical(step$shell, "pwsh")) next
    k <- k + 1L
    script <- tempfile(fileext = ".sh")
    writeLines(step$run, script)
    out <- suppressWarnings(system2("bash", c("-n", script), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
      bad <- bad + 1L
      cat(sprintf("%s, job %s, step %s:\n  %s\n", basename(f), job,
        if (is.null(step$name)) k else step$name,
        paste(sub("^[^:]*: ", "", out), collapse = "\n  ")))
    }
    unlink(script)
  }
}
if (bad) stop(bad, " workflow step(s) do not parse as shell", call. = FALSE)
cat("every workflow step parses\n")
