# Flag the habits of machine-written prose in what users read: help pages
# (roxygen in R/), vignettes, NEWS.md, NEW-IN-4.md and README.Rmd.
#
#   Rscript dev/check-prose.R                     # counts per file
#   Rscript dev/check-prose.R --lines             # and every flagged line
#   Rscript dev/check-prose.R --changed juliaFit  # only lines added since a commit
#   Rscript dev/check-prose.R --lines R/ctFit.R vignettes/uncertainty.qmd
#
# A hit is a place to reread, not an error: "not" is often the right word.
# Rewrite the sentence when it reads as written for effect, and leave it when
# a plain writer would have put it the same way. The script always exits 0.
# The rule it supports is "Writing for users" in CT-SEM/CLAUDE.md.

args <- commandArgs(trailingOnly = TRUE)
show_lines <- "--lines" %in% args
since <- NULL
if ("--changed" %in% args) {
  i <- match("--changed", args)
  since <- args[i + 1]
  args <- args[-c(i, i + 1)]
}
files <- setdiff(args, "--lines")

setwd(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")))

if (!length(files)) {
  files <- c("NEWS.md", "NEW-IN-4.md", "README.Rmd",
             Sys.glob("vignettes/*.qmd"), Sys.glob("R/*.R"))
}
files <- files[file.exists(files)]

# Each pattern: what it catches, and what to do instead.
patterns <- list(
  list(name = "X, not Y", re = ",\\s+not\\s+(?!only\\b)",
       fix = "say what it is; mention what it is not only if a reader would expect that"),
  list(name = "not X but Y", re = "\\bnot\\b[^.;:]{1,60}\\bbut\\b",
       fix = "state the positive case"),
  list(name = "rather than", re = "\\brather than\\b",
       fix = "drop the contrast unless the alternative is what a reader would assume"),
  list(name = "not just/only", re = "\\bnot (just|only|merely|simply)\\b",
       fix = "list the things plainly"),
  list(name = "never/always", re = "\\b(never|always)\\b",
       fix = "absolute claims: keep only where literally true and worth saying"),
  list(name = "dash aside", re = "(\\s--\\s|—|\\s–\\s)",
       fix = "use a comma, brackets or a new sentence"),
  list(name = "colon reveal", re = "\\b[a-z]{3,}:\\s+(it|this|that|they|the|a|an|so|which|each|every|one|no)\\b",
       fix = "write it as one ordinary sentence"),
  list(name = "signpost", re = "(?i)\\b(note that|importantly|crucially|in other words|put differently|in short|notably|essentially|simply put|in effect|that is,|in practice)\\b",
       fix = "delete the signpost"),
  list(name = "So/Where opener", re = "(^|[.!?]\\s+)(So|Where|Whatever)\\s+[a-z]",
       fix = "these openers are a house tic of the notes; start with the subject"),
  list(name = "whatever/itself", re = "\\b(whatever|by itself|on its own|in itself)\\b",
       fix = "usually padding"),
  list(name = "provenance", re = "(\\b20[0-9]{2}-[01][0-9]-[0-3][0-9]\\b|\\buntil (recently|version|now)\\b|\\b(it|this|that|they) used to\\b|\\bpreviously\\b|\\bjuliaFit\\b|review/|scripts/)",
       fix = "history and sources belong in commit messages and notes, not user docs"),
  list(name = "internal term", re = "(\\.ct[A-Z][A-Za-z]+\\(|\\b_[a-z][a-z0-9_]+!?|\\b(coreset|coresets|hybrid|anchors?|anchored|certif(y|ied|ication)|route|routes)\\b)",
       fix = "a user does not know this word; name what it does for them"),
  list(name = "semicolon", re = ";\\s+[a-z]",
       fix = "two sentences are usually plainer"))

# Lines a user reads, with their line numbers.
prose_lines <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  keep <- rep(FALSE, length(x))
  ext <- tolower(tools::file_ext(path))
  if (ext == "r") {
    rox <- grepl("^\\s*#'", x)
    body <- sub("^\\s*#'\\s?", "", x)
    # A block marked @noRd is a developer's note, not a help page.
    block <- cumsum(c(TRUE, rox[-1] != rox[-length(rox)]))
    nord <- tapply(grepl("^@noRd\\b", body) & rox, block, any)
    rox <- rox & !nord[as.character(block)]
    in_code <- FALSE
    for (i in seq_along(x)) {
      if (!rox[i]) { in_code <- FALSE; next }
      tag <- regmatches(body[i], regexpr("^@[A-Za-z]+", body[i]))
      if (length(tag)) {
        in_code <- tag %in% c("@examples", "@examplesIf", "@usage")
        if (tag %in% c("@importFrom", "@import", "@export", "@rdname",
                       "@aliases", "@keywords", "@name", "@family",
                       "@inheritParams", "@method", "@useDynLib", "@docType",
                       "@format", "@md", "@noRd", "@include", "@encoding",
                       "@concept", "@exportS3Method", "@examples",
                       "@examplesIf", "@usage")) next
      }
      if (!in_code) keep[i] <- TRUE
    }
    x <- body
  } else {
    in_fence <- FALSE; in_yaml <- FALSE; in_comment <- FALSE
    for (i in seq_along(x)) {
      l <- x[i]
      if (i == 1 && grepl("^---\\s*$", l) && ext %in% c("qmd", "rmd")) {
        in_yaml <- TRUE; next
      }
      if (in_yaml) { if (grepl("^---\\s*$", l)) in_yaml <- FALSE; next }
      if (grepl("^\\s*(```|~~~)", l)) { in_fence <- !in_fence; next }
      if (in_fence) next
      if (grepl("<!--", l)) in_comment <- TRUE
      if (in_comment) { if (grepl("-->", l)) in_comment <- FALSE; next }
      keep[i] <- TRUE
    }
  }
  # Code, links and maths are not prose.
  y <- x
  y <- gsub("`[^`]*`", " CODE ", y)
  y <- gsub("\\\\(code|link|eqn|deqn|href|url|emph|strong|item)\\{[^{}]*\\}", " CODE ", y)
  y <- gsub("\\$[^$]+\\$", " MATH ", y)
  y <- gsub("https?://\\S+", " URL ", y)
  data.frame(line = which(keep), text = y[keep], raw = x[keep],
             stringsAsFactors = FALSE)
}

added_lines <- function(path, rev) {
  d <- suppressWarnings(system2("git", c("diff", "-U0", rev, "--", shQuote(path)),
                                stdout = TRUE))
  hunks <- regmatches(d, regexpr("^@@ [^@]*\\+[0-9]+(,[0-9]+)? @@", d))
  out <- integer()
  for (h in hunks) {
    m <- regmatches(h, regexec("\\+([0-9]+)(,([0-9]+))?", h))[[1]]
    start <- as.integer(m[2]); n <- if (nzchar(m[4])) as.integer(m[4]) else 1L
    if (n > 0) out <- c(out, seq(start, length.out = n))
  }
  out
}

summary_rows <- list(); hits_all <- list()
for (f in files) {
  p <- prose_lines(f)
  if (!is.null(since)) p <- p[p$line %in% added_lines(f, since), , drop = FALSE]
  if (!nrow(p)) next
  words <- sum(lengths(strsplit(trimws(p$text), "\\s+")))
  counts <- integer(length(patterns))
  for (k in seq_along(patterns)) {
    hit <- grepl(patterns[[k]]$re, p$text, perl = TRUE)
    counts[k] <- sum(hit)
    if (any(hit)) hits_all[[length(hits_all) + 1]] <- data.frame(
      file = f, line = p$line[hit], pattern = patterns[[k]]$name,
      text = trimws(p$raw[hit]), stringsAsFactors = FALSE)
  }
  summary_rows[[f]] <- data.frame(file = f, words = words, hits = sum(counts),
    per1000 = round(1000 * sum(counts) / max(words, 1), 1),
    t(setNames(counts, vapply(patterns, `[[`, "", "name"))),
    check.names = FALSE, stringsAsFactors = FALSE)
}

if (!length(summary_rows)) {
  cat("No prose lines found in", length(files), "files.\n")
  quit(status = 0)
}
tab <- do.call(rbind, summary_rows)
tab <- tab[order(-tab$hits), ]
rownames(tab) <- NULL
tot <- colSums(tab[, -1])
cat(sprintf("%d files, %d prose words, %d flagged lines (%.1f per 1000 words)%s\n\n",
            nrow(tab), tot[["words"]], tot[["hits"]],
            1000 * tot[["hits"]] / max(tot[["words"]], 1),
            if (is.null(since)) "" else paste0(", lines added since ", since)))
cat("By pattern:\n")
for (k in seq_along(patterns)) {
  nm <- patterns[[k]]$name
  cat(sprintf("  %-15s %5d  %s\n", nm, tot[[nm]], patterns[[k]]$fix))
}
cat("\nBy file (most flagged first):\n")
print(tab[, c("file", "words", "hits", "per1000")], row.names = FALSE)

if (show_lines && length(hits_all)) {
  h <- do.call(rbind, hits_all)
  h <- h[order(h$file, h$line), ]
  cat("\n")
  for (i in seq_len(nrow(h)))
    cat(sprintf("%s:%d [%s] %s\n", h$file[i], h$line[i], h$pattern[i],
                substr(h$text[i], 1, 160)))
}
cat("\nDONE\n")
