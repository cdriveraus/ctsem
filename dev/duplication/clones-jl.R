# Same token-shingle clone detector, lexing Julia rather than R.
#
# No Julia needed: this is a regex lexer, which is enough because the detector
# only requires a token stream that is stable under renaming and blind to
# comments. Docstrings ("""...""") and # comments are stripped first -- in this
# codebase they are the bulk of the text and they would otherwise dominate.

lex_julia <- function(path) {
  src <- readLines(path, warn = FALSE, encoding = "UTF-8")
  n <- length(src)
  keep <- rep(TRUE, n)
  # Strip triple-quoted docstrings by line. Crude but safe here: the house
  # style opens and closes them on their own lines.
  indoc <- FALSE
  for (i in seq_len(n)) {
    q <- length(gregexpr('"""', src[i], fixed = TRUE)[[1]])
    has <- grepl('"""', src[i], fixed = TRUE)
    if (indoc) { keep[i] <- FALSE; if (has) indoc <- FALSE; next }
    if (has) {
      keep[i] <- FALSE
      if (q %% 2 == 1) indoc <- TRUE
    }
  }
  lines <- ifelse(keep, src, "")
  lines <- sub("#.*$", "", lines)          # line comments
  toks <- list(); lns <- integer(0)
  pat <- "[A-Za-z_][A-Za-z0-9_!]*|[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?|\"[^\"]*\"|[-+*/^%<>=!~&|:;,.()\\[\\]{}]"
  for (i in seq_along(lines)) {
    if (!nzchar(trimws(lines[i]))) next
    m <- regmatches(lines[i], gregexpr(pat, lines[i], perl = TRUE))[[1]]
    if (!length(m)) next
    toks[[length(toks) + 1L]] <- m
    lns <- c(lns, rep(i, length(m)))
  }
  tk <- unlist(toks)
  if (!length(tk)) return(NULL)
  # Normalise: numbers, strings and lowercase identifiers that are not Julia
  # keywords collapse to placeholders; keep keywords and anything with an
  # underscore or capital (function and type names carry the meaning).
  kw <- c("function","end","for","while","if","elseif","else","return","struct",
    "mutable","const","let","do","begin","local","global","in","where","import",
    "using","export","module","abstract","type","try","catch","finally","throw",
    "true","false","nothing","Int","Float64","Bool","Vector","Matrix")
  isnum <- grepl("^[0-9]", tk)
  isstr <- grepl('^"', tk)
  isid  <- grepl("^[A-Za-z_]", tk)
  out <- tk
  out[isnum] <- "N"; out[isstr] <- "S"
  plain <- isid & !isnum & !(tk %in% kw) & !grepl("[_A-Z]", tk)
  out[plain] <- "v"
  data.frame(tok = out, line = lns, stringsAsFactors = FALSE)
}

fun_index_jl <- function(path) {
  src <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hit <- grep("^\\s*function\\s+([A-Za-z_][A-Za-z0-9_!]*)", src)
  if (!length(hit)) return(NULL)
  nm <- sub("^\\s*function\\s+([A-Za-z_][A-Za-z0-9_!]*).*$", "\\1", src[hit])
  data.frame(name = nm, start = hit,
    end = c(hit[-1] - 1L, length(src)), stringsAsFactors = FALSE)
}

which_fun <- function(idx, line) {
  if (is.null(idx)) return(NA_character_)
  hit <- which(idx$start <= line & idx$end >= line)
  if (!length(hit)) NA_character_ else idx$name[hit[1]]
}

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else "inst/julia/ContinuousTimeSEM/src"
W   <- if (length(args) >= 2) as.integer(args[2]) else 60L
files <- list.files(dir, pattern = "\\.jl$", full.names = TRUE)

alltok <- list(); allfun <- list()
for (f in files) {
  n <- lex_julia(f)
  if (is.null(n) || nrow(n) < W) next
  n$file <- basename(f); alltok[[f]] <- n
  allfun[[basename(f)]] <- fun_index_jl(f)
}

keys <- new.env(hash = TRUE, parent = emptyenv())
for (f in names(alltok)) {
  n <- alltok[[f]]; tk <- n$tok
  m <- length(tk) - W + 1L
  if (m < 1) next
  for (i in seq_len(m)) {
    win <- tk[i:(i + W - 1L)]
    # A window of repeated literals (a captured data table, a long numeric
    # vector) normalises to a handful of distinct tokens and then matches
    # every other such window trivially. Require real variety before a window
    # can count as evidence of a copy.
    if (length(unique(win)) < 12L) next
    k <- paste(win, collapse = "")
    keys[[k]] <- rbind(keys[[k]], data.frame(file = n$file[i], tokstart = i,
      line1 = n$line[i], line2 = n$line[i + W - 1L], stringsAsFactors = FALSE))
  }
}

hits <- list()
for (k in ls(keys)) { d <- keys[[k]]; if (nrow(d) > 1) hits[[length(hits)+1L]] <- d }
cat("windows with >1 occurrence:", length(hits), "\n")
flat <- do.call(rbind, lapply(seq_along(hits), function(i) { d <- hits[[i]]; d$group <- i; d }))
if (is.null(flat)) { cat("no clones at window", W, "\n"); quit(save = "no") }

flat <- flat[order(flat$file, flat$tokstart), ]
run <- cumsum(c(1L, !(flat$file[-1] == flat$file[-nrow(flat)] &
  flat$tokstart[-1] == flat$tokstart[-nrow(flat)] + 1L)))
flat$run <- run
runs <- do.call(rbind, lapply(split(flat, flat$run), function(d) data.frame(
  file = d$file[1], line1 = min(d$line1), line2 = max(d$line2),
  ntok = W + nrow(d) - 1L, group = d$group[1], stringsAsFactors = FALSE)))
runs <- runs[order(-runs$ntok), ]
seen <- character(0); n_rep <- 0L
cat(sprintf("\n%-36s %-36s %6s\n", "location A", "location B", "tokens"))
cat(strrep("-", 82), "\n")
for (i in seq_len(nrow(runs))) {
  g <- runs$group[i]; mates <- runs[runs$group == g, ]
  if (nrow(mates) < 2) next
  a <- mates[1, ]; b <- mates[2, ]
  key <- paste(a$file, a$line1, b$file, b$line1)
  if (key %in% seen) next
  seen <- c(seen, key)
  cat(sprintf("%-36s %-36s %6d\n",
    sprintf("%s:%d-%d", a$file, a$line1, a$line2),
    sprintf("%s:%d-%d", b$file, b$line1, b$line2), a$ntok))
  cat(sprintf("   %-33s %-33s\n",
    paste0("in ", which_fun(allfun[[a$file]], a$line1)),
    paste0("in ", which_fun(allfun[[b$file]], b$line1))))
  n_rep <- n_rep + 1L
  if (n_rep >= 30L) break
}
cat("\nDONE, reported", n_rep, "clone pairs at window", W, "tokens\n")
