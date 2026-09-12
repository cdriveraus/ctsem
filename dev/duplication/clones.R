# Token-shingle clone detector for R sources. No dependencies beyond base R.
#
# Normalises away comments, whitespace, literals and local variable names, so
# that a copied block still matches after a rename. Keeps function-call names,
# which is what stops it reporting every `for` loop as a clone of every other.

normalise_file <- function(path) {
  exprs <- try(parse(path, keep.source = TRUE), silent = TRUE)
  if (inherits(exprs, "try-error")) return(NULL)
  pd <- utils::getParseData(exprs)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  pd <- pd[pd$terminal & pd$token != "COMMENT", ]
  pd <- pd[order(pd$line1, pd$col1), ]
  tok <- pd$token
  txt <- pd$text
  out <- ifelse(tok %in% c("SYMBOL", "SLOT"), "v",
    ifelse(tok == "NUM_CONST", "N",
      ifelse(tok == "STR_CONST", "S",
        ifelse(tok %in% c("SYMBOL_FORMALS", "SYMBOL_SUB"), "a", txt))))
  data.frame(tok = out, line = pd$line1, stringsAsFactors = FALSE)
}

# Top-level function definitions, so a clone can be reported as "in f and g".
fun_index <- function(path) {
  exprs <- try(parse(path, keep.source = TRUE), silent = TRUE)
  if (inherits(exprs, "try-error")) return(NULL)
  refs <- attr(exprs, "srcref")
  nm <- character(0); s <- integer(0); e <- integer(0)
  for (i in seq_along(exprs)) {
    ex <- exprs[[i]]
    if (is.call(ex) && length(ex) >= 3 &&
        as.character(ex[[1]])[1] %in% c("<-", "=", "assign") &&
        is.call(ex[[3]]) && as.character(ex[[3]][[1]])[1] == "function") {
      nm <- c(nm, paste(deparse(ex[[2]]), collapse = ""))
      s <- c(s, refs[[i]][1]); e <- c(e, refs[[i]][3])
    }
  }
  if (!length(nm)) return(NULL)
  data.frame(name = nm, start = s, end = e, stringsAsFactors = FALSE)
}

which_fun <- function(idx, line) {
  if (is.null(idx)) return(NA_character_)
  hit <- which(idx$start <= line & idx$end >= line)
  if (!length(hit)) return(NA_character_) else idx$name[hit[1]]
}

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else "R"
W   <- if (length(args) >= 2) as.integer(args[2]) else 60L   # window, in tokens
files <- list.files(dir, pattern = "\\.[Rr]$", full.names = TRUE)

alltok <- list(); allfun <- list()
for (f in files) {
  n <- normalise_file(f)
  if (is.null(n) || nrow(n) < W) next
  n$file <- basename(f)
  alltok[[f]] <- n
  allfun[[basename(f)]] <- fun_index(f)
}

# Hash every window of W tokens.
keys <- new.env(hash = TRUE, parent = emptyenv())
for (f in names(alltok)) {
  n <- alltok[[f]]
  tk <- n$tok
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
    prev <- keys[[k]]
    keys[[k]] <- rbind(prev, data.frame(file = n$file[i], tokstart = i,
      line1 = n$line[i], line2 = n$line[i + W - 1L], stringsAsFactors = FALSE))
  }
}

# Keep windows seen in more than one place, then merge overlapping windows
# into maximal runs so one copied block is reported once, not W times.
hits <- list()
for (k in ls(keys)) {
  d <- keys[[k]]
  if (nrow(d) > 1) hits[[length(hits) + 1L]] <- d
}
cat("windows with >1 occurrence:", length(hits), "\n")

flat <- do.call(rbind, lapply(seq_along(hits), function(i) {
  d <- hits[[i]]; d$group <- i; d }))
if (is.null(flat)) { cat("no clones at window", W, "\n"); quit(save = "no") }

# Merge: within a file, consecutive matching windows (tokstart differing by 1)
# that belong to groups of the same size are one run.
flat <- flat[order(flat$file, flat$tokstart), ]
run <- cumsum(c(1L, !(flat$file[-1] == flat$file[-nrow(flat)] &
  flat$tokstart[-1] == flat$tokstart[-nrow(flat)] + 1L)))
flat$run <- run
runs <- do.call(rbind, lapply(split(flat, flat$run), function(d) data.frame(
  file = d$file[1], line1 = min(d$line1), line2 = max(d$line2),
  ntok = W + nrow(d) - 1L, group = d$group[1], stringsAsFactors = FALSE)))

# Pair runs that share a group -- i.e. runs that are copies of each other.
runs <- runs[order(-runs$ntok), ]
seen <- character(0)
cat(sprintf("\n%-34s %-34s %6s\n", "location A", "location B", "tokens"))
cat(strrep("-", 78), "\n")
n_rep <- 0L
for (i in seq_len(nrow(runs))) {
  g <- runs$group[i]
  mates <- runs[runs$group == g, ]
  if (nrow(mates) < 2) next
  a <- mates[1, ]; b <- mates[2, ]
  key <- paste(a$file, a$line1, b$file, b$line1)
  if (key %in% seen) next
  seen <- c(seen, key)
  fa <- which_fun(allfun[[a$file]], a$line1)
  fb <- which_fun(allfun[[b$file]], b$line1)
  cat(sprintf("%-34s %-34s %6d\n",
    sprintf("%s:%d-%d", a$file, a$line1, a$line2),
    sprintf("%s:%d-%d", b$file, b$line1, b$line2), a$ntok))
  cat(sprintf("   %-31s %-31s\n", paste0("in ", fa), paste0("in ", fb)))
  n_rep <- n_rep + 1L
  if (n_rep >= 40L) break
}
cat("\nDONE, reported", n_rep, "clone pairs at window", W, "tokens\n")
