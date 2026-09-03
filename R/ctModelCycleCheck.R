# Circular dependencies between t0 matrices, PARS and the latent states.
#
# The initial state is set from T0MEANS, PARS is computed from the state, and
# any system matrix may reference either. That makes a loop expressible in the
# model specification, and the generated Stan program resolves it by evaluation
# order rather than by a rule: T0MEANS is built once from the cells that depend
# on nothing, the state is seeded from that, PARS is computed, then the t0
# matrices are rebuilt to pick up their PARS-dependent cells. A cell caught in a
# genuine loop is read before anything has written it, and what it reads is a
# number rather than an error.
#
# So the loop has to be refused at specification time, where the user can still
# see which cells they wrote. Two rules, and the messages name the cells.

# Latent state indices referenced by a cell's text, by name or as state[k].
# Names are matched longest-first so eta1 cannot match inside eta10, and on a
# word boundary so it cannot match inside a longer identifier.
.ctCycleStateRefs <- function(txt, latentNames) {
  if (is.null(txt) || is.na(txt) || !nzchar(txt)) return(integer(0))
  txt <- strsplit(txt, "|", fixed = TRUE)[[1]][1]
  if (is.na(txt)) return(integer(0))
  out <- integer(0)
  bracketed <- regmatches(txt, gregexpr("\\bstate\\s*\\[\\s*\\d+\\s*\\]", txt))[[1]]
  if (length(bracketed)) {
    out <- c(out, as.integer(gsub("\\D", "", bracketed)))
  }
  if (length(latentNames)) {
    ord <- order(nchar(latentNames), decreasing = TRUE)
    for (i in ord) {
      nm <- latentNames[i]
      if (!nzchar(nm)) next
      if (grepl(paste0("\\b", nm, "\\b"), txt)) out <- c(out, i)
    }
  }
  sort(unique(out))
}

# PARS cells referenced by a cell's text, as "row,col" keys.
.ctCycleParsRefs <- function(txt) {
  if (is.null(txt) || is.na(txt) || !nzchar(txt)) return(character(0))
  txt <- strsplit(txt, "|", fixed = TRUE)[[1]][1]
  if (is.na(txt)) return(character(0))
  m <- regmatches(txt, gregexpr("\\bPARS\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]", txt))[[1]]
  if (!length(m)) return(character(0))
  unique(vapply(m, function(x) {
    d <- as.integer(regmatches(x, gregexpr("\\d+", x))[[1]])
    paste0(d[1], ",", d[2])
  }, character(1), USE.NAMES = FALSE))
}

.ctCycleCellName <- function(matrix, row, col) {
  paste0(matrix, "[", row, ",", col, "]")
}

# Rule 1. A t0 matrix defines the initial state, so a state reference in one is
# circular by construction: the cell's value would have to be known before the
# state it reads exists. Today such a cell is not refused but silently
# reinterpreted, because it is stored with when = 0 and its state index lands in
# the column that otherwise holds a parameter number, so `mcalc` materialises it
# from the parameter vector instead of from the state.
.ctCheckT0StateRefs <- function(pars, latentNames, t0matrices) {
  rows <- which(pars$matrix %in% t0matrices & !is.na(pars$param))
  bad <- character(0)
  for (ri in rows) {
    refs <- .ctCycleStateRefs(pars$param[ri], latentNames)
    if (length(refs)) {
      bad <- c(bad, paste0(
        .ctCycleCellName(pars$matrix[ri], pars$row[ri], pars$col[ri]),
        " references ", paste0(latentNames[refs], collapse = " and ")))
    }
  }
  if (length(bad)) {
    stop("A t0 matrix cannot reference a latent state:\n  ",
      paste0(bad, collapse = "\n  "),
      "\nThe t0 matrices set the initial state, so a state reference in one is ",
      "circular: the cell would have to be evaluated before the state it reads ",
      "exists. Put the state dependent quantity in PARS and reference that from ",
      "the matrix that needs it, or give this cell a value that does not depend ",
      "on a state.", call. = FALSE)
  }
  invisible(TRUE)
}

# Rule 2. A loop that runs through PARS. state k comes from T0MEANS[k,1], so
# state k depends on whatever that cell references; a PARS cell depends on
# whatever its own text references. A cycle in that graph cannot be evaluated in
# any order.
.ctCheckStatePARSCycles <- function(pars, latentNames, nlatent) {
  celltext <- function(matrix, row, col) {
    ri <- which(pars$matrix %in% matrix & pars$row == row & pars$col == col)
    if (!length(ri)) return(NA_character_)
    pars$param[ri[1]]
  }

  # Edges out of each node, as node keys.
  edges <- list()
  refsToNodes <- function(txt) {
    c(paste0("state:", .ctCycleStateRefs(txt, latentNames)),
      paste0("pars:", .ctCycleParsRefs(txt)))
  }
  if (nlatent >= 1) for (k in seq_len(nlatent)) {
    edges[[paste0("state:", k)]] <- refsToNodes(celltext("T0MEANS", k, 1))
  }
  parsrows <- which(pars$matrix %in% "PARS" & !is.na(pars$param))
  for (ri in parsrows) {
    key <- paste0("pars:", pars$row[ri], ",", pars$col[ri])
    edges[[key]] <- refsToNodes(pars$param[ri])
  }

  label <- function(node) {
    if (startsWith(node, "state:")) {
      k <- as.integer(sub("state:", "", node, fixed = TRUE))
      nm <- if (k >= 1 && k <= length(latentNames)) latentNames[k] else paste0("state ", k)
      paste0("state ", nm, " (from T0MEANS[", k, ",1])")
    } else {
      paste0("PARS[", sub("pars:", "", node, fixed = TRUE), "]")
    }
  }

  # Depth first search reporting the first cycle it closes, as a path.
  state <- setNames(rep("new", length(edges)), names(edges))
  found <- NULL
  visit <- function(node, path) {
    if (!is.null(found)) return(invisible(NULL))
    if (!node %in% names(edges)) return(invisible(NULL))
    if (identical(state[[node]], "open")) {
      start <- match(node, path)
      found <<- c(path[start:length(path)], node)
      return(invisible(NULL))
    }
    if (identical(state[[node]], "closed")) return(invisible(NULL))
    state[[node]] <<- "open"
    for (nxt in edges[[node]]) visit(nxt, c(path, node))
    state[[node]] <<- "closed"
    invisible(NULL)
  }
  for (node in names(edges)) visit(node, character(0))

  if (!is.null(found)) {
    stop("Circular dependency in the model specification:\n  ",
      paste0(vapply(found, label, character(1)), collapse = "\n    -> "),
      "\nThe initial state is set from T0MEANS and PARS is computed from the ",
      "state, so this chain has no order in which every cell can be evaluated. ",
      "Break it by making one of these cells a fixed value, or independent of ",
      "the others.", call. = FALSE)
  }
  invisible(TRUE)
}

# Both rules. Called from the model specification, on the finalised par table.
.ctCheckModelCycles <- function(pars, latentNames, nlatent,
  t0matrices = c("T0MEANS", "T0VAR", "J0")) {
  if (is.null(pars) || !nrow(pars)) return(invisible(TRUE))
  if (!all(c("matrix", "row", "col", "param") %in% colnames(pars))) {
    return(invisible(TRUE))
  }
  t0matrices <- intersect(t0matrices, unique(pars$matrix))
  .ctCheckT0StateRefs(pars, latentNames, t0matrices)
  .ctCheckStatePARSCycles(pars, latentNames, nlatent)
  invisible(TRUE)
}
