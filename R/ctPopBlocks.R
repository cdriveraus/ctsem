# Shared dimensions in blocks, on a reduced-rank level: `ctFit(..., popblocks =)`.
#
# Under `poprank = r` every effect at a level loads on all r shared dimensions,
# so any effect may covary with any other. `popblocks` says which effects
# share which dimensions. Each block is a number of dimensions and the effects
# that load on them; an effect loads on the blocks that name it and on no
# others. Two effects covary only through the blocks they share, so effects
# with no block in common are uncorrelated. An effect named in no block has no
# shared dimension at all and needs an own variance (`popunique`).
#
# Why it exists: in the affect pilot each item's person-level intercept
# residual sees the latent level only through the item's loading. If the
# residuals could load on the level's dimensions, they could copy any change
# in the level's loadings, and the level's variance would not be identified.
# Blocks keep the residuals off the level's dimensions while letting them
# covary with each other, and with the dynamics, on dimensions of their own.
# Zero covariance with the level does not by itself identify the level: the
# residuals' own dimensions must not be able to reproduce a pattern in
# proportion to the item loadings, which a low rank for their block keeps
# them from doing generically.
#
# In the engine a patterned level's loadings carry the effect and dimension of
# each entry (`load_row`, `load_col`, `_laplace_loading_entries`). Each block
# is triangular in its own first effects, which removes the rotation within
# it; the zeros between blocks remove the rest. A model without blocks keeps
# the triangle and sends no pattern, so it crosses to the engine exactly as
# before. Like `popunique`, this is built where the engine builds loadings, so
# it reaches `'laplace'` and `'none'` and nothing else.

# TRUE when `popblocks` asks for anything.
.ctPopblocksRequested <- function(popblocks) {
  !is.null(popblocks) && length(popblocks) > 0L
}

# One block as given: a list with `rank` and `effects`.
.ctPopblocksIsBlock <- function(x) is.list(x) && !is.null(x[["rank"]])

# `popblocks` per grouping level. A list of blocks applies to every level, as
# an unnamed `poprank` or `popunique` does; a list named by level gives each
# level its own.
.ctPopblocksByLevel <- function(popblocks, levelnames) {
  if (!is.list(popblocks)) stop("popblocks must be a list of blocks, each a ",
    "list with rank and effects, or such lists named by level.", call. = FALSE)
  if (all(vapply(popblocks, .ctPopblocksIsBlock, logical(1)))) {
    return(stats::setNames(rep(list(popblocks), length(levelnames)), levelnames))
  }
  named <- names(popblocks)
  if (is.null(named) || any(is.na(named) | !nzchar(named))) stop(
    "popblocks must be a list of blocks, each a list with rank and effects, ",
    "or such lists named by level: ", paste(levelnames, collapse = ', '), ".",
    call. = FALSE)
  unknown <- setdiff(named, levelnames)
  if (length(unknown)) stop("popblocks names no level called ",
    paste(unknown, collapse = ', '), ". The levels are ",
    paste(levelnames, collapse = ', '), ".", call. = FALSE)
  for (level in named) {
    if (!all(vapply(popblocks[[level]], .ctPopblocksIsBlock, logical(1)))) stop(
      "popblocks for level '", level, "' must be a list of blocks, each a list ",
      "with rank and effects.", call. = FALSE)
  }
  popblocks
}

# Resolve and check `popblocks`: a named list, level -> its blocks, each with
# `name`, `rank` and `param`, the parameters that load on it. `effects` selects
# as `popunique` does: TRUE for every effect varying at the level, or names of
# matrices or parameters.
.ctPopblocksResolve <- function(ctm, popblocks, levelnames, levelcolumns) {
  asked <- .ctPopblocksByLevel(popblocks, levelnames)
  out <- list()
  for (level in names(asked)) {
    blocks <- asked[[level]]
    labels <- names(blocks)
    resolved <- list()
    for (b in seq_along(blocks)) {
      label <- if (!is.null(labels) && nzchar(labels[b])) labels[b] else
        paste0("block", b)
      rank <- suppressWarnings(as.numeric(blocks[[b]][["rank"]]))
      if (length(rank) != 1L || !is.finite(rank) || rank < 1 || rank != round(rank)) {
        stop("popblocks: block '", label, "' at level '", level, "' needs a ",
          "whole rank of 1 or more.", call. = FALSE)
      }
      effects <- blocks[[b]][["effects"]]
      if (is.null(effects)) stop("popblocks: block '", label, "' at level '",
        level, "' names no effects.", call. = FALSE)
      param <- .ctPopuniqueSelect(ctm$pars, levelcolumns[[level]], effects,
        level, what = paste0("popblocks block '", label, "'"))
      if (length(param) < rank) stop("popblocks: block '", label, "' at level '",
        level, "' has rank ", rank, " over ", length(param), " effect",
        if (length(param) != 1L) "s", ": its first effects form its basis, so ",
        "the rank can be at most their number.", call. = FALSE)
      resolved[[b]] <- list(name = label, rank = as.integer(rank), param = param)
    }
    sets <- vapply(resolved, function(x) paste(sort(x$param), collapse = "\r"),
      character(1))
    if (anyDuplicated(sets)) stop("popblocks: two blocks at level '", level,
      "' name the same effects, so rotations between them are not ",
      "identified. Give them as one block with the summed rank.", call. = FALSE)
    out[[level]] <- resolved
  }
  out
}

# A level's loading pattern from its blocks: the effect (`row`, an index into
# `names`, the level's varying parameters in order) and dimension (`col`) of
# each loading, dimension by dimension, and each loading's prior scale. Block
# `b` takes the next `rank` dimensions; its effects in level order, the first
# `rank` of them its basis: dimension `j` of the block is loaded by its effects
# from the `j`-th on. With one block naming every effect that is the triangle
# a level without blocks has, in the same order. The prior scale is
# `1/sqrt(n)`, `n` the dimensions the loading's effect may load on, which keeps
# an effect's implied prior variance at 1 as `1/sqrt(rank)` does for the
# triangle.
.ctJuliaBlockPattern <- function(blocks, names) {
  rows <- integer(); cols <- integer()
  span <- numeric(length(names))
  for (b in blocks) span[names %in% b$param] <- span[names %in% b$param] + b$rank
  start <- 0L
  for (b in blocks) {
    members <- which(names %in% b$param)
    if (length(members) < b$rank) stop("Internal error: block '", b$name,
      "' has fewer varying effects at this level than its rank.", call. = FALSE)
    for (j in seq_len(b$rank)) {
      take <- members[seq_along(members) >= j]
      rows <- c(rows, take)
      cols <- c(cols, rep(start + j, length(take)))
    }
    start <- start + b$rank
  }
  list(row = as.integer(rows), col = as.integer(cols),
    prior = 1 / sqrt(span[rows]))
}

# The effect and dimension of each loading of a reduced-rank level of a
# specification, in `load_index` order: its pattern, or the triangle a level
# without one has. The R side's one reader of the layout, as
# `_laplace_loading_entries` is the engine's.
.ctJuliaLoadingEntries <- function(level) {
  if (length(level$load_row)) {
    return(list(row = as.integer(level$load_row), col = as.integer(level$load_col)))
  }
  k <- as.integer(level$nrandom)
  rank <- as.integer(if (is.null(level$rank)) k else level$rank)
  rows <- integer(); cols <- integer()
  for (q in seq_len(rank)) {
    rows <- c(rows, seq(q, k))
    cols <- c(cols, rep(q, k - q + 1L))
  }
  n <- length(level$load_index)
  list(row = rows[seq_len(n)], col = cols[seq_len(n)])
}

# What the summary says about a level's blocks, or NULL when none has any.
.ctBackendPopBlocksNote <- function(spec) {
  levels <- spec$laplace$levels
  out <- character()
  for (lv in levels) {
    if (!length(lv$blocks)) next
    described <- vapply(lv$blocks, function(b) paste0(b$name, " (",
      b$rank, " dimension", if (b$rank > 1L) "s", ": ",
      paste(b$param, collapse = ", "), ")"), character(1))
    alone <- setdiff(as.character(lv$param),
      unlist(lapply(lv$blocks, function(b) b$param)))
    out <- c(out, paste0("Shared dimensions at level '", lv$name, "' come in ",
      "blocks: ", paste(described, collapse = "; "), ". Effects with no block ",
      "in common are uncorrelated",
      if (length(alone)) paste0(", and ", paste(alone, collapse = ", "),
        if (length(alone) > 1L) " have" else " has",
        " an own variance alone"), "."))
  }
  if (length(out)) paste(out, collapse = " ") else NULL
}
