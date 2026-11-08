# The bound on correlations: random effects', and the model's own covariance
# matrices'.
#
# Unpenalised maximum likelihood often puts a random-effects covariance on the
# singular boundary: an effect with no variance of its own, which the others
# predict exactly. The sd-and-correlation constructions reach that boundary
# only as their coordinates go to infinity, so without a prior a fit walked
# along a ridge, stopped wherever its stopping rule said (the augmented route
# 0.56 nats short of the laplace one on the ?ctFit model), and reported the
# boundary as non-identification. A correlation matrix can go singular without
# any one correlation nearing one -- that model's maximum-likelihood covariance
# had rank 2 of 4 and no correlation above 0.96 -- so the bound is on each
# effect's multiple correlation with the others at its level: at most 0.999,
# which for two effects is the correlation itself. The engine holds it with a
# wall in the objective (`_ctsem_corbound`, ctsem_backend.jl), on every julia
# route; review/PRIORS-default-boundary-2026-10-19.md has the evidence.
#
# T0VAR, DIFFUSION and MANIFESTVAR are built from the same correlation
# constructions and reach the same boundary the same way. On the ctsem book's
# linear trend model with random intercepts, T0VAR's correlation walked to
# one (raw -8.7 to -10) once the prior warm-up placed the start, and all eight
# seeds ended unconverged within 0.013 nats of the maximum, where Laplace's
# term jumps; bounded, all eight converged at the bound, 0.0004 nats below
# the unbounded maximum (2026-11-08, review/bench/2026-11-08-warmup.md). So
# their correlations are bounded too, since 2026-11-08.
#
# This file says which blocks the wall applies to, and reports a fit that
# ended on it.

# The blocks the bound applies to, as the engine reads them: one per full-rank
# population covariance of two or more effects built from correlations (`'z'`
# or the correlation square root), with each cell below the diagonal, in
# column-major order, as the raw index it reads, or 0 and the value it is fixed
# at. NULL when there are none. A reduced-rank level is singular by design and
# built from loadings, and `'cholesky'` has no correlation cells, so neither is
# bounded.
#
# Both routes read the coordinates untransformed -- a Laplace level as its
# `cor_index`, the augmented route as `RAWPOPVAR`'s cells -- which is what lets
# one description serve both; an augmented cell with any other transform is not
# this construction's, and its block is left unbounded rather than guessed at.
#
# The population's blocks come first, then the model's own matrices' (each
# with `matrix` naming it), so a reader of the population's takes those whose
# `matrix` is NULL.
.ctCorBoundSpec <- function(spec,
    table = as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)) {
  code <- if (is.null(spec$covmatcode)) 0L else as.integer(spec$covmatcode)[1L]
  if (!code %in% c(0L, 2L)) return(NULL)
  blocks <- c(.ctCorBoundPopulation(spec, table, code),
    .ctCorBoundMatrices(table, code))
  if (length(blocks)) blocks else NULL
}

# A k x k matrix's cells below the diagonal in column-major order, as the
# engine reads a block's coordinates.
.ctCorBoundLower <- function(k) {
  cells <- which(lower.tri(diag(k)), arr.ind = TRUE)
  cells[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
}

# The population covariance's blocks: one per full-rank level, as above.
.ctCorBoundPopulation <- function(spec, table, code) {
  lower <- .ctCorBoundLower
  if (!is.null(spec$laplace)) {
    blocks <- list()
    for (lv in spec$laplace$levels) {
      k <- as.integer(lv$nrandom)
      rank <- if (is.null(lv$rank)) k else as.integer(lv$rank)
      index <- as.integer(lv$cor_index)
      if (length(k) != 1L || k < 2L || rank < k || length(lv$load_index) ||
          length(index) != k * (k - 1L) / 2L) next
      fixed <- if (length(lv$cor_fixed) == length(index)) as.numeric(lv$cor_fixed) else
        rep(NA_real_, length(index))
      blocks[[length(blocks) + 1L]] <- list(k = k, code = code, index = index,
        fixed = fixed, names = as.character(lv$param),
        level = if (is.null(lv$name)) "" else as.character(lv$name))
    }
    return(blocks)
  }
  if (!any(table$matrix %in% "RAWPOPVAR")) return(list())
  if (!is.null(spec$model) && .ctPopFactorConstruction(spec$model)) return(list())
  pop <- table[table$matrix %in% "RAWPOPVAR", , drop = FALSE]
  row <- as.integer(pop$row); col <- as.integer(pop$col)
  k <- max(row)
  if (k < 2L) return(list())
  cells <- lower(k)
  index <- integer(nrow(cells)); fixed <- rep(NA_real_, nrow(cells))
  for (r in seq_len(nrow(cells))) {
    hit <- which(row == cells[r, 1L] & col == cells[r, 2L])
    if (length(hit) != 1L) return(list())
    p <- suppressWarnings(as.integer(pop$parnumber[hit]))
    if (!is.na(p) && p > 0L) {
      if (!identical(gsub(" ", "", as.character(pop$transform[hit])),
        paste0("param[", p, "]"))) return(list())
      index[r] <- p
    } else {
      fixed[r] <- as.numeric(pop$value[hit])
    }
  }
  diagonal <- which(row == col)
  names <- sub("^popsd_", "", as.character(pop$param[diagonal][order(row[diagonal])]))
  list(list(k = k, code = code, index = index, fixed = fixed, names = names,
    level = if (is.null(spec$subjectIDname)) "" else as.character(spec$subjectIDname)[1L]))
}

# The model's own covariance matrices' blocks: for each of T0VAR, DIFFUSION
# and MANIFESTVAR with a free correlation, the states its free correlations
# join, and any a fixed nonzero correlation joins to those, since the
# construction couples them. Read from the parameter table, which lays the
# cells out alike on every route. A free correlation under any transform but
# the identity is not the construction's coordinate, and its matrix is left
# unbounded rather than guessed at, as the population's is.
.ctCorBoundMatrices <- function(table, code) {
  blocks <- list()
  for (m in c("T0VAR", "DIFFUSION", "MANIFESTVAR")) {
    rows <- table[table$matrix %in% m, , drop = FALSE]
    if (!nrow(rows)) next
    row <- as.integer(rows$row); col <- as.integer(rows$col)
    p <- suppressWarnings(as.integer(rows$parnumber))
    p[is.na(p)] <- 0L
    value <- suppressWarnings(as.numeric(rows$value))
    # Below the diagonal: ctsem reads a covariance matrix's lower triangle.
    off <- row > col
    free <- off & p > 0L
    if (!any(free)) next
    identity <- gsub(" ", "", as.character(rows$transform[free])) ==
      paste0("param[", p[free], "]")
    if (!all(identity %in% TRUE)) next
    linked <- off & (free | (!is.na(value) & value != 0))
    states <- sort(unique(c(row[free], col[free])))
    repeat {
      more <- sort(unique(c(states, row[linked & col %in% states],
        col[linked & row %in% states])))
      if (length(more) == length(states)) break
      states <- more
    }
    k <- length(states)
    if (k < 2L) next
    cells <- .ctCorBoundLower(k)
    index <- integer(nrow(cells)); fixed <- rep(0, nrow(cells))
    for (r in seq_len(nrow(cells))) {
      a <- states[cells[r, 1L]]; b <- states[cells[r, 2L]]
      hit <- which(off & row == a & col == b)
      if (!length(hit)) next
      hit <- hit[1L]
      if (p[hit] > 0L) index[r] <- p[hit] else
        if (!is.na(value[hit])) fixed[r] <- value[hit]
    }
    names <- vapply(states, function(i) {
      d <- which(row == i & col == i)
      if (length(d) && !is.na(rows$param[d[1L]]) && nzchar(rows$param[d[1L]]))
        as.character(rows$param[d[1L]]) else paste0(m, "[", i, ",", i, "]")
    }, character(1))
    blocks[[length(blocks) + 1L]] <- list(k = k, code = code, index = index,
      fixed = fixed, names = names, level = "", matrix = m)
  }
  blocks
}

# The bound, read from the engine rather than duplicated.
.ctCorBoundValue <- function() {
  value <- try(.ctJuliaEval("ContinuousTimeSEM._CTSEM_CORBOUND"), silent = TRUE)
  if (inherits(value, "try-error") || !length(value)) 0.999 else as.numeric(value)[1L]
}

# A bounded block's correlation matrix at raw parameters `raw`, as the engine
# builds it: its cells below the diagonal in column-major order, each the raw
# coordinate its index reads or the value it is fixed at.
.ctCorBoundMatrix <- function(b, raw) {
  cells <- diag(b$k)
  low <- which(lower.tri(cells), arr.ind = TRUE)
  low <- low[order(low[, 2L], low[, 1L]), , drop = FALSE]
  cells[low] <- ifelse(b$index > 0L, raw[pmax(b$index, 1L)], b$fixed)
  .ctCovForward(cells, if (b$code == 2L) "z" else "rawcorr")
}

# Each bounded effect's multiple correlation with the others at the fit's
# estimate, whether it is on the bound, and the rank the block's correlation
# matrix points to (`.ctImproveRankFromCov`); NULL when nothing is bounded.
.ctCorBoundReport <- function(fit) {
  blocks <- try(.ctCorBoundSpec(fit$model_spec), silent = TRUE)
  raw <- as.numeric(fit$estimate$raw)
  if (inherits(blocks, "try-error") || !length(blocks) || !length(raw)) return(NULL)
  bound <- .ctCorBoundValue()
  rows <- lapply(blocks, function(b) {
    R <- .ctCorBoundMatrix(b, raw)
    inverse <- try(solve(R), silent = TRUE)
    mc <- if (inherits(inverse, "try-error")) rep(1, b$k) else
      sqrt(pmax(0, 1 - 1 / (diag(R) * diag(inverse))))
    # The wall lets an estimate sit a hair past the bound; within 5e-4 of it
    # is on it.
    data.frame(matrix = if (is.null(b$matrix)) "population" else b$matrix,
      level = b$level, effect = b$names, multiple_correlation = mc,
      at_bound = mc >= bound - 5e-4, rank = .ctImproveRankFromCov(R), k = b$k,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "bound") <- bound
  out
}

# The rank each bounded block's correlation matrix points to across a
# sample's draws, as `.ctCorBoundReport()` reads it at an estimate: the median
# over up to `maxdraws` draws, evenly spaced, of `.ctImproveRankFromCov()`,
# rounded up so as not to suggest less than the draws show. Draw by draw
# rather than at the draws' mean, since near-singular draws turned different
# ways average to a matrix of higher rank. NULL when nothing is bounded or
# there are no draws.
.ctCorBoundDrawRank <- function(fit, maxdraws = 400L) {
  blocks <- try(.ctCorBoundSpec(fit$model_spec), silent = TRUE)
  if (!inherits(blocks, "try-error"))
    blocks <- Filter(function(b) is.null(b$matrix), blocks)
  draws <- fit$estimate$rawposterior
  if (inherits(blocks, "try-error") || !length(blocks) || is.null(draws) ||
      !NROW(draws)) return(NULL)
  draws <- as.matrix(draws)
  use <- unique(round(seq(1, nrow(draws), length.out = min(maxdraws, nrow(draws)))))
  rows <- lapply(blocks, function(b) {
    ranks <- vapply(use, function(i) {
      r <- .ctImproveRankFromCov(.ctCorBoundMatrix(b, draws[i, ]))
      if (length(r) && is.finite(r)) as.numeric(r) else NA_real_
    }, numeric(1))
    if (all(is.na(ranks))) return(NULL)
    data.frame(level = b$level, k = b$k,
      rank = as.integer(ceiling(stats::median(ranks, na.rm = TRUE))),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# What a sample whose transitions stopped at the bound says about its
# covariance's rank (`$corbound_rank`, `.ctCorBoundDrawRank()`): for each level
# whose median draw points below its number of effects, that rank and the
# `poprank` that would fit it, as the fit's warning suggests at an estimate.
# Empty where there is none.
.ctSampleBoundRankNote <- function(diagnostics) {
  r <- diagnostics$corbound_rank
  if (is.null(r) || !NROW(r)) return("")
  r <- r[is.finite(r$rank) & r$rank < r$k, , drop = FALSE]
  if (!NROW(r)) return("")
  levels <- length(unique(diagnostics$corbound_rank$level)) > 1L
  paste0(" The population covariance is about rank ", r$rank, " of ", r$k,
    if (levels) paste0(" at level ", r$level) else "",
    " in the median draw: consider poprank = ", r$rank,
    " or fewer varying parameters.", collapse = "")
}

# One warning per level with an effect on the bound.
.ctCorBoundWarn <- function(report) {
  if (is.null(report) || !any(report$at_bound)) return(invisible(NULL))
  if (!is.null(report$matrix)) {
    own <- report$matrix != "population"
    for (m in unique(report$matrix[own & report$at_bound])) {
      r <- report[report$matrix == m, , drop = FALSE]
      warning(m, " correlations at their bound (multiple correlation ",
        attr(report, "bound"), ") for ", paste(r$effect[r$at_bound], collapse = ", "),
        ": the matrix is near singular, about rank ", r$rank[1L], " of ", r$k[1L],
        ". See fit$corbound.", call. = FALSE)
    }
    report <- report[!own, , drop = FALSE]
    if (!any(report$at_bound)) return(invisible(NULL))
  }
  levels <- unique(report$level[report$at_bound])
  for (lv in levels) {
    r <- report[report$level == lv, , drop = FALSE]
    warning("Random-effect correlations at their bound (multiple correlation ",
      attr(report, "bound"), ") for ", paste(r$effect[r$at_bound], collapse = ", "),
      if (length(unique(report$level)) > 1L) paste0(" at level ", lv),
      ": the population covariance is near singular, about rank ", r$rank[1L],
      " of ", r$k[1L], ". Consider poprank = ", r$rank[1L],
      " or fewer varying parameters. See fit$corbound.", call. = FALSE)
  }
  invisible(NULL)
}
