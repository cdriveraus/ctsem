# The bound on random-effect correlations.
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
# effect's multiple correlation with the others at its level: at most 0.99,
# which for two effects is the correlation itself. The engine holds it with a
# wall in the objective (`_ctsem_corbound`, ctsem_backend.jl), on every julia
# route; review/PRIORS-default-boundary-2026-10-19.md has the evidence.
#
# This file says which blocks the wall applies to and at what bound, and
# reports a fit that ended on it.

# `optimcontrol$corbound`, resolved: the largest multiple correlation, 0.99 when
# not given, 0 (or FALSE) for no bound.
.ctCorBoundResolve <- function(x) {
  if (is.null(x)) return(0.99)
  if (isFALSE(x)) return(0)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x >= 1)
    stop("optimcontrol$corbound is the largest multiple correlation a random ",
      "effect may have with the others: a number in [0, 1), 0 for none.",
      call. = FALSE)
  as.numeric(x)
}

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
.ctCorBoundSpec <- function(spec,
    table = as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)) {
  bound <- .ctCorBoundResolve(spec$corbound)
  if (!(bound > 0)) return(NULL)
  code <- if (is.null(spec$covmatcode)) 0L else as.integer(spec$covmatcode)[1L]
  if (!code %in% c(0L, 2L)) return(NULL)
  lower <- function(k) {
    cells <- which(lower.tri(diag(k)), arr.ind = TRUE)
    cells[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
  }
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
        fixed = fixed, bound = bound, names = as.character(lv$param),
        level = if (is.null(lv$name)) "" else as.character(lv$name))
    }
    return(if (length(blocks)) blocks else NULL)
  }
  if (!any(table$matrix %in% "RAWPOPVAR")) return(NULL)
  if (!is.null(spec$model) && .ctPopFactorConstruction(spec$model)) return(NULL)
  pop <- table[table$matrix %in% "RAWPOPVAR", , drop = FALSE]
  row <- as.integer(pop$row); col <- as.integer(pop$col)
  k <- max(row)
  if (k < 2L) return(NULL)
  cells <- lower(k)
  index <- integer(nrow(cells)); fixed <- rep(NA_real_, nrow(cells))
  for (r in seq_len(nrow(cells))) {
    hit <- which(row == cells[r, 1L] & col == cells[r, 2L])
    if (length(hit) != 1L) return(NULL)
    p <- suppressWarnings(as.integer(pop$parnumber[hit]))
    if (!is.na(p) && p > 0L) {
      if (!identical(gsub(" ", "", as.character(pop$transform[hit])),
        paste0("param[", p, "]"))) return(NULL)
      index[r] <- p
    } else {
      fixed[r] <- as.numeric(pop$value[hit])
    }
  }
  diagonal <- which(row == col)
  names <- sub("^popsd_", "", as.character(pop$param[diagonal][order(row[diagonal])]))
  list(list(k = k, code = code, index = index, fixed = fixed, bound = bound,
    names = names,
    level = if (is.null(spec$subjectIDname)) "" else as.character(spec$subjectIDname)[1L]))
}

# Each bounded effect's multiple correlation with the others at the fit's
# estimate, whether it is on the bound, and the rank the block's correlation
# matrix points to (`.ctImproveRankFromCov`); NULL when nothing is bounded.
.ctCorBoundReport <- function(fit) {
  blocks <- try(.ctCorBoundSpec(fit$model_spec), silent = TRUE)
  raw <- as.numeric(fit$estimate$raw)
  if (inherits(blocks, "try-error") || !length(blocks) || !length(raw)) return(NULL)
  bound <- blocks[[1L]]$bound
  rows <- lapply(blocks, function(b) {
    cells <- diag(b$k)
    low <- which(lower.tri(cells), arr.ind = TRUE)
    low <- low[order(low[, 2L], low[, 1L]), , drop = FALSE]
    cells[low] <- ifelse(b$index > 0L, raw[pmax(b$index, 1L)], b$fixed)
    R <- .ctCovForward(cells, if (b$code == 2L) "z" else "rawcorr")
    inverse <- try(solve(R), silent = TRUE)
    mc <- if (inherits(inverse, "try-error")) rep(1, b$k) else
      sqrt(pmax(0, 1 - 1 / (diag(R) * diag(inverse))))
    # The wall lets an estimate sit a hair past the bound; within 5e-4 of it
    # is on it.
    data.frame(level = b$level, effect = b$names, multiple_correlation = mc,
      at_bound = mc >= bound - 5e-4, rank = .ctImproveRankFromCov(R), k = b$k,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "bound") <- bound
  out
}

# One warning per level with an effect on the bound.
.ctCorBoundWarn <- function(report) {
  if (is.null(report) || !any(report$at_bound)) return(invisible(NULL))
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
