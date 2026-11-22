# Writing a covariance matrix into a model as the covariance it is meant to be.
#
# ctsem's covariance matrices -- T0VAR, DIFFUSION, MANIFESTVAR and the
# population spreads RAWPOPVAR / RAWPOPVAR_<level> -- are specified as a lower
# triangle: standard deviations on the diagonal and, below it, coordinates
# that the model's `covmattransform` turns into correlations. Only for two
# variables under 'z' is a coordinate a familiar quantity (Fisher's z); in
# general what one cell yields depends on the others, so a fixed matrix could
# not be written from a covariance at all. `ctCov()` inverts the construction.
#
# The cells depend on the construction, and `covmattransform` can change after
# the matrix is written, so a model remembers the covariance behind each
# `ctCov()` matrix (`model$covinput`) and `.ctCovRefresh()` rewrites the cells
# for whichever construction reads them: the model's own when fitting or
# generating on julia, a Cholesky factor for the R generator. A cell edited by
# hand since is the user's statement and is left alone.

#
# Under 'z' a cell of T0VAR, DIFFUSION or MANIFESTVAR that is fixed at zero is
# a zero covariance. Between variables that the matrix's other cells connect, a
# coordinate held at 0 would not be one (three variables, the other two cells
# at 0.5: a correlation of 0.115), so the construction solves those cells: the
# engine's cov_zero.jl says why there is always exactly one solution, and
# `.ctCovZeroSolve()` is the same computation for R's own readers. The
# population matrices keep every cell as the coordinate it is, since no syntax
# fixes one of their covariances at zero.

# Of the cells below the diagonal marked in the logical matrix `zero`, those
# whose two variables the other cells below the diagonal connect. Those are
# the ones the construction solves; a zero between separate blocks is a zero
# covariance as it stands.
#' @keywords internal
.ctCovLinkedZeros <- function(zero) {
  d <- nrow(zero)
  zero <- zero & lower.tri(zero)
  if (!any(zero)) return(zero)
  reach <- diag(d) > 0
  reach[lower.tri(zero) & !zero] <- TRUE
  reach <- reach | t(reach)
  repeat {
    more <- reach | (reach %*% reach > 0)
    if (identical(more, reach)) break
    reach <- more
  }
  zero & reach
}

# The hollow symmetric matrix `A` with the cells `zero` (below the diagonal)
# solved so that exp(A) is zero there, or NULL where there is no solution to
# find in double precision. Newton's method on log tr exp(A), which is convex
# in those cells, as the engine's `_zcov_solve!`.
#' @keywords internal
.ctCovZeroSolve <- function(A, zero) {
  d <- nrow(A)
  Z <- which(zero & lower.tri(zero), arr.ind = TRUE)
  z <- nrow(Z)
  if (!z) return(A)
  basis <- lapply(seq_len(z), function(k) {
    E <- matrix(0, d, d)
    E[Z[k, 1L], Z[k, 2L]] <- E[Z[k, 2L], Z[k, 1L]] <- 1
    E
  })
  evaluate <- function(A) {
    e <- eigen(A, symmetric = TRUE)
    w <- exp(e$values - e$values[1L])
    Y <- e$vectors %*% (w * t(e$vectors))
    list(e = e, w = w, g = Y[Z], total = sum(w),
      phi = e$values[1L] + log(sum(w)),
      r = max(abs(Y[Z]) / sqrt(diag(Y)[Z[, 1L]] * diag(Y)[Z[, 2L]])))
  }
  at <- evaluate(A)
  if (!is.finite(at$phi) || !is.finite(at$r)) return(NULL)
  for (iteration in seq_len(60L)) {
    if (at$r <= 1e-13) return(A)
    l <- at$e$values - at$e$values[1L]
    delta <- outer(l, l, "-")
    psi <- ifelse(abs(delta) < 1e-3, exp(outer(l, l, "+") / 2) * (1 + delta^2 / 24),
      outer(exp(l), exp(l), "-") / delta)
    V <- at$e$vectors
    H <- vapply(basis, function(E) (V %*% ((crossprod(V, E %*% V)) * psi) %*% t(V))[Z],
      numeric(z))
    step <- tryCatch(-solve(matrix(H, z, z), at$g), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) return(if (at$r <= 1e-8) A else NULL)
    slope <- 2 * sum(at$g * step) / at$total
    size <- 1
    accepted <- FALSE
    while (size > 1e-10) {
      trial <- A
      trial[Z] <- trial[Z] + size * step
      trial[Z[, 2:1, drop = FALSE]] <- trial[Z]
      new <- evaluate(trial)
      if (is.finite(new$phi) && is.finite(new$r) && new$phi <= at$phi +
          1e-4 * size * slope + 8 * .Machine$double.eps * (1 + abs(at$phi))) {
        accepted <- TRUE
        break
      }
      size <- size / 2
    }
    if (!accepted) return(if (at$r <= 1e-8) A else NULL)
    stalled <- at$r <= 1e-8 && new$r > at$r / 2
    A <- trial
    at <- new
    if (stalled) return(A)
  }
  if (at$r <= 1e-8) A else NULL
}

# Each construction the cells are read through, by name: the forward map from
# cells to covariance, as `sdcovsqrt2cov` (engine and stan) computes it.
# `zero` marks, below the diagonal, the cells that are zero covariances under
# 'z' (above), and is not read under the other constructions. NULL is none:
# every cell the coordinate it holds. Where the solve has no solution the
# result is NaN, as the engine's is.
#' @keywords internal
.ctCovForward <- function(cells, covmattransform, zero = NULL) {
  d <- nrow(cells)
  if (identical(covmattransform, "cholesky")) {
    low <- cells
    low[upper.tri(low)] <- 0
    return(tcrossprod(low))
  }
  sd <- diag(as.matrix(cells))
  if (identical(covmattransform, "z")) {
    A <- matrix(0, d, d)
    A[lower.tri(A)] <- cells[lower.tri(cells)]
    linked <- if (is.null(zero)) NULL else .ctCovLinkedZeros(zero)
    A <- A + t(A)
    if (any(linked)) {
      A[linked | t(linked)] <- 0
      A <- .ctCovZeroSolve(A, linked)
      if (is.null(A)) return(matrix(NaN, d, d))
    }
    Y <- expm::expm(A)
    out <- diag(sd, nrow = d) %*% stats::cov2cor(Y) %*% diag(sd, nrow = d)
    if (any(linked)) out[linked | t(linked)] <- 0
    return(out)
  }
  # 'rawcorr' and 'rawcorr_indep': the row-normalised correlation square root.
  O <- constraincorsqrt1(cells)
  tcrossprod(diag(sd, nrow = d) %*% O)
}

# The cells of a matrix of fixed numbers that are zero covariances under 'z':
# every one below the diagonal that is exactly zero.
#' @keywords internal
.ctCovZeroCells <- function(cells) {
  cells <- as.matrix(cells)
  lower.tri(cells) & !is.na(cells) & cells == 0
}

# The cells that give `cov` under `covmattransform`, checked by putting them
# back through the forward map.
#
# `zeros`: whether a cell of 0 is read as a zero covariance under 'z', as it
# is in T0VAR, DIFFUSION and MANIFESTVAR. Then a zero of `cov` is written as a
# cell of 0 and the other cells are found given that, and a covariance that is
# not zero is never written as exactly 0, which would be read as one. FALSE
# for the population matrices, whose cells are all coordinates.
#' @keywords internal
.ctCovCells <- function(cov, covmattransform = "z", what = "the covariance",
    zeros = TRUE) {
  cov <- as.matrix(cov)
  storage.mode(cov) <- "double"
  d <- nrow(cov)
  if (ncol(cov) != d || !all(is.finite(cov)) ||
      max(abs(cov - t(cov))) > 1e-10 * max(1, abs(cov))) {
    stop(what, " must be a finite symmetric matrix.", call. = FALSE)
  }
  cov <- (cov + t(cov)) / 2
  if (inherits(try(chol(cov), silent = TRUE), "try-error")) {
    stop(what, " is not positive definite.", call. = FALSE)
  }
  tf <- if (is.null(covmattransform)) "rawcorr" else as.character(covmattransform)
  cells <- matrix(0, d, d, dimnames = dimnames(cov))
  # The cells read as zero covariances, which stay at 0, and the cells to find.
  readzeros <- identical(tf, "z") && isTRUE(zeros)
  zero <- if (readzeros) lower.tri(cov) & cov == 0 else NULL
  open <- if (readzeros) lower.tri(cov) & !zero else lower.tri(cov)
  if (identical(tf, "cholesky")) {
    cells[] <- t(chol(cov))
  } else {
    if (identical(tf, "z")) {
      diag(cells) <- sqrt(diag(cov))
      start <- function(R) {
        e <- eigen(R, symmetric = TRUE)
        (e$vectors %*% diag(log(e$values), nrow = d) %*% t(e$vectors))[open]
      }
    } else {
      # Each row of the square root has norm 1 + 1e-5 rather than 1, so the
      # correlations the cells give are scaled by that, and the sds absorb it.
      diag(cells) <- sqrt(diag(cov) / (1 + 1e-5))
      start <- function(R) 2 * atanh(pmax(pmin(R[open], 0.99), -0.99))
    }
    if (any(open)) {
      target <- cov[open]
      residual <- function(x) {
        trial <- cells
        trial[open] <- x
        .ctCovForward(trial, tf, zero)[open] - target
      }
      x <- start(stats::cov2cor(cov))
      for (iteration in seq_len(60)) {
        r <- residual(x)
        if (max(abs(r)) < 1e-13 * max(1, abs(target))) break
        J <- vapply(seq_along(x), function(j) {
          h <- 1e-7 * max(1, abs(x[j]))
          e <- x
          e[j] <- e[j] + h
          (residual(e) - r) / h
        }, numeric(length(x)))
        step <- tryCatch(solve(J, r), error = function(e) NULL)
        if (is.null(step)) break
        # Halve a step that does not reduce the residual: the rawcorr map
        # saturates, and a full Newton step can overshoot into it.
        size <- 1
        while (size > 1e-6 && sum(residual(x - size * step)^2) >= sum(r^2)) {
          size <- size / 2
        }
        x <- x - size * step
      }
      # A cell of exactly 0 would be read as a zero covariance, which this
      # one is not: the smallest number is the same coordinate to any reader.
      if (readzeros) x[x == 0] <- .Machine$double.xmin
      cells[open] <- x
    }
  }
  # Through the forward map as a reader of these cells applies it.
  back <- .ctCovForward(cells, tf, if (readzeros) .ctCovZeroCells(cells))
  if (max(abs(back - cov)) > 1e-8 * max(1, abs(cov))) {
    stop(what, " cannot be written under covmattransform='", tf, "': the ",
      "nearest the construction reaches differs by ",
      format(signif(max(abs(back - cov)), 3)), ".",
      if (!identical(tf, "z")) " covmattransform='z' reaches every positive definite matrix.",
      call. = FALSE)
  }
  cells[upper.tri(cells)] <- 0
  cells
}

#' Write a covariance matrix as the cells a ctsem model reads
#'
#' ctsem's covariance matrices (\code{T0VAR}, \code{DIFFUSION},
#' \code{MANIFESTVAR}, and the population spreads \code{RAWPOPVAR} and
#' \code{RAWPOPVAR_<idname>}) are specified as a lower triangle, with standard
#' deviations on the diagonal and, below it, values that the model's
#' \code{covmattransform} turns into correlations. Except for two variables
#' under \code{'z'}, where each value is Fisher's z of the correlation, the
#' correlations depend on all the values together, so the cells for a given
#' covariance are hard to work out by hand. \code{ctCov()} computes them.
#'
#' Under \code{'z'}, a cell of \code{T0VAR}, \code{DIFFUSION} or
#' \code{MANIFESTVAR} that is fixed to 0 gives a covariance of zero, and
#' \code{ctCov()} writes a zero covariance as a cell of 0. Under
#' \code{'rawcorr'} and \code{'cholesky'} a cell of 0 gives a zero covariance
#' only in some positions. The cells of \code{RAWPOPVAR} are not treated this
#' way, and \code{backend = 'stan'} stops for a model with a cell fixed to 0
#' between variables that other covariances connect.
#'
#' An estimated cell that is 0, or near 0, does not show that the covariance
#' is zero when there are three or more variables. Read covariances and
#' correlations from the matrices that \code{\link{summary.ctStanFit}} and
#' \code{\link{ctSummaryMatrices}} report (\code{DIFFUSIONcov},
#' \code{DIFFUSIONcor}, \code{T0cov}, \code{T0cor}, \code{MANIFESTcov} and
#' \code{MANIFESTcor}), and not from the cells.
#'
#' Use the result wherever a fixed matrix can go, such as
#' \code{ctModel(DIFFUSION = ctCov(S))}, \code{m$matrices$T0VAR <- ctCov(S)},
#' or \code{m$RAWPOPVAR <- ctCov(S)}, always for the whole matrix. The model
#' also stores the covariance, so the matrix keeps its meaning if
#' \code{covmattransform} is changed later, and the R generator of
#' \code{\link{ctGenerate}}, which reads covariance cells as a Cholesky
#' factor, generates from the same covariance as the julia generator. A cell
#' edited by hand afterwards replaces what \code{ctCov()} wrote.
#'
#' For \code{RAWPOPVAR}, the covariance is that of the raw, untransformed
#' parameters across subjects (or groups, for \code{RAWPOPVAR_<idname>}),
#' which is the scale a fit reports it on.
#'
#' @param cov A symmetric positive definite covariance matrix.
#' @param covmattransform The construction to write the cells for:
#' \code{'z'} (the default, and what \code{\link{ctModel}} uses),
#' \code{'rawcorr'} or \code{'cholesky'}. A \code{ctCov()} matrix assigned into
#' a model is rewritten for that model's own construction, so this matters only
#' when the cells are read directly.
#' @return The cells, a lower triangular numeric matrix of class \code{ctCov},
#' with the covariance as attribute \code{covariance}. \code{'rawcorr'} cannot
#' reach every correlation matrix, and a covariance it cannot represent is
#' refused.
#' @examples
#' S <- matrix(c(1, .6, .3,
#'               .6, 2, .5,
#'               .3, .5, 1.5), 3)
#' ctCov(S)
#'
#' m <- ctModel(LAMBDA = diag(3), DIFFUSION = ctCov(S), Tpoints = 5)
#' @export
ctCov <- function(cov, covmattransform = "z") {
  covmattransform <- match.arg(covmattransform, c("z", "rawcorr", "cholesky"))
  cov <- as.matrix(cov)
  cells <- .ctCovCells(cov, covmattransform)
  structure(cells, class = c("ctCov", class(cells)), covariance = cov,
    covmattransform = covmattransform)
}

#' @export
print.ctCov <- function(x, ...) {
  cat("Cells for a covariance under covmattransform='",
    attr(x, "covmattransform"), "':\n", sep = "")
  print(unclass(structure(x, covariance = NULL, covmattransform = NULL)), ...)
  invisible(x)
}

# The fields a model's covariance cells can live in.
.ctCovSystemMatrices <- c("T0VAR", "DIFFUSION", "MANIFESTVAR")

# Whether a cell of 0 in the matrix `name` is read as a zero covariance under
# 'z': in the model's own matrices, not in the population's.
.ctCovReadsZeros <- function(name) name %in% .ctCovSystemMatrices

# The cells of a model's covariance matrix `name` that the construction solves
# as zero covariances: below the diagonal, fixed at zero, and joined by the
# matrix's other cells below the diagonal. A logical matrix, all FALSE under
# any construction but 'z'. From the model's `pars`, where a fixed cell has a
# value and no parameter.
#' @keywords internal
.ctCovSolvedCells <- function(model, name) {
  pars <- model$pars
  rows <- which(pars$matrix %in% name)
  d <- if (length(rows)) max(pars$row[rows], pars$col[rows]) else 0L
  none <- matrix(FALSE, d, d)
  tf <- if (is.null(model$covmattransform)) "rawcorr" else
    as.character(model$covmattransform)[1L]
  if (!identical(tf, "z") || d < 3L) return(none)
  zero <- none
  fixed <- rows[is.na(pars$param[rows]) & pars$value[rows] %in% 0 &
    pars$row[rows] > pars$col[rows]]
  zero[cbind(pars$row[fixed], pars$col[fixed])] <- TRUE
  .ctCovLinkedZeros(zero)
}

# The Stan program holds every cell as the coordinate it is, so a model with a
# cell the julia engine solves as a zero covariance is another model there:
# refused by name, with what to use.
#' @keywords internal
.ctCovRefuseSolved <- function(model) {
  for (name in .ctCovSystemMatrices) {
    cells <- which(.ctCovSolvedCells(model, name), arr.ind = TRUE)
    if (!nrow(cells)) next
    stop(name, "[", cells[1L, 1L], ",", cells[1L, 2L], "] is fixed at zero and ",
      "other covariances connect its two variables. backend = 'julia' fits ",
      "that as a zero covariance, which the stan backend cannot. Use backend = ",
      "'julia', or set model$covmattransform <- 'rawcorr', under which the ",
      "cell is a coordinate of the correlation matrix held at 0.", call. = FALSE)
  }
  invisible(NULL)
}

# The cells a model currently holds for `name`, as a numeric matrix, or NULL
# when any of them is not a fixed number.
#' @keywords internal
.ctCovCurrentCells <- function(model, name, d) {
  if (grepl("^RAWPOPVAR", name)) {
    held <- model[[name]]
    if (is.null(held) || nrow(held) != d) return(NULL)
    value <- suppressWarnings(matrix(as.numeric(held), d, d))
    return(if (anyNA(value[lower.tri(value, diag = TRUE)])) NULL else value)
  }
  pars <- model$pars
  rows <- which(pars$matrix %in% name & pars$row >= pars$col & pars$row <= d &
    pars$col <= d)
  if (!length(rows) || anyNA(pars$value[rows])) return(NULL)
  out <- matrix(0, d, d)
  out[cbind(pars$row[rows], pars$col[rows])] <- pars$value[rows]
  out
}

# Write `cells` into the model as fixed values of `name`.
#' @keywords internal
.ctCovWriteCells <- function(model, name, cells) {
  d <- nrow(cells)
  if (grepl("^RAWPOPVAR", name)) {
    held <- model[[name]]
    value <- matrix(sprintf("%.17g", cells), d, d)
    value[upper.tri(value)] <- "0"
    held[seq_len(d), seq_len(d)] <- value
    model[[name]] <- held
    return(model)
  }
  pars <- model$pars
  for (i in seq_len(d)) for (j in seq_len(i)) {
    row <- which(pars$matrix %in% name & pars$row == i & pars$col == j)
    if (!length(row)) next
    pars$value[row] <- cells[i, j]
    pars$param[row] <- NA
    pars$transform[row] <- NA
    if (!is.null(pars$indvarying)) pars$indvarying[row] <- FALSE
  }
  model$pars <- pars
  model
}

# A ctCov() matrix arriving at `name`: the cells for the model's construction,
# with the covariance remembered so `.ctCovRefresh()` can rewrite them.
#' @keywords internal
.ctCovAccept <- function(model, name, value) {
  cov <- attr(value, "covariance")
  tf <- if (is.null(model$covmattransform)) "rawcorr" else model$covmattransform
  cells <- .ctCovCells(cov, tf, what = name, zeros = .ctCovReadsZeros(name))
  dimnames(cells) <- dimnames(cov)
  inputs <- model$covinput
  inputs[[name]] <- list(cov = cov, cells = unname(cells), covmattransform = tf)
  model$covinput <- inputs
  list(model = model, cells = cells)
}

# Rewrite every remembered covariance for `covmattransform`, dropping any whose
# cells have been edited since.
#' @keywords internal
.ctCovRefresh <- function(model, covmattransform = model$covmattransform) {
  inputs <- model$covinput
  if (!length(inputs)) return(model)
  tf <- if (is.null(covmattransform)) "rawcorr" else as.character(covmattransform)
  for (name in names(inputs)) {
    entry <- inputs[[name]]
    d <- nrow(entry$cov)
    current <- .ctCovCurrentCells(model, name, d)
    if (is.null(current) || !isTRUE(all.equal(current, entry$cells,
        tolerance = 1e-10, check.attributes = FALSE))) {
      inputs[[name]] <- NULL
      next
    }
    if (identical(entry$covmattransform, tf)) next
    cells <- unname(.ctCovCells(entry$cov, tf, what = name,
      zeros = .ctCovReadsZeros(name)))
    model <- .ctCovWriteCells(model, name, cells)
    inputs[[name]] <- list(cov = entry$cov, cells = cells, covmattransform = tf)
  }
  model$covinput <- if (length(inputs)) inputs else NULL
  model
}
