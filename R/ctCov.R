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

# Each construction the cells are read through, by name: the forward map from
# cells to covariance, as `sdcovsqrt2cov` (engine and stan) computes it.
#' @keywords internal
.ctCovForward <- function(cells, covmattransform) {
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
    A <- A + t(A)
    Y <- expm::expm(A)
    return(diag(sd, nrow = d) %*% stats::cov2cor(Y) %*% diag(sd, nrow = d))
  }
  # 'rawcorr' and 'rawcorr_indep': the row-normalised correlation square root.
  O <- constraincorsqrt1(cells)
  tcrossprod(diag(sd, nrow = d) %*% O)
}

# The cells that give `cov` under `covmattransform`, checked by putting them
# back through the forward map.
#' @keywords internal
.ctCovCells <- function(cov, covmattransform = "z", what = "the covariance") {
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
  if (identical(tf, "cholesky")) {
    cells[] <- t(chol(cov))
  } else {
    if (identical(tf, "z")) {
      diag(cells) <- sqrt(diag(cov))
      start <- function(R) {
        e <- eigen(R, symmetric = TRUE)
        (e$vectors %*% diag(log(e$values), nrow = d) %*% t(e$vectors))[lower.tri(R)]
      }
    } else {
      # Each row of the square root has norm 1 + 1e-5 rather than 1, so the
      # correlations the cells give are scaled by that, and the sds absorb it.
      diag(cells) <- sqrt(diag(cov) / (1 + 1e-5))
      start <- function(R) 2 * atanh(pmax(pmin(R[lower.tri(R)], 0.99), -0.99))
    }
    if (d > 1L) {
      target <- cov[lower.tri(cov)]
      residual <- function(x) {
        trial <- cells
        trial[lower.tri(trial)] <- x
        .ctCovForward(trial, tf)[lower.tri(cov)] - target
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
      cells[lower.tri(cells)] <- x
    }
  }
  back <- .ctCovForward(cells, tf)
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
#' ctsem's covariance matrices -- \code{T0VAR}, \code{DIFFUSION},
#' \code{MANIFESTVAR}, and the population spreads \code{RAWPOPVAR} and
#' \code{RAWPOPVAR_<idname>} -- are specified as a lower triangle: standard
#' deviations on the diagonal, and below it coordinates that the model's
#' \code{covmattransform} turns into correlations. Except for two variables
#' under \code{'z'}, where a coordinate is Fisher's z, what one coordinate
#' yields depends on the others, so the cells for a given covariance are not
#' something to work out by hand. \code{ctCov()} computes them.
#'
#' Pass the result wherever a fixed matrix goes: \code{ctModel(DIFFUSION =
#' ctCov(S))}, \code{m$matrices$T0VAR <- ctCov(S)}, or \code{m$RAWPOPVAR <-
#' ctCov(S)}, in every case for the whole matrix. The model keeps the
#' covariance as well as the cells, so the matrix keeps its meaning if
#' \code{covmattransform} is changed afterwards, and the R generator of
#' \code{\link{ctGenerate}} -- which reads covariance cells as a Cholesky
#' factor -- generates from the same covariance as the julia one. A cell edited
#' by hand afterwards replaces what \code{ctCov()} wrote.
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
  cells <- .ctCovCells(cov, tf, what = name)
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
    cells <- unname(.ctCovCells(entry$cov, tf, what = name))
    model <- .ctCovWriteCells(model, name, cells)
    inputs[[name]] <- list(cov = entry$cov, cells = cells, covmattransform = tf)
  }
  model$covinput <- if (length(inputs)) inputs else NULL
  model
}
