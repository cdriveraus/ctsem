# Network representations of a continuous time model --------------------------
#
# Network psychometrics reads a multivariate time series as two graphs: a
# *temporal* one, the directed lag-1 regression matrix, and a *contemporaneous*
# one, the partial correlations of the innovations. graphicalVAR, mlVAR and
# qgraph all draw those two.
#
# A continuous time model has both, but neither is a fixed object: both are
# functions of the interval you ask about. The temporal network is
# expm(DRIFT * dt) -- one matrix per interval, and its *shape* changes with dt,
# not only its magnitude, because a path that runs a -> b -> c needs time to
# arrive and is invisible at short dt and dominant at long dt. The
# contemporaneous network is the partial correlation structure of the
# innovation covariance accumulated over dt,
#
#     Q(dt) = int_0^dt expm(DRIFT s) DIFFUSIONcov expm(DRIFT s)' ds,
#
# which tends to the partial correlations of DIFFUSIONcov itself as dt -> 0
# (partial correlations are scale free, so the dt factor drops out) and to
# those of the asymptotic covariance as dt -> Inf. So a discrete time VAR's
# two networks are one slice through a one-parameter family, and the slice
# depends on how often that study happened to measure.
#
# That is the whole point of drawing these from a ctsem model rather than from
# a VAR fit, and it is why `dt` is the first argument after the model.
#
# Q(dt) is computed by the Van Loan block-exponential rather than the
# stationary identity asymDIFFUSIONcov - dtDRIFT asymDIFFUSIONcov dtDRIFT',
# because the block form needs no stationarity: a model with an unstable or
# singular DRIFT still has a perfectly well defined innovation covariance over
# a finite interval, and refusing to draw its contemporaneous network would
# refuse exactly the models where the picture is most wanted.

# Covariance from ctsem's DIFFUSION / T0VAR / MANIFESTVAR parameterisation.
#
# The lower triangle is an unconstrained correlation square root and the
# diagonal holds standard deviations, so this mirrors the stan function
# `sdcovsqrt2cov()` in R/ctModelWriter.R -- including its epsilon, which is why
# an sd of 0.5 gives a variance of 0.2500025 rather than 0.25. Used only on the
# model-specification path; a fit already reports DIFFUSIONcov.
.ctSdcovsqrt2cov <- function(mat, cholesky = FALSE){
  mat <- as.matrix(mat)
  if(nrow(mat) == 0) return(mat)
  if(cholesky) return(tcrossprod(mat))
  tcrossprod(diag(diag(mat), nrow = nrow(mat)) %*% constraincorsqrt1(mat))
}

# Partial correlations of a covariance matrix -- the Gaussian graphical model
# every contemporaneous network in the literature is.
#
# A process with no noise of its own (a zero DIFFUSION row, as in the second
# order models where one latent is the derivative of another) makes the
# covariance singular. Those processes are dropped from the inversion and get
# zero edges, rather than the whole matrix failing.
.ctNetworkPcor <- function(cov){
  cov <- as.matrix(cov)
  n <- nrow(cov)
  out <- matrix(0, n, n, dimnames = dimnames(cov))
  if(n == 0) return(out)
  d <- diag(cov)
  keep <- which(is.finite(d) & d > max(c(d, 0), na.rm = TRUE) * 1e-10)
  if(length(keep) < 2) return(out)
  sub <- cov[keep, keep, drop = FALSE]
  prec <- try(solve(sub), silent = TRUE)
  if(inherits(prec, 'try-error')) return(out)
  p <- -prec / sqrt(outer(diag(prec), diag(prec)))
  diag(p) <- 0
  out[keep, keep] <- p
  out
}

# Asymptotic (stationary) within-subject covariance from DRIFT and
# DIFFUSIONcov, for the model-specification path only -- a fit reports
# asymDIFFUSIONcov directly. Returns NULL when the system is not stationary,
# because then there is nothing to return rather than a number to distrust.
.ctNetworkAsym <- function(drift, diffusioncov, continuoustime = TRUE){
  drift <- as.matrix(drift); q <- as.matrix(diffusioncov)
  n <- nrow(drift)
  if(n == 0) return(NULL)
  ev <- try(eigen(drift, only.values = TRUE)$values, silent = TRUE)
  if(inherits(ev, 'try-error')) return(NULL)
  stable <- if(continuoustime) max(Re(ev)) < 0 else max(Mod(ev)) < 1
  if(!stable) return(NULL)
  # vec(DRIFT X + X DRIFT') = (I kron DRIFT + DRIFT kron I) vec(X) in column
  # major order, so the continuous time Lyapunov equation is one linear solve.
  I <- diag(n)
  sol <- if(continuoustime){
    try(solve(kronecker(I, drift) + kronecker(drift, I), -as.vector(q)), silent = TRUE)
  } else try(solve(diag(n * n) - kronecker(drift, drift), as.vector(q)), silent = TRUE)
  if(inherits(sol, 'try-error')) return(NULL)
  x <- matrix(sol, n, n)
  x <- (x + t(x)) / 2
  if(any(!is.finite(diag(x))) || any(diag(x) < 0)) return(NULL)
  dimnames(x) <- dimnames(q)
  x
}

# Innovation (process error) covariance accumulated over an interval.
#
# Continuous time uses Van Loan's block exponential: with
# B = [[DRIFT, Q], [0, -DRIFT']] * dt, the top blocks of expm(B) give
# F11 = expm(DRIFT dt) and F12, and F12 F11' is the integral. No stationarity
# assumption, and correct for a singular DRIFT.
.ctNetworkInnovation <- function(drift, diffusioncov, dt, continuoustime = TRUE){
  drift <- as.matrix(drift); q <- as.matrix(diffusioncov)
  n <- nrow(drift)
  if(n == 0) return(q)
  if(continuoustime){
    B <- rbind(cbind(drift, q), cbind(matrix(0, n, n), -t(drift))) * dt
    E <- as.matrix(expm::expm(B))
    out <- E[1:n, n + (1:n), drop = FALSE] %*% t(E[1:n, 1:n, drop = FALSE])
  } else {
    steps <- as.integer(round(dt))
    if(steps < 1) stop(call. = FALSE,
      'A discrete time model needs dt to be a positive whole number of steps.')
    out <- matrix(0, n, n)
    ak <- diag(n)
    for(k in seq_len(steps)){
      out <- out + ak %*% q %*% t(ak)
      ak <- ak %*% drift
    }
  }
  out <- (out + t(out)) / 2
  dimnames(out) <- dimnames(q)
  out
}

# DRIFT, DIFFUSIONcov, asymDIFFUSIONcov and LAMBDA as point values, from either
# a model specification or a fit on either backend.
#
# The fit path is `ctSummaryMatrices()`, which is one implementation for stan
# and julia and reports the posterior median. The model path has no posterior,
# so a fixed cell is its own value and a free cell is the value the parameter
# takes at a raw value of zero -- which is where ctsem's own optimiser starts,
# and is the only defensible reading of "the model as specified" for something
# nobody has estimated yet.
.ctNetworkInputs <- function(x, quiet = FALSE, ...){
  if(inherits(x, 'ctStanModel')){
    m <- ctModelTransformsToNum(x)
    pars <- m$pars
    # A cell whose label is an expression rather than a plain parameter name is
    # a specification, not a number waiting for one -- a state-dependent DRIFT
    # cell is the usual case. Filling those would draw a graph of a model that
    # was not written, so they are reported and left out.
    reserved <- c(m$latentNames, m$manifestNames, m$TDpredNames, m$TIpredNames)
    label <- !is.na(pars$param) &
      grepl('^[A-Za-z.][A-Za-z0-9._]*$', as.character(pars$param)) &
      !as.character(pars$param) %in% reserved
    # Only the three matrices the networks are built from. Reporting a free
    # T0MEANS or MANIFESTVAR as filled in would be true and useless: it cannot
    # move an edge, and it is most of the parameters in a typical model.
    free <- which(is.na(pars$value) & label &
        as.character(pars$matrix) %in% c('DRIFT', 'DIFFUSION', 'LAMBDA'))
    for(i in free) pars$value[i] <- tform(0, pars$transform[i], pars$multiplier[i],
      pars$meanscale[i], pars$offset[i], pars$inneroffset[i])
    m$pars <- pars
    mats <- listOfMatrices(m$pars)
    tonum <- function(name, nrow, ncol){
      raw <- mats[[name]]
      if(is.null(raw)) return(matrix(0, nrow, ncol))
      v <- suppressWarnings(as.numeric(raw))
      bad <- sum(is.na(v))
      if(bad > 0 && !quiet) message(bad, ' ', name,
        ' cell(s) depend on the latent state or another parameter and are shown as absent; ',
        'fit the model, or see ctPhasePortrait() for what a state dependent system does.')
      v[is.na(v)] <- 0
      matrix(v, nrow(raw), ncol(raw))[seq_len(nrow), seq_len(ncol), drop = FALSE]
    }
    nl <- m$n.latent; nm <- m$n.manifest
    drift <- tonum('DRIFT', nl, nl)
    diffusion <- tonum('DIFFUSION', nl, nl)
    lambda <- tonum('LAMBDA', nm, nl)
    diffusioncov <- .ctSdcovsqrt2cov(diffusion,
      cholesky = identical(as.character(m$covmattransform)[1], 'cholesky'))
    if(length(free) && !quiet) message(length(free), ' free ',
      paste(sort(unique(as.character(pars$matrix[free]))), collapse = '/'),
      ' parameter(s) were taken at their untransformed value of zero, which ',
      'is where estimation starts. Fix them in the model, or fit it, for a ',
      'graph of anything else.')
    out <- list(DRIFT = drift, DIFFUSIONcov = diffusioncov,
      asymDIFFUSIONcov = .ctNetworkAsym(drift, diffusioncov, m$continuoustime),
      LAMBDA = lambda, latentNames = m$latentNames[seq_len(nl)],
      manifestNames = m$manifestNames[seq_len(nm)],
      continuoustime = isTRUE(m$continuoustime), source = 'model')
  } else if(inherits(x, 'ctStanFit') || inherits(x, 'ctJuliaFit') || inherits(x, 'ctFit')){
    m <- .ctFitModelObject(x)
    mats <- suppressMessages(ctSummaryMatrices(x, ...))
    nl <- .ctFitNlatent(x)
    nm <- length(m$manifestNames)
    trim <- function(mat, nrow, ncol) as.matrix(mat)[seq_len(nrow), seq_len(ncol), drop = FALSE]
    asym <- if(is.null(mats$asymDIFFUSIONcov)) NULL else trim(mats$asymDIFFUSIONcov, nl, nl)
    if(!is.null(asym) && (any(!is.finite(diag(asym))) || any(diag(asym) < 0))) asym <- NULL
    out <- list(DRIFT = trim(mats$DRIFT, nl, nl),
      DIFFUSIONcov = trim(mats$DIFFUSIONcov, nl, nl),
      asymDIFFUSIONcov = asym, LAMBDA = trim(mats$LAMBDA, nm, nl),
      latentNames = m$latentNames[seq_len(nl)],
      manifestNames = m$manifestNames[seq_len(nm)],
      continuoustime = isTRUE(m$continuoustime), source = 'fit')
  } else stop(call. = FALSE, paste0('ctNetwork() needs a ctModel(type="ct"/"dt") ',
    'specification or a fit from ctFit(); got ', paste(class(x), collapse = '/'), '.'))

  ln <- out$latentNames
  dimnames(out$DRIFT) <- dimnames(out$DIFFUSIONcov) <- list(ln, ln)
  if(!is.null(out$asymDIFFUSIONcov)) dimnames(out$asymDIFFUSIONcov) <- list(ln, ln)
  dimnames(out$LAMBDA) <- list(out$manifestNames, ln)
  out
}

# Matrix -> tidy edges. `mat[i,j]` is always the effect of column j on row i, so
# an edge runs from the column name to the row name.
.ctNetworkEdges <- function(mat, network, directed, threshold, selfloops = TRUE){
  if(is.null(mat) || !length(mat)) return(NULL)
  idx <- which(abs(mat) > threshold & is.finite(mat), arr.ind = TRUE)
  if(!nrow(idx)) return(NULL)
  keep <- rep(TRUE, nrow(idx))
  if(!selfloops) keep <- keep & rownames(mat)[idx[, 1]] != colnames(mat)[idx[, 2]]
  # An undirected network is symmetric; one edge per pair, not two.
  if(!directed) keep <- keep & idx[, 1] >= idx[, 2]
  idx <- idx[keep, , drop = FALSE]
  if(!nrow(idx)) return(NULL)
  data.frame(network = network, from = colnames(mat)[idx[, 2]],
    to = rownames(mat)[idx[, 1]], weight = mat[idx], directed = directed,
    self = rownames(mat)[idx[, 1]] == colnames(mat)[idx[, 2]],
    stringsAsFactors = FALSE)
}

#' Temporal and contemporaneous networks of a continuous time model
#'
#' @description
#' The two graphs network psychometrics draws from a multivariate time series --
#' the directed temporal network and the undirected contemporaneous network --
#' computed from a ctsem model specification or from a fit on either backend.
#'
#' Both are functions of the time interval, which is what a continuous time
#' model has to say here that a discrete time VAR does not. Call it at several
#' \code{dt} and the networks change shape, not merely scale.
#'
#' @details
#' Three networks over the latent processes, plus the measurement model:
#'
#' \describe{
#'   \item{\code{temporal}}{\code{expm(DRIFT * dt)}, the discrete time
#'     autoregression and cross-lagged matrix at interval \code{dt}. Directed,
#'     with the diagonal as self-loops (the autoregressions). Cell
#'     \code{[i,j]} is the effect of process \code{j} now on process \code{i}
#'     one interval later, so an edge runs from the column to the row. This is
#'     \code{\link{ctDiscretePars}} at a single interval, and shares its
#'     \code{observational} and \code{standardise} arguments.}
#'   \item{\code{contemporaneous}}{Partial correlations of the innovation
#'     covariance accumulated over \code{dt},
#'     \code{Q(dt) = int_0^dt expm(DRIFT s) DIFFUSIONcov t(expm(DRIFT s)) ds}.
#'     This is the continuous time counterpart of the graphicalVAR
#'     contemporaneous network. Partial correlations are scale free, so as
#'     \code{dt} shrinks this tends to the partial correlations of
#'     \code{DIFFUSIONcov} itself, and as \code{dt} grows to those of the
#'     asymptotic covariance.}
#'   \item{\code{asymptotic}}{Partial correlations of
#'     \code{asymDIFFUSIONcov}, the long run within-subject covariance -- the
#'     network a cross-sectional study of a stationary system would see.
#'     \code{NULL} if the system is not stationary.}
#'   \item{\code{measurement}}{\code{LAMBDA}, as a bipartite manifest by latent
#'     graph. Directed from latent to manifest.}
#' }
#'
#' From a fit, the matrices are the posterior medians
#' \code{\link{ctSummaryMatrices}} reports. From an unfitted specification, a
#' fixed cell is its own value and a free cell takes the value its transform
#' gives at a raw parameter of zero, which is where estimation starts; how many
#' DRIFT, DIFFUSION and LAMBDA cells were treated that way is reported, those
#' being the only matrices an edge can come from.
#'
#' @param x A \code{ctStanModel} from \code{\link{ctModel}}, or a fit from
#'   \code{\link{ctFit}} on either backend.
#' @param dt Time interval for the temporal and contemporaneous networks. A
#'   whole number of steps for a discrete time model. Length one; to see the
#'   interval dependence, call \code{\link{ctNetworkPlot}} with a vector.
#' @param standardise If TRUE (the default), temporal edges are in standard
#'   deviation units of the processes themselves, using
#'   \code{asymDIFFUSIONcov}, so edges between processes on different scales are
#'   comparable. Falls back to FALSE, with a message, when the system has no
#'   stationary variance to standardise by.
#' @param observational What a one unit change in a process brings with it; see
#'   \code{\link{ctDiscretePars}}, whose argument this is. \code{FALSE} (the
#'   default) is the partial regression, a property of the dynamics alone.
#' @param threshold Edges with \code{abs(weight)} at or below this are dropped
#'   from \code{$edges} and from any plot. The matrices are never thresholded.
#' @param networks Which networks appear in \code{$edges} and in the plot. Any
#'   of \code{'temporal'}, \code{'contemporaneous'}, \code{'asymptotic'},
#'   \code{'measurement'}. All four matrices are computed regardless.
#' @param plot If TRUE, returns the \code{\link{ctNetworkPlot}} of the result
#'   instead of the result.
#' @param quiet Suppress the messages about how the numbers were arrived at.
#' @param ... Passed to \code{\link{ctSummaryMatrices}} when \code{x} is a fit,
#'   and to \code{\link{ctNetworkPlot}} when \code{plot=TRUE}.
#'
#' @return A list of class \code{ctNetwork}:
#' \describe{
#'   \item{\code{temporal}, \code{contemporaneous}, \code{asymptotic},
#'     \code{measurement}}{Weight matrices with the model's own names, ready to
#'     hand to \code{qgraph::qgraph()}. \code{asymptotic} may be \code{NULL}.}
#'   \item{\code{edges}}{Tidy data frame: \code{network}, \code{from},
#'     \code{to}, \code{weight}, \code{directed}, \code{self}. Thresholded, and
#'     restricted to \code{networks}.}
#'   \item{\code{nodes}}{Data frame of \code{name} and \code{type}
#'     (\code{'latent'} or \code{'manifest'}).}
#'   \item{\code{DRIFT}, \code{DIFFUSIONcov}, \code{asymDIFFUSIONcov},
#'     \code{innovation}}{The matrices the networks were derived from, so every
#'     number above can be checked.}
#' }
#'
#' @seealso \code{\link{ctNetworkPlot}} to draw it,
#'   \code{\link{ctDiscretePars}} for the temporal network across a continuum of
#'   intervals, \code{\link{ctPhasePortrait}} for a nonlinear system's vector
#'   field.
#'
#' @examples
#' # Two processes, b driven by a, drawn at one interval.
#' m <- ctModel(type='ct', n.latent=2, n.manifest=2,
#'   manifestNames=c('y1','y2'), latentNames=c('a','b'),
#'   LAMBDA=diag(2),
#'   DRIFT=matrix(c(-1, 0, .5, -2), 2, 2, byrow=TRUE),
#'   DIFFUSION=matrix(c(1, 0, .3, 1), 2, 2, byrow=TRUE))
#'
#' net <- ctNetwork(m, dt=1, standardise=FALSE)
#' round(net$temporal, 4)          # exp(-1) and exp(-2) on the diagonal
#' round(net$contemporaneous, 4)
#' net$edges
#'
#' @export
ctNetwork <- function(x, dt = 1, standardise = TRUE, observational = FALSE,
  threshold = 0, networks = c('temporal', 'contemporaneous'), plot = FALSE,
  quiet = FALSE, ...){

  if(length(dt) != 1 || !is.finite(dt) || dt <= 0) stop(call. = FALSE,
    'dt must be a single positive number. ctNetworkPlot() takes a vector of them.')
  networks <- match.arg(networks,
    c('temporal', 'contemporaneous', 'asymptotic', 'measurement'), several.ok = TRUE)

  dots <- list(...)
  plotargs <- dots[names(dots) %in% names(formals(ctNetworkPlot))]
  summaryargs <- dots[!names(dots) %in% names(formals(ctNetworkPlot))]

  inputs <- do.call(.ctNetworkInputs, c(list(x, quiet = quiet), summaryargs))
  if(!inputs$continuoustime && abs(dt - round(dt)) > 1e-8) stop(call. = FALSE,
    'This is a discrete time model, so dt must be a whole number of steps.')

  if(standardise && is.null(inputs$asymDIFFUSIONcov)){
    if(!quiet) message('standardise=TRUE needs a stationary variance to standardise by, ',
      'and this system has none; reporting unstandardised temporal edges.')
    standardise <- FALSE
  }

  # The temporal network goes through ctDiscreteParsDrift() rather than a local
  # expm() call, so that `observational` and `standardise` mean exactly what
  # they mean everywhere else in the package, with one implementation.
  as4d <- function(mat) array(mat, dim = c(1, 1, dim(mat)))
  ctpars <- list(DRIFT = as4d(inputs$DRIFT), DIFFUSIONcov = as4d(inputs$DIFFUSIONcov),
    asymDIFFUSIONcov = as4d(if(is.null(inputs$asymDIFFUSIONcov))
      inputs$DIFFUSIONcov else inputs$asymDIFFUSIONcov))
  temporal <- ctDiscreteParsDrift(ctpars, times = dt, observational = observational,
    standardise = standardise, cov = FALSE, discreteInput = !inputs$continuoustime,
    quiet = TRUE)
  temporal <- matrix(temporal[1, 1, 1, , ], nrow(inputs$DRIFT),
    dimnames = dimnames(inputs$DRIFT))

  innovation <- .ctNetworkInnovation(inputs$DRIFT, inputs$DIFFUSIONcov, dt,
    inputs$continuoustime)
  contemporaneous <- .ctNetworkPcor(innovation)
  asymptotic <- if(is.null(inputs$asymDIFFUSIONcov)) NULL else
    .ctNetworkPcor(inputs$asymDIFFUSIONcov)

  out <- list(temporal = temporal, contemporaneous = contemporaneous,
    asymptotic = asymptotic, measurement = inputs$LAMBDA,
    DRIFT = inputs$DRIFT, DIFFUSIONcov = inputs$DIFFUSIONcov,
    asymDIFFUSIONcov = inputs$asymDIFFUSIONcov, innovation = innovation)

  edges <- do.call(rbind, c(
    if('temporal' %in% networks) list(.ctNetworkEdges(temporal, 'temporal',
      directed = TRUE, threshold = threshold)),
    if('contemporaneous' %in% networks) list(.ctNetworkEdges(contemporaneous,
      'contemporaneous', directed = FALSE, threshold = threshold, selfloops = FALSE)),
    if('asymptotic' %in% networks) list(.ctNetworkEdges(asymptotic, 'asymptotic',
      directed = FALSE, threshold = threshold, selfloops = FALSE)),
    if('measurement' %in% networks) list(.ctNetworkEdges(inputs$LAMBDA,
      'measurement', directed = TRUE, threshold = threshold))))
  if(is.null(edges)) edges <- data.frame(network = character(), from = character(),
    to = character(), weight = numeric(), directed = logical(), self = logical(),
    stringsAsFactors = FALSE)
  rownames(edges) <- NULL
  out$edges <- edges

  out$nodes <- rbind(
    data.frame(name = inputs$latentNames, type = 'latent', stringsAsFactors = FALSE),
    data.frame(name = inputs$manifestNames, type = 'manifest', stringsAsFactors = FALSE))

  attributes(out)$dt <- dt
  attributes(out)$standardise <- standardise
  attributes(out)$observational <- observational
  attributes(out)$threshold <- threshold
  attributes(out)$networks <- networks
  attributes(out)$continuoustime <- inputs$continuoustime
  attributes(out)$source <- inputs$source
  class(out) <- c('ctNetwork', 'list')

  if(plot) return(do.call(ctNetworkPlot, c(list(out), plotargs)))
  out
}

#' @export
print.ctNetwork <- function(x, ...){
  cat('ctNetwork from a ', attributes(x)$source, ', dt = ',
    format(attributes(x)$dt), ', ',
    if(attributes(x)$standardise) 'standardised' else 'unstandardised', '\n',
    sep = '')
  cat(sum(x$nodes$type == 'latent'), ' latent processes: ',
    paste(x$nodes$name[x$nodes$type == 'latent'], collapse = ', '), '\n', sep = '')
  if(nrow(x$edges)){
    tab <- table(x$edges$network)
    cat('edges past threshold ', format(attributes(x)$threshold), ': ',
      paste0(names(tab), ' ', as.integer(tab), collapse = ', '), '\n', sep = '')
  } else cat('no edges past threshold ', format(attributes(x)$threshold), '\n', sep = '')
  cat('matrices: temporal, contemporaneous',
    if(!is.null(x$asymptotic)) ', asymptotic' else '',
    ', measurement; edge list in $edges\n', sep = '')
  invisible(x)
}


# Drawing ---------------------------------------------------------------------
#
# ggplot2 is already a ctsem dependency and qgraph is not, so the default
# renderer is written here rather than adding a dependency for a picture. The
# matrices are returned in the shape qgraph wants, so `engine='qgraph'` is a
# thin passthrough for anyone who has it, and
# `qgraph::qgraph(net$temporal, directed=TRUE)` works without this function at
# all.

# Node coordinates. Latents on a unit circle; each manifest, when the
# measurement network is drawn, just outside the latent it loads most strongly
# on, fanned out when several share one -- so the measurement panel reads as the
# factor structure rather than as two unrelated rings.
.ctNetworkLayout <- function(nodes, lambda = NULL){
  polar <- function(names, radius, ang) data.frame(name = names,
    x = radius * cos(ang), y = radius * sin(ang), angle = ang,
    stringsAsFactors = FALSE)
  latents <- nodes$name[nodes$type == 'latent']
  n <- length(latents)
  # Two nodes go left and right rather than top and bottom, which uses a panel's
  # width instead of fighting it; three or more start at the top.
  start <- if(n == 2) pi else pi / 2
  ang <- if(n == 1) pi / 2 else start - 2 * pi * (seq_len(n) - 1) / n
  out <- polar(latents, 1, ang)
  manifests <- nodes$name[nodes$type == 'manifest']
  if(!length(manifests)) return(out)
  parent <- rep(NA_integer_, length(manifests))
  if(!is.null(lambda) && nrow(lambda)){
    hit <- match(manifests, rownames(lambda))
    for(j in seq_along(manifests)){
      if(is.na(hit[j])) next
      load <- abs(lambda[hit[j], ])
      if(any(load > 0)) parent[j] <- which.max(load)
    }
  }
  mang <- numeric(length(manifests))
  loose <- which(is.na(parent))
  for(k in unique(stats::na.omit(parent))){
    group <- which(parent == k)
    spread <- if(length(group) > 1) seq(-0.4, 0.4, length.out = length(group)) else 0
    mang[group] <- ang[k] + spread
  }
  if(length(loose)) mang[loose] <- pi / 2 - 2 * pi * (seq_along(loose) - 1) / length(loose)
  rbind(out, polar(manifests, 1.7, mang))
}

# Edge endpoints, pulled back from the node centres so an arrowhead lands on the
# node boundary rather than under the label. Self-loops become a short segment
# offset around the node, which geom_curve's curvature closes into a loop.
.ctNetworkEdgeGeom <- function(edges, layout, gap = 0.16, loop = 0.2){
  if(!nrow(edges)) return(edges)
  pos <- layout[match(edges$from, layout$name), c('x', 'y', 'angle')]
  end <- layout[match(edges$to, layout$name), c('x', 'y')]
  out <- edges
  dx <- end$x - pos$x; dy <- end$y - pos$y
  len <- sqrt(dx^2 + dy^2); len[len == 0] <- 1
  out$x <- pos$x + dx / len * gap
  out$y <- pos$y + dy / len * gap
  out$xend <- end$x - dx / len * gap
  out$yend <- end$y - dy / len * gap
  if(any(edges$self)){
    s <- which(edges$self)
    a <- pos$angle[s]
    r <- sqrt(pos$x[s]^2 + pos$y[s]^2) + loop
    spread <- 0.115
    out$x[s] <- r * cos(a - spread); out$y[s] <- r * sin(a - spread)
    out$xend[s] <- r * cos(a + spread); out$yend[s] <- r * sin(a + spread)
  }
  out
}

.ctNetworkGG <- function(edges, nodes, lambda, maxwidth, title, arrowsize, nodesize,
  labelsize, poscolour, negcolour, ncol){
  layout <- .ctNetworkLayout(nodes, lambda)
  geom <- .ctNetworkEdgeGeom(edges, layout)
  panels <- levels(droplevels(edges$panel))
  layout$type <- nodes$type[match(layout$name, nodes$name)]
  # A latent process belongs in every panel even with no edges past the
  # threshold -- an isolated node is a finding. A manifest belongs only in a
  # panel that draws the measurement model, or it is decoration.
  nodedf <- do.call(rbind, lapply(panels, function(p){
    used <- unlist(edges[edges$panel %in% p, c('from', 'to')])
    keep <- layout$type %in% 'latent' | layout$name %in% used
    if(!any(keep)) return(NULL)
    cbind(layout[keep, , drop = FALSE], panel = p)
  }))
  nodedf$panel <- factor(nodedf$panel, levels = panels)

  # Three edge layers, because geom_curve takes curvature as a parameter rather
  # than an aesthetic: straight for undirected, gently curved for directed
  # (which separates a -> b from b -> a, since curvature is relative to each
  # segment's own direction), and tightly curved for self-loops.
  head <- ggplot2::arrow(length = ggplot2::unit(arrowsize, 'mm'), type = 'closed')
  layer <- function(rows, curvature, arrow){
    d <- geom[rows, , drop = FALSE]
    if(!nrow(d)) return(NULL)
    ggplot2::geom_curve(data = d, curvature = curvature, arrow = arrow,
      lineend = 'round', alpha = .85, ggplot2::aes(x = .data$x, y = .data$y,
        xend = .data$xend, yend = .data$yend, linewidth = abs(.data$weight),
        colour = .data$weight > 0))
  }
  ggplot2::ggplot() +
    layer(!geom$directed, 0, NULL) +
    layer(geom$directed & !geom$self, 0.13, head) +
    layer(geom$self, -3, head) +
    ggplot2::geom_point(data = nodedf, ggplot2::aes(x = .data$x, y = .data$y,
      shape = .data$type), size = nodesize, fill = 'white', colour = 'grey30',
      stroke = .6, show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = c(latent = 21, manifest = 22)) +
    ggplot2::geom_text(data = nodedf, ggplot2::aes(x = .data$x, y = .data$y,
      label = .data$name), size = labelsize) +
    ggplot2::scale_linewidth_continuous(range = c(.25, maxwidth), name = '|weight|') +
    ggplot2::scale_colour_manual(values = c(`TRUE` = poscolour, `FALSE` = negcolour),
      labels = c(`TRUE` = 'positive', `FALSE` = 'negative'), name = 'sign') +
    ggplot2::facet_wrap(~ panel, ncol = ncol) +
    ggplot2::coord_equal(clip = 'off') +
    ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(6, 4), 'pt'))
}

#' Draw the temporal and contemporaneous networks of a continuous time model
#'
#' @description
#' Plots what \code{\link{ctNetwork}} computes, in the style of the network
#' psychometrics literature: signed edges, width by magnitude, arrowheads and
#' self-loops on the directed temporal network, and no arrowheads on the
#' undirected contemporaneous one.
#'
#' Give \code{dt} more than one value and the same networks are drawn side by
#' side at each interval, which is the continuous time point: the graph a study
#' would have recovered depends on how often it measured.
#'
#' @details
#' Drawn with ggplot2, which ctsem already depends on, so no extra package is
#' needed. \code{engine='qgraph'} hands the weight matrix to
#' \code{qgraph::qgraph()} instead, if that package is installed --
#' \code{\link{ctNetwork}}'s matrices are in the shape qgraph expects, so
#' \code{qgraph::qgraph(ctNetwork(fit)$temporal, directed=TRUE)} needs nothing
#' from this function.
#'
#' @param x A \code{ctStanModel}, a fit, or a \code{ctNetwork} from
#'   \code{\link{ctNetwork}}.
#' @param dt One or more time intervals. Ignored, with a message, when \code{x}
#'   is already a \code{ctNetwork}.
#' @param networks Which networks to draw; see \code{\link{ctNetwork}}. Each
#'   becomes a panel, as does each \code{dt}.
#' @param threshold Edges with \code{abs(weight)} at or below this are not drawn.
#' @param engine \code{'ggplot'} (the default) or \code{'qgraph'}.
#' @param maxwidth Line width of the largest edge drawn.
#' @param arrowsize Arrowhead length in mm.
#' @param nodesize,labelsize Node and label size.
#' @param poscolour,negcolour Colours for positive and negative edges. The
#'   defaults are a colour-blind safe blue and red; the literature's green and
#'   red are \code{'darkgreen'} and \code{'red'}.
#' @param title Panel title, or NULL for none. \code{'auto'} names the source
#'   and the interval.
#' @param ... Passed to \code{\link{ctNetwork}}, or to \code{qgraph::qgraph()}
#'   when \code{engine='qgraph'}.
#'
#' @return A ggplot2 object, which can be modified further or printed. With
#'   \code{engine='qgraph'}, the qgraph object, invisibly, having drawn.
#'
#' @examples
#' m <- ctModel(type='ct', n.latent=3, n.manifest=3,
#'   manifestNames=c('y1','y2','y3'), latentNames=c('a','b','c'),
#'   LAMBDA=diag(3),
#'   DRIFT=matrix(c(-.4,0,0, .3,-.5,0, 0,.4,-.6), 3, 3, byrow=TRUE),
#'   DIFFUSION=matrix(c(1,0,0, .3,1,0, 0,-.4,1), 3, 3, byrow=TRUE))
#'
#' # The indirect path a -> b -> c is invisible at a short interval and present
#' # at a long one -- the same model, two networks.
#' ctNetworkPlot(m, dt=c(.2, 2), networks='temporal')
#'
#' @export
ctNetworkPlot <- function(x, dt = 1, networks = c('temporal', 'contemporaneous'),
  threshold = 0, engine = c('ggplot', 'qgraph'), maxwidth = 3, arrowsize = 2.6,
  nodesize = 10, labelsize = 3.2, poscolour = '#2166AC', negcolour = '#B2182B',
  title = 'auto', ...){

  engine <- match.arg(engine)

  if(inherits(x, 'ctNetwork')){
    if(!missing(dt)) message('dt is ignored when a ctNetwork object is supplied; ',
      'it was built at dt = ', format(attributes(x)$dt), '.')
    nets <- list(x)
    dt <- attributes(x)$dt
    if(missing(networks)) networks <- attributes(x)$networks
    if(missing(threshold)) threshold <- attributes(x)$threshold
  } else {
    if(!length(dt) || any(!is.finite(dt)) || any(dt <= 0)) stop(call. = FALSE,
      'dt must be positive.')
    nets <- lapply(seq_along(dt), function(i) ctNetwork(x, dt = dt[i],
      networks = networks, threshold = threshold, quiet = i > 1, ...))
  }
  networks <- match.arg(networks,
    c('temporal', 'contemporaneous', 'asymptotic', 'measurement'), several.ok = TRUE)

  if(engine == 'qgraph'){
    if(!requireNamespace('qgraph', quietly = TRUE)) stop(call. = FALSE,
      paste0("engine='qgraph' needs the qgraph package, which ctsem does not ",
        "depend on: install.packages('qgraph'). The default engine='ggplot' ",
        "needs nothing extra, and ctNetwork(x)$temporal can be passed to ",
        "qgraph::qgraph() directly."))
    mat <- nets[[1]][[networks[1]]]
    if(is.null(mat)) stop(call. = FALSE, 'The ', networks[1],
      ' network is not available for this model.')
    if(length(nets) > 1 || length(networks) > 1) message(
      "engine='qgraph' draws one network; showing ", networks[1], ' at dt = ',
      format(dt[1]), '. Call it once per panel, or use engine="ggplot".')
    return(invisible(qgraph::qgraph(mat, directed = networks[1] %in%
        c('temporal', 'measurement'), diag = networks[1] %in% 'temporal',
      labels = colnames(mat), ...)))
  }

  edges <- do.call(rbind, lapply(seq_along(nets), function(i){
    e <- nets[[i]]$edges
    e <- e[e$network %in% networks & abs(e$weight) > threshold, , drop = FALSE]
    if(!nrow(e)) return(NULL)
    e$panel <- if(length(nets) > 1) paste0(e$network, ', dt = ', format(dt[i])) else
      e$network
    e
  }))
  if(is.null(edges) || !nrow(edges)) stop(call. = FALSE,
    'No edges past threshold = ', format(threshold), ', so there is nothing to draw.')
  # Panels in the order asked for, not alphabetically.
  wanted <- unlist(lapply(seq_along(nets), function(i) vapply(networks, function(n)
    if(length(nets) > 1) paste0(n, ', dt = ', format(dt[i])) else n, character(1))))
  edges$panel <- factor(edges$panel, levels = wanted[wanted %in% edges$panel])

  nodes <- nets[[1]]$nodes
  nodes <- nodes[nodes$type == 'latent' | 'measurement' %in% networks, , drop = FALSE]

  if(!is.null(title) && identical(title, 'auto')) title <- paste0(
    'Model implied networks from a ', attributes(nets[[1]])$source,
    if(length(nets) == 1) paste0(', dt = ', format(dt)) else '')

  # Panels are ordered interval-major, so one row per interval means as many
  # columns as networks -- unless there is only one network, when the intervals
  # themselves should run across the page.
  .ctNetworkGG(edges, nodes, lambda = nets[[1]]$measurement, maxwidth = maxwidth,
    title = title, arrowsize = arrowsize, nodesize = nodesize,
    labelsize = labelsize, poscolour = poscolour, negcolour = negcolour,
    ncol = if(length(networks) > 1) length(networks) else length(nets))
}
