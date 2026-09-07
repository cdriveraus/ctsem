# Which parameters the data can inform, before spending a fit finding out.
#
# `ctModel()` checks that a specification is well formed. Nothing checks whether
# its free parameters can be estimated from a particular dataset, because that
# is a property of the model and the data together. So a specification can pass
# every structural check and still be unidentified, and the way that surfaces is
# an optimiser that will not settle or intervals that span everything -- long
# after the point where the user could have changed the model cheaply.
#
# A direction of the parameter space the likelihood is flat in is a direction
# the data does not determine, and its eigenvector says which parameters are
# involved. `.ctBackendIdentifiability()` reads that off the Hessian of a
# finished fit, which is the right object *at a maximum*.
#
# It is the wrong object here, and measurably so. Away from the mode the Hessian
# is indefinite: the curvature of the log likelihood at an arbitrary point mixes
# the information with a residual term carrying nothing about identification.
# Evaluated at starting values on a plainly identified one-latent model it
# reported diff_eta1 and mvarY1 as undetermined and found a direction of
# negative curvature -- both facts about standing away from the maximum rather
# than about the data.
#
# The outer product of the per-subject scores is the right object. It is
# positive semi-definite wherever it is evaluated, so there is no spurious
# negative curvature to explain away, and its null space says exactly the thing
# worth saying: a direction in which *no subject's* log likelihood changes is a
# direction no amount of this data can distinguish. One traced pass computes it,
# which is cheaper than a Hessian as well as more meaningful.

# Evaluated away from an optimum, which is the whole point and also the
# limitation, so it is worth being exact about what carries.
#
# Structural non-identification -- two parameters that only ever appear as a
# product, a latent nothing loads on, a scale fixed nowhere -- is a rank
# deficiency that holds at every point in the parameter space, so an arbitrary
# point finds it as well as an optimum would. That is the case this catches, and
# it is the case worth catching early.
#
# Weak identification is not like that. A likelihood can be flat somewhere near
# the start and perfectly curved at the mode, or the reverse. A single point
# could therefore report a flat direction that a fit would not have, which would
# be worse than not checking -- a pre-fit warning people learn to ignore. So
# several points are evaluated, and a parameter is only reported as uninformed
# if it is implicated at *all* of them; the ones implicated at some are reported
# separately, as the weaker statement it is.
#
# Those points have to be genuinely different, which is why they are not
# `.ctJuliaInitialValues()`. That draws from rnorm(npar, 0, 0.01) -- correct for
# starting an optimiser near the origin, useless here, because three such draws
# are three names for the same point and agreement between them says nothing.
#
# One thing this cannot see past: the outer product of n subject scores has rank
# at most n, so a model with more free parameters than subjects has flat
# directions for that reason alone. It is reported, because it is true and worth
# knowing, but it is a statement about the sample size rather than about the
# specification.

# The information the data carries about each direction at a given point.
#' @keywords internal
.ctIdentifyInformation <- function(spec, at) {
  # The engine call directly rather than `.ctBackendScoreMatrix()`, which
  # expects a fitted object: here there is a prepared model and no fit, and
  # wrapping one to look like the other would be the more fragile of the two.
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  result <- try(JuliaConnectoR::juliaGet(module$ctsem_subject_gradients(
    .ctJuliaObjective(handle), .ctJuliaVector(as.numeric(at)))), silent = TRUE)
  if (inherits(result, "try-error") || is.null(result$scores)) return(NULL)
  scores <- as.matrix(result$scores)
  if (any(!is.finite(scores))) return(NULL)
  # Whichever axis is not the parameter axis is the subject axis; orienting on
  # the parameter count rather than on an assumption keeps this correct if the
  # engine's layout changes.
  if (nrow(scores) == length(at) && ncol(scores) != length(at)) {
    scores <- t(scores)
  }
  if (ncol(scores) != length(at)) return(NULL)
  information <- crossprod(scores)

  # Scaled to unit diagonal before anything reads its eigenvalues.
  #
  # The raw parameters are on incomparable scales -- each has its own transform,
  # so a unit step in raw `lambda` and a unit step in raw `T0var` move the log
  # likelihood by different orders of magnitude. Unscaled, the information
  # matrix of a perfectly well identified one-latent model has a condition
  # number of 3.3e7, and a relative-eigenvalue threshold applied to that is
  # measuring the units rather than the identification.
  #
  # In correlation form the question becomes the right one: is this direction
  # uninformed *relative to how well its parameters are informed individually*.
  # A parameter whose own diagonal entry is zero is a separate and stronger
  # finding -- the data says nothing about it at all -- and is reported as
  # such rather than being hidden by the division.
  diagonal <- diag(information)
  flat <- !is.finite(diagonal) | diagonal <= 0
  scale <- sqrt(ifelse(flat, 1, diagonal))
  scaled <- information / outer(scale, scale)
  # `scale` goes out with it: an eigenvector of the scaled matrix is a
  # direction in scaled coordinates, and anything comparing it against a
  # gradient taken in raw ones has to divide by this first.
  list(information = scaled, uninformed = flat, nsubjects = nrow(scores),
    metric = scale)
}

# The evaluation points, all of them reproducible.
#
# `nstart - 1` dispersed draws used to come straight from the session stream,
# so two identical calls could report different parameters -- observed:
# `always = (empty)` from one and `always = popsd_df11` from the next. Only the
# first point, the origin, was reproducible, and the docstring claimed no more
# than that. A named seed fixes the whole function, and the session stream is
# put back afterwards so that calling this neither depends on where the user's
# RNG had got to nor moves it.
#' @keywords internal
.ctIdentifyPoints <- function(npar, nstart, spread, inits, seed) {
  draw <- function() {
    if (!is.null(inits)) return(list(.ctJuliaInitialValues(npar, inits)))
    c(list(numeric(npar)),
      if (nstart > 1L) lapply(seq_len(as.integer(nstart) - 1L),
        function(i) stats::rnorm(npar, 0, spread)) else NULL)
  }
  if (is.null(seed)) return(draw())
  existing <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else NULL
  set.seed(as.integer(seed)[1L])
  on.exit({
    if (is.null(existing)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else assign(".Random.seed", existing, envir = globalenv())
  }, add = TRUE)
  draw()
}

# Which coordinates lie in the flat subspace, rather than on one of the axes an
# eigendecomposition happened to return.
#
# Intersecting parameter *names* across evaluation points is what hid the model
# this function exists for. On it there is a flat direction at every point --
# relative eigenvalues 1.8e-19, 1.2e-16, 3.9e-18 -- but the loading rotates
# along the ridge, so the names implicated at one point are not the names
# implicated at the next, the intersection came out empty, and everything
# landed in `$sometimes` under text saying a fit may still settle. It will not.
#
# The invariant is the *dimension* of the flat subspace, not any basis for it.
# So `k`, the smallest number of flat directions seen at any point, is the rank
# deficiency the data has everywhere, and the `k` flattest directions at each
# point are a basis for it. A coordinate's involvement is then the length of
# its projection onto that subspace, which is what an orthonormal basis gives
# as the norm of its row -- unchanged by any rotation within the subspace.
# Directions beyond the first `k` at a point that had more are the weaker
# statement, and are reported as one.
#' @keywords internal
.ctIdentifySubspaceLoading <- function(directions, npar, loading = 0.25) {
  if (!length(directions)) return(logical(npar))
  vectors <- vapply(directions, function(d) {
    v <- d$vector
    if (length(v) != npar) rep(NA_real_, npar) else v
  }, numeric(npar))
  vectors <- matrix(vectors, nrow = npar)
  if (any(!is.finite(vectors))) return(logical(npar))
  sqrt(rowSums(vectors^2)) >= loading
}

#' Check which parameters a dataset can inform, before fitting
#'
#' Reports the directions of the parameter space this data carries no
#' information about, and the parameters involved in them, without fitting the
#' model.
#'
#' This catches structural non-identification -- parameters that only appear
#' together, a latent process nothing measures, a scale that is never fixed --
#' which is a property of the model and holds wherever it is evaluated. It is
#' not a verdict on the model: a direction that is flat at one point may be
#' perfectly informed at the estimate, which is why several points are used and
#' only parameters lying in the flat subspace at every one of them are reported
#' as uninformed.
#'
#' Not every flat direction means the same thing. A random effect on a variance
#' cell (DIFFUSION, MANIFESTVAR) under \code{intoverpop = 'augmented'} is
#' partially identified: the population covariances it generates with the other
#' individually varying parameters are determined, while their decomposition
#' into a scale and correlations is a ridge. The printed summary distinguishes
#' the two cases, because the advice differs -- fixing a value discards
#' information in the first case and is the fix in the second.
#'
#' @param datalong Long format data, as for \code{\link{ctFit}}.
#' @param model Model from \code{\link{ctModel}}.
#' @param ctstanmodel Deprecated. Use \code{model}.
#' @param inits Optional evaluation point. Supplying one uses that single point
#'   and disables the comparison across points described above.
#' @param nstart Number of points to evaluate. The first is the origin of the
#'   unconstrained space; the rest are dispersed draws around it.
#' @param spread Standard deviation of the additional evaluation points on the
#'   raw scale. Large enough that they are genuinely different points, which
#'   \code{\link{ctFit}}'s own \code{rnorm(npar, 0, 0.01)} initialisation is
#'   not.
#' @param seed Seed for the dispersed evaluation points, so that two identical
#'   calls give identical output. The session's random stream is restored
#'   afterwards, so this neither depends on it nor disturbs it. \code{NULL}
#'   draws from the session stream instead, which makes the result depend on
#'   where that had got to.
#' @param priors Whether to include the prior. \code{FALSE} by default, because
#'   the question is what the \emph{data} can inform: a prior informs every
#'   direction, so a model unidentified by its data looks fine with priors on
#'   and the check stops being able to say anything.
#' @param intoverpop Integration approach, as for \code{\link{ctFit}}.
#' @param cores Passed to the engine.
#' @param verbose Passed to the engine.
#' @param rtol A direction counts as uninformed when the information along it
#'   falls below this fraction of the best-informed direction's. The default
#'   separates zero from merely small: on the scaled information a structurally
#'   unidentified direction sits at machine precision (measured: 7e-17), while
#'   a weakly identified one on a comparable model sat at 3e-12. Raise it to
#'   catch weak directions too, at the cost of reporting some a fit would have
#'   settled.
#'
#' @return An object of class \code{ctIdentify} with \code{$parameters} (names
#'   in the flat subspace at every evaluation point), \code{$partial} and
#'   \code{$structural} splitting those into the partially and completely
#'   unidentified, \code{$sometimes} (in a flat direction beyond that subspace
#'   at some points), \code{$rotating} (whether the basis of the flat subspace
#'   turns between points), \code{$nweak}, \code{$condition}, and
#'   \code{$starts} holding the per-point detail. Printing it summarises what
#'   was found.
#'
#' @details Requires the Julia backend, which is where the per-subject scores
#'   come from. The statistic has rank at most the number of subjects, so a
#'   model with more free parameters than subjects reports flat directions for
#'   that reason alone; the printed summary says when that applies.
#'
#' @seealso \code{\link{ctFit}} runs the equivalent check on a finished fit,
#'   evaluated at the estimate and stored as \code{fit$identifiability}.
#'
#' @examples
#' \donttest{
#' # ctIdentify() always uses the julia backend internally, whatever the
#' # eventual fit will use, so the example is inert where Julia is absent --
#' # including on CRAN, whose check machines have none.
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   gen <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
#'     manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1),
#'     DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
#'     T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
#'     MANIFESTMEANS = matrix(0), Tpoints = 8))
#'   datalong <- ctGenerate(gen, n.subjects = 40, Tpoints = 8, backend = 'r')
#'
#'   model <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
#'     manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1),
#'     T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0)))
#'
#'   result <- ctIdentify(datalong, model, cores = 1)
#'   print(result)
#' }
#' }
#'
#' @export
ctIdentify <- function(datalong, model, inits = NULL, nstart = 3L,
  spread = 0.5, priors = FALSE, intoverpop = "augmented", cores = 1L,
  verbose = 0L, rtol = 1e-13, seed = 1L, ctstanmodel) {

  if(missing(model)){
    if(missing(ctstanmodel)) stop('model must be supplied')
    warning('ctstanmodel argument is deprecated, use model instead')
    model <- ctstanmodel
  } else if(!missing(ctstanmodel)) {
    stop('Use only one of model or deprecated ctstanmodel')
  }
  ctstanmodel <- model

  spec <- .ctFitJuliaBackend(datalong, ctstanmodel, fit = FALSE,
    priors = priors, intoverpop = intoverpop, cores = cores, verbose = verbose)
  npar <- .ctBackendNpar(spec)
  if (!is.finite(npar) || npar < 1L) {
    stop("This model has no free parameters, so there is nothing to identify.",
      call. = FALSE)
  }
  parnames <- .ctBackendRawParameterNames(list(model_spec = spec), npar)

  # The origin first, then the dispersed points; see `.ctIdentifyPoints()` for
  # why the whole set is seeded rather than only the first being fixed.
  points <- .ctIdentifyPoints(npar, nstart, spread, inits, seed)

  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  nsubjects <- NA_integer_
  starts <- lapply(points, function(at) {
    info <- .ctIdentifyInformation(spec, at)
    if (is.null(info)) return(NULL)
    nsubjects <<- info$nsubjects
    # `.ctBackendIdentifiability()` takes a Hessian and negates it, so the
    # information is passed negated to arrive the right way up. `metric` is
    # what the information was scaled by, so the partial-identification check
    # can bring a scaled eigenvector back to raw coordinates.
    result <- .ctBackendIdentifiability(-info$information, parnames, rtol = rtol,
      fit = handle, at = at, metric = info$metric, vectors = TRUE)
    result$at <- at
    result$uninformed <- parnames[info$uninformed]
    spectrum <- eigen(info$information, symmetric = TRUE,
      only.values = TRUE)$values
    result$relative <- sort(spectrum / max(spectrum))
    result
  })
  usable <- !vapply(starts, is.null, logical(1))
  if (!any(usable)) {
    stop("The engine could not compute per-subject scores for this model at ",
      "any evaluation point, so identifiability cannot be assessed here.",
      call. = FALSE)
  }
  starts <- starts[usable]

  # Aggregated by subspace, not by name. `k` is the rank deficiency the data
  # has at every point, and the `k` flattest directions at each point span it;
  # see `.ctIdentifySubspaceLoading()` for what intersecting names did instead.
  nweak <- vapply(starts, function(s) as.integer(s$nweak), integer(1))
  k <- min(nweak)
  flattest <- lapply(starts, function(s) {
    if (!length(s$directions)) return(list())
    s$directions[order(vapply(s$directions, function(d) as.numeric(d$relative),
      numeric(1)))]
  })
  core <- lapply(flattest, function(d) d[seq_len(min(k, length(d)))])
  extra <- lapply(flattest, function(d)
    if (length(d) > k) d[seq.int(k + 1L, length(d))] else list())
  ispartial <- function(directions) !vapply(directions,
    function(d) is.null(d$partial), logical(1))
  named <- function(directions) parnames[.ctIdentifySubspaceLoading(directions,
    npar)]
  corenames <- lapply(core, named)
  extranames <- lapply(extra, named)
  ridgenames <- lapply(core, function(d) named(d[ispartial(d)]))
  # Union rather than intersection over the core. Every core direction is flat
  # at the point it came from, so a coordinate lying in one is part of the
  # structural finding there; a coordinate that appears at one point and not
  # another is the ridge turning, which is a fact about the ridge's curvature
  # and not about identification. `rotating` says that happened, because a
  # reader comparing this against the per-point detail in `$starts` will
  # otherwise see two different answers and not know which to trust.
  always <- as.character(sort(unique(unlist(corenames))))
  rotating <- length(always) > 0L &&
    !all(vapply(corenames, function(n) setequal(n, always), logical(1)))
  partition <- .ctIdentifyPartition(unlist(core, recursive = FALSE))
  ridge <- as.character(sort(unique(unlist(ridgenames))))
  structure(list(
    parameters = always,
    # Model parameter names, for the sentence about population sds; the raw
    # coordinates that ridge occupies are `$ridge`.
    partial = partition$partial,
    partners = partition$partners,
    ridge = ridge,
    structural = setdiff(always, ridge),
    rotating = rotating,
    sometimes = setdiff(as.character(sort(unique(unlist(extranames)))), always),
    nweak = k,
    nweakmax = max(nweak),
    condition = stats::median(vapply(starts, function(s)
      as.numeric(s$condition), numeric(1))),
    smallest = max(vapply(starts, function(s) s$relative[1L], numeric(1))),
    npar = npar, parnames = parnames, nstart = length(starts),
    nsubjects = nsubjects, rankLimited = isTRUE(nsubjects < npar),
    priors = priors, starts = starts), class = "ctIdentify")
}

#' @export
print.ctIdentify <- function(x, ...) {
  cat("Identifiability before fitting\n")
  cat("  ", x$npar, " free parameters, ", x$nsubjects, " subjects, ",
    x$nstart, " evaluation point", if (x$nstart > 1L) "s" else "",
    if (x$priors) ", prior included" else ", data only", "\n", sep = "")
  # The number, not just the verdict. A threshold has to fall somewhere, and a
  # reader who can see 7e-17 against 3e-12 can tell a structural flat direction
  # from a merely weak one without trusting where this function put the line.
  cat("  Weakest direction carries ", signif(x$smallest, 2),
    " of the information of the strongest.\n", sep = "")
  if (!length(x$parameters) && !length(x$sometimes)) {
    cat("  No uninformed directions found: this data carries information about\n",
      "  every direction of this model. That is what identification looks like\n",
      "  here, and it is not a guarantee about the estimate.\n", sep = "")
    return(invisible(x))
  }
  if (length(x$parameters)) {
    cat("  ", x$nweak, " direction", if (x$nweak != 1L) "s" else "",
      " the data carries no information about, at every point checked.\n",
      sep = "")
    if (isTRUE(x$rotating)) {
      cat("  It turns between points, so what follows names every parameter\n",
        "  lying in it somewhere rather than only those it involves at each.\n",
        "  See $starts for the per-point detail.\n", sep = "")
    }
    # Partially and completely unidentified read the same in the eigenvalues
    # and need opposite advice, so `.ctIdentifyAdvice()` -- the same wording a
    # finished fit warns with -- says which is which. The fallback is for a
    # partition that classifies everything and so has nothing left to list:
    # printing nothing here would read as a clean result.
    advice <- .ctIdentifyAdvice(x)
    if (!length(advice)) advice <- paste0("Parameters involved: ",
      paste(x$parameters, collapse = ", "), ".")
    writeLines(strwrap(advice, indent = 2, exdent = 2, width = 78))
  }
  if (length(x$sometimes)) {
    cat("  Uninformed at some points but not all: ",
      paste(x$sometimes, collapse = ", "), "\n", sep = "")
    cat("  Weakly rather than structurally determined; a fit may still settle.\n",
      sep = "")
  }
  if (isTRUE(x$rankLimited)) {
    cat("  Note: ", x$nsubjects, " subjects against ", x$npar,
      " free parameters. At least ", x$npar - x$nsubjects, " direction",
      if (x$npar - x$nsubjects != 1L) "s follow" else " follows",
      " from\n  that alone, whatever the specification says.\n", sep = "")
  }
  invisible(x)
}
