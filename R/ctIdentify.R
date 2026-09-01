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
  list(information = scaled, uninformed = flat, nsubjects = nrow(scores))
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
#' only parameters implicated at every one of them are reported as uninformed.
#'
#' @param datalong Long format data, as for \code{\link{ctFit}}.
#' @param ctstanmodel Model from \code{\link{ctModel}}.
#' @param inits Optional evaluation point. Supplying one uses that single point
#'   and disables the comparison across points described above.
#' @param nstart Number of points to evaluate. The first is the origin of the
#'   unconstrained space, so a result is reproducible; the rest are dispersed
#'   draws around it.
#' @param spread Standard deviation of the additional evaluation points on the
#'   raw scale. Large enough that they are genuinely different points, which
#'   \code{\link{ctFit}}'s own \code{rnorm(npar, 0, 0.01)} initialisation is
#'   not.
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
#'   implicated at every evaluation point), \code{$sometimes} (implicated at
#'   some), \code{$nweak}, \code{$condition}, and \code{$starts} holding the
#'   per-point detail. Printing it summarises what was found.
#'
#' @details Requires the Julia backend, which is where the per-subject scores
#'   come from. The statistic has rank at most the number of subjects, so a
#'   model with more free parameters than subjects reports flat directions for
#'   that reason alone; the printed summary says when that applies.
#'
#' @seealso \code{\link{ctFit}} runs the equivalent check on a finished fit,
#'   evaluated at the estimate and stored as \code{fit$identifiability}.
#'
#' @export
ctIdentify <- function(datalong, ctstanmodel, inits = NULL, nstart = 3L,
  spread = 0.5, priors = FALSE, intoverpop = "augmented", cores = 1L,
  verbose = 0L, rtol = 1e-13) {

  spec <- ctFitJuliaBackend(datalong, ctstanmodel, fit = FALSE,
    priors = priors, intoverpop = intoverpop, cores = cores, verbose = verbose)
  npar <- .ctBackendNpar(spec)
  if (!is.finite(npar) || npar < 1L) {
    stop("This model has no free parameters, so there is nothing to identify.",
      call. = FALSE)
  }
  parnames <- .ctBackendRawParameterNames(list(model_spec = spec), npar)

  # The origin first: it is the one point that is the same on every run, so the
  # result is reproducible rather than a function of the seed.
  points <- list(numeric(npar))
  if (!is.null(inits)) {
    points <- list(.ctJuliaInitialValues(npar, inits))
  } else if (nstart > 1L) {
    points <- c(points, lapply(seq_len(as.integer(nstart) - 1L),
      function(i) stats::rnorm(npar, 0, spread)))
  }

  nsubjects <- NA_integer_
  starts <- lapply(points, function(at) {
    info <- .ctIdentifyInformation(spec, at)
    if (is.null(info)) return(NULL)
    nsubjects <<- info$nsubjects
    # `.ctBackendIdentifiability()` takes a Hessian and negates it, so the
    # information is passed negated to arrive the right way up.
    result <- .ctBackendIdentifiability(-info$information, parnames, rtol = rtol)
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

  implicated <- lapply(starts, function(s) s$parameters)
  # Implicated everywhere, which is the structural statement, against
  # implicated somewhere, which is not one.
  always <- Reduce(intersect, implicated)
  ever <- Reduce(union, implicated)
  structure(list(
    parameters = always,
    sometimes = setdiff(ever, always),
    nweak = min(vapply(starts, function(s) as.integer(s$nweak), integer(1))),
    nweakmax = max(vapply(starts, function(s) as.integer(s$nweak), integer(1))),
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
    cat("  Parameters involved: ", paste(x$parameters, collapse = ", "), "\n",
      sep = "")
    cat("  These are not estimable from this data as the model stands. Fix one\n",
      "  of each set to a value, or remove it.\n", sep = "")
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
