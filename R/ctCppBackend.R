# Hand-written C++ likelihood backend --------------------------------------
#
# The engine lives in `inst/include/ctsemcpp/` and is compiled into ctsem's own
# shared library, so `backend='cpp'` needs no external toolchain at fit time and
# no per-model Stan compile.
#
# Deliberately, this file is thin. The C++ engine consumes *exactly* the model
# specification `.ctJuliaPrepare()` already builds -- one row per model-matrix
# cell, with each cell's transform rendered as an arithmetic expression string
# -- and interprets those strings at runtime rather than compiling them. That
# is the whole reason a third backend costs so little R-side code: the
# expensive, error-prone part (turning a `ctModel` into a canonical, augmented
# parameter table) is already shared with the Julia backend, and re-deriving it
# is exactly the seam where the Stan/Julia backends have historically drifted
# apart.

.ct_cpp_cache <- new.env(parent = emptyenv())
.ct_cpp_cache$objectives <- new.env(parent = emptyenv())

.ctCppUnsupported <- function(model, optimize, priors, intoverpop, vb, gendata,
  stanmodeltext, compileArgs, forcerecompile) {
  failures <- character()
  if (!isTRUE(optimize)) failures <- c(failures, "optimize=FALSE (HMC)")
  if (any(model$manifesttype > 0)) failures <- c(failures, "non-Gaussian manifest variables")
  if (isTRUE(vb)) failures <- c(failures, "variational Bayes")
  if (isTRUE(gendata)) failures <- c(failures, "generation")
  if (!is.na(stanmodeltext)[1] || length(compileArgs) > 0L || isTRUE(forcerecompile)) {
    failures <- c(failures, "Stan compilation controls")
  }
  if (length(failures)) {
    stop("C++ backend v1 does not support: ", paste(failures, collapse = ", "), ".", call. = FALSE)
  }
}

.ctCppPrepare <- function(datalong, model, prepared_data = NULL, priors = FALSE) {
  # Same canonical specification as backend='julia'; the C++ and Julia engines
  # differ only in what they do with the transform strings (interpret an AST
  # vs. `Meta.parse` + `eval`), not in what they are given.
  spec <- .ctJuliaPrepare(datalong, model, prepared_data = prepared_data,
    project = NULL, priors = priors)
  spec$project <- NULL
  spec$engine <- NULL
  spec$class <- "ctCppModel"
  spec
}

# The cache only saves rebuilding the objective (parsing every transform and
# copying the data), which is milliseconds. Without digest installed the backend
# still works, it just rebuilds; that is a better failure mode than refusing to
# run over a Suggests-level dependency.
.ctCppObjectiveKey <- function(spec) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NULL)
  digest::digest(list(spec$parameter_table, spec$subject_starts, spec$times,
    spec$manifest_data, spec$tdpred_data, spec$tipred_data,
    spec$ti_effects, spec$priors, spec$max_timestep, spec$dynamic_state_indices),
    algo = "sha256")
}

.ctCppSpecForEngine <- function(spec) {
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  for (name in c("transform", "predicttransform", "updatetransform", "tdtransform")) {
    if (is.null(table[[name]])) table[[name]] <- NA_character_
    table[[name]] <- as.character(table[[name]])
  }
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) {
    effects <- data.frame(parameter = integer(), predictor = integer(), coefficient = integer())
  }
  list(
    parameter_table = list(
      matrix = as.character(table$matrix), row = as.integer(table$row),
      col = as.integer(table$col), parnumber = as.integer(table$parnumber),
      value = as.numeric(table$value), transform = table$transform,
      predicttransform = table$predicttransform,
      updatetransform = table$updatetransform, tdtransform = table$tdtransform),
    ti_effects = list(parameter = as.integer(effects$parameter),
      predictor = as.integer(effects$predictor),
      coefficient = as.integer(effects$coefficient)),
    dynamic_state_indices = as.integer(spec$dynamic_state_indices),
    subject_starts = as.integer(spec$subject_starts),
    times = as.numeric(spec$times),
    manifest_data = matrix(as.numeric(spec$manifest_data),
      nrow = nrow(spec$manifest_data), ncol = ncol(spec$manifest_data)),
    tdpred_data = matrix(as.numeric(spec$tdpred_data),
      nrow = nrow(spec$tdpred_data), ncol = ncol(spec$tdpred_data)),
    tipred_data = matrix(as.numeric(spec$tipred_data),
      nrow = nrow(spec$tipred_data), ncol = ncol(spec$tipred_data)),
    max_timestep = as.numeric(spec$max_timestep)[1L],
    continuous_time = as.integer(isTRUE(spec$continuoustime)),
    prior_index = if (is.null(spec$priors)) integer() else as.integer(spec$priors$index),
    prior_scale = if (is.null(spec$priors)) numeric() else as.numeric(spec$priors$scale),
    prior_weight = if (is.null(spec$priors)) 1 else as.numeric(spec$priors$weight)
  )
}

.ctCppObjective <- function(object) {
  stopifnot(inherits(object, "ctCppModel") || inherits(object, "ctCppFit"))
  spec <- if (inherits(object, "ctCppFit")) object$model_spec else object
  key <- .ctCppObjectiveKey(spec)
  # The handle is an external pointer, so it cannot survive saveRDS or a new
  # session; it is rebuilt on demand, exactly as the Julia proxy is.
  if (!is.null(key) && exists(key, envir = .ct_cpp_cache$objectives, inherits = FALSE)) {
    handle <- get(key, envir = .ct_cpp_cache$objectives, inherits = FALSE)
    valid <- tryCatch({ .ctsemCppNpars(handle); TRUE }, error = function(e) FALSE)
    if (valid) return(handle)
  }
  handle <- .ctsemCppBuild(.ctCppSpecForEngine(spec))
  if (!is.null(key)) assign(key, handle, envir = .ct_cpp_cache$objectives)
  handle
}

.ctCppNpar <- function(spec) {
  max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient), na.rm = TRUE)
}

#' Evaluate a prepared C++ ctsem likelihood
#'
#' @param object A \code{ctCppModel} (from \code{ctFit(..., backend='cpp', fit=FALSE)})
#'   or a \code{ctCppFit}.
#' @param pars Unconstrained parameters; defaults to the fitted estimate.
#' @param gradient Return the reverse-mode (adjoint) gradient.
#' @param contributions Also return per-subject log likelihoods.
#' @return A list with \code{value} and, when requested, \code{gradient}.
#' @details The gradient is a hand-written reverse-mode adjoint over doubles and
#'   Eigen: one traced forward filter plus one reverse sweep per subject, so its
#'   cost does not grow with the number of free parameters.
#' @export
ctCppEvaluate <- function(object, pars = NULL, gradient = TRUE, contributions = FALSE) {
  if (!inherits(object, c("ctCppModel", "ctCppFit"))) {
    stop("object must be a ctCppModel or ctCppFit", call. = FALSE)
  }
  if (is.null(pars)) {
    if (inherits(object, "ctCppFit")) pars <- object$estimate$raw
    else stop("pars must be supplied for a prepared ctCppModel", call. = FALSE)
  }
  .ctsemCppEvaluate(.ctCppObjective(object), as.numeric(pars), gradient = isTRUE(gradient),
    contributions = isTRUE(contributions))
}

ctFitCppBackend <- function(datalong, model, prepared_data = NULL, inits = NULL, cores = 1L,
  backendcontrol = list(), optimcontrol = list(), verbose = 0L, fit = TRUE,
  priors = FALSE) {
  model_spec <- .ctCppPrepare(datalong, model, prepared_data = prepared_data,
    priors = priors)
  if (!fit) return(structure(model_spec, class = c("ctCppModel", "ctFitModel")))

  handle <- .ctCppObjective(structure(model_spec, class = c("ctCppModel", "ctFitModel")))
  npar <- .ctCppNpar(model_spec)
  # Match stanoptimis() and the Julia backend: absent initial values are small,
  # R-seeded draws in unconstrained space rather than an exact all-zero vector.
  start <- if (is.null(inits) || identical(inits, "random")) stats::rnorm(npar, 0, .01) else {
    values <- as.numeric(inits)
    if (length(values) != npar || anyNA(values)) {
      stop("C++ initial values must be one finite number per free parameter.", call. = FALSE)
    }
    values
  }
  result <- .ctsemCppOptimize(handle, start,
    maxiter = as.integer(.ctJuliaOr(backendcontrol$maxiter, 1000L)),
    gtol = .ctJuliaOr(backendcontrol$g_tol, 1e-8))

  out <- list(backend = "cpp", model = model, model_spec = model_spec,
    data = datalong, estimate = list(raw = as.numeric(result$minimizer),
      loglik = as.numeric(result$maximum_loglik), gradient = as.numeric(result$gradient),
      subject_loglik = as.numeric(result$subject_loglik),
      converged = isTRUE(result$converged), iterations = as.integer(result$iterations)),
    args = list(backend = "cpp", backendcontrol = backendcontrol, cores = cores))
  class(out) <- c("ctCppFit", "ctFit")
  out
}

#' @export
print.ctCppFit <- function(x, ...) {
  cat("ctsem C++ fit\n")
  cat("  log likelihood:", format(x$estimate$loglik), "\n")
  cat("  converged:", x$estimate$converged, " iterations:", x$estimate$iterations, "\n")
  invisible(x)
}

#' @export
coef.ctCppFit <- function(object, ...) object$estimate$raw

#' @export
logLik.ctCppFit <- function(object, ...) {
  structure(object$estimate$loglik, df = length(object$estimate$raw),
    nobs = nrow(object$data), class = "logLik")
}

#' @export
summary.ctCppFit <- function(object, timeinterval = 1, digits = 3, parmatrices = TRUE, ...) {
  .ctBackendSummary(object, timeinterval = timeinterval, digits = digits,
    parmatrices = parmatrices, ...)
}

#' @export
ctExtract.ctCppFit <- function(object, subjectMatrices = FALSE, cores = 2,
  nsamples = "all", subjects = "all", ...) {
  .ctBackendExtract(object, subjectMatrices = subjectMatrices, nsamples = nsamples,
    subjects = subjects, ...)
}

#' @export
ctSummaryMatrices.ctCppFit <- function(fit, calcfunc = quantile,
  calcfuncargs = list(probs = 0.5), timeinterval = 1, ...) {
  .ctBackendSummaryMatrices(fit, calcfunc = calcfunc, calcfuncargs = calcfuncargs,
    timeinterval = timeinterval, ...)
}
