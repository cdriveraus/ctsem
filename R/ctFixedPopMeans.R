# Stated population means: ctPar(mean=), 'name=value|...' in a cell, or
# `model$fixedpopmeans` set directly.
#
# A parameter with random effects normally has its population mean estimated
# beside the effects' spread. Stating the mean instead -- zero for a
# person-level latent level, say -- leaves only the spread to fit. On the julia
# route the stated value has no raw coordinate: the optimiser, the samplers, the
# Hessian and every raw vector a fit reports cover the free parameters alone,
# and the engine appends the stated values after them (`CTSEMAnchors`)
# wherever a subject's parameters are built, so a cell reads one and the random
# effects shift it like any other raw parameter.
#
# `.ctJuliaPrepare()` lays the raw vector out as it would for an estimated mean
# and then moves each stated one past the end (`.ctJuliaAnchorLayout()`),
# renumbering every index that names a raw position. A reader of a position in
# the specification then finds it where the engine does, and `spec$anchors`
# says what follows the free parameters. R code that indexes a raw vector by
# such a position extends the vector first (`.ctBackendRawFull()`); positions
# up to the free count mean what they always did.

# The model's stated means, checked, with the raw value each corresponds to:
# one row a parameter, `param`, `value` on the scale summary() reports
# population means on, and `raw`, the value its transform maps there. NULL
# when there are none.
#' @keywords internal
.ctFixedPopMeans <- function(model) {
  means <- model[['fixedpopmeans']]
  if (!length(means)) return(NULL)
  if (!is.numeric(means) || is.null(names(means)) || any(!nzchar(names(means))) ||
      any(!is.finite(means)) || anyDuplicated(names(means))) {
    stop("fixedpopmeans must be a named numeric vector: one finite population ",
      "mean per parameter name.", call. = FALSE)
  }
  pars <- model$pars
  unknown <- setdiff(names(means), pars$param)
  if (length(unknown)) {
    stop("fixedpopmeans names ", paste(unknown, collapse = ", "), ", which ",
      if (length(unknown) > 1L) "are not parameters" else "is not a parameter",
      " of this model.", call. = FALSE)
  }
  varying <- .ctVaryingParams(model)
  raw <- numeric(length(means))
  for (i in seq_along(means)) {
    name <- names(means)[i]
    row <- match(name, pars$param)
    transform <- pars$transform[row]
    if (is.na(transform) || !nzchar(transform)) {
      stop("A population mean is stated for ", name, ", which is not a free ",
        "parameter with a transform of its own.", call. = FALSE)
    }
    if (!name %in% varying) {
      stop("A population mean is stated for ", name, ", which does not vary ",
        "(indvarying). A stated mean is for a parameter with random effects ",
        "about it; to fix the parameter itself, write its value in the cell.",
        call. = FALSE)
    }
    raw[i] <- .ctGenerateRawFor(transform, means[[i]],
      paste0("the population mean of ", name))
  }
  data.frame(param = names(means), value = as.numeric(means), raw = raw,
    stringsAsFactors = FALSE)
}

# The specification's pieces with every stated mean's raw position moved past
# the end of the vector, and `anchors` describing them. `nfull` is the
# vector's length with them in place. Free positions keep their order, so the
# renumbering closes the gaps the stated means leave.
#' @keywords internal
.ctJuliaAnchorLayout <- function(means, nfull, parameter_table, laplace,
  ti_effects, ti_fixed, ti_missing, random_effects, priors) {
  table <- as.data.frame(parameter_table, stringsAsFactors = FALSE)
  old <- as.integer(table$parnumber[match(means$param, table$param)])
  if (anyNA(old) || any(old < 1L) || any(old > nfull)) {
    stop("Internal layout error: a stated population mean has no raw position. ",
      "This is a bug in ctsem rather than in the model; please report it.",
      call. = FALSE)
  }
  free <- setdiff(seq_len(nfull), old)
  map <- integer(nfull)
  map[free] <- seq_along(free)
  map[old] <- length(free) + seq_along(old)
  remap <- function(x) {
    if (is.null(x)) return(x)
    x <- as.integer(x)
    use <- !is.na(x) & x > 0L
    x[use] <- map[x[use]]
    x
  }
  # The rendered transforms read raw positions by number, `param[k]`.
  retext <- function(x) {
    if (is.null(x)) return(x)
    out <- as.character(x)
    for (i in which(!is.na(out) & grepl("param[", out, fixed = TRUE))) {
      at <- gregexpr("\\bparam\\[([0-9]+)\\]", out[i], perl = TRUE)
      hits <- regmatches(out[i], at)[[1]]
      if (!length(hits)) next
      k <- as.integer(sub("^param\\[([0-9]+)\\]$", "\\1", hits))
      regmatches(out[i], at) <- list(paste0("param[", map[k], "]"))
    }
    out
  }
  parameter_table$parnumber <- remap(parameter_table$parnumber)
  for (column in c("transform", "predicttransform", "updatetransform",
      "tdtransform")) {
    if (!is.null(parameter_table[[column]])) {
      parameter_table[[column]] <- retext(parameter_table[[column]])
    }
  }
  if (!is.null(laplace)) {
    fields <- c("re_index", "sd_index", "cor_index", "load_index", "unique_index")
    for (l in seq_along(laplace$levels)) {
      for (f in fields) {
        if (!is.null(laplace$levels[[l]][[f]])) {
          laplace$levels[[l]][[f]] <- remap(laplace$levels[[l]][[f]])
        }
      }
    }
    for (f in c("re_index", "sd_index", "cor_index")) {
      if (!is.null(laplace[[f]])) laplace[[f]] <- remap(laplace[[f]])
    }
  }
  if (!is.null(ti_effects) && nrow(ti_effects)) {
    ti_effects$parameter <- remap(ti_effects$parameter)
    ti_effects$coefficient <- remap(ti_effects$coefficient)
  }
  if (!is.null(ti_fixed) && nrow(ti_fixed)) ti_fixed$parameter <- remap(ti_fixed$parameter)
  if (!is.null(ti_missing) && nrow(ti_missing)) {
    ti_missing$parameter <- remap(ti_missing$parameter)
  }
  if (!is.null(random_effects) && length(random_effects) && nrow(random_effects)) {
    random_effects$parameter <- remap(random_effects$parameter)
  }
  # A stated value takes no prior: it is not a parameter.
  if (!is.null(priors)) {
    keep <- !as.integer(priors$index) %in% old
    priors$index <- remap(priors$index[keep])
    priors$scale <- priors$scale[keep]
  }
  list(parameter_table = parameter_table, laplace = laplace,
    ti_effects = ti_effects, ti_fixed = ti_fixed, ti_missing = ti_missing,
    random_effects = random_effects, priors = priors,
    anchors = list(param = means$param, value = means$value, raw = means$raw,
      index = as.integer(length(free) + seq_along(old))))
}

# How many stated population means a fit or specification carries.
#' @keywords internal
.ctBackendNanchors <- function(x) {
  spec <- if (!is.null(x$model_spec)) x$model_spec else x
  length(spec$anchors$index)
}

# `raw` -- a vector of the free parameters, or a matrix of draws of them, one
# row a draw -- with the stated population means after them, for R code that
# indexes it by a position in the specification. Returned as it is for a model
# with none, or when it already carries them.
#' @keywords internal
.ctBackendRawFull <- function(x, raw) {
  spec <- if (!is.null(x$model_spec)) x$model_spec else x
  values <- as.numeric(spec$anchors$raw)
  na <- length(values)
  if (!na) return(raw)
  nfree <- min(spec$anchors$index) - 1L
  width <- if (is.matrix(raw)) ncol(raw) else length(raw)
  if (width == nfree + na) return(raw)
  if (width != nfree) {
    stop("Internal layout error: a raw vector of ", width, " for a model with ",
      nfree, " free parameters and ", na, " stated population means. This is ",
      "a bug in ctsem; please report it.", call. = FALSE)
  }
  if (is.matrix(raw)) {
    return(cbind(raw, matrix(values, nrow(raw), na, byrow = TRUE,
      dimnames = list(NULL, spec$anchors$param))))
  }
  c(raw, values)
}
