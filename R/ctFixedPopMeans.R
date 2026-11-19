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

# The model's stated means, checked: one row a parameter, `param`, its
# `value` on the scale summary() reports population means on, and `raw`, the
# value its transform maps there. NULL when there are none.
#
# The transform is the cell's own, text as ctStanModel() writes it. The model
# setup (`.ctFitModelSetup()`, which ctFit() and the generator run) converts
# before anything rewrites the model -- the augmented route moves a
# varying intercept into a carrier state that holds it raw, the scaling into
# the cell that reads the state -- and keeps the result as
# `fixedpopmeansraw`, used while it still describes the stated means.
#' @keywords internal
.ctFixedPopMeans <- function(model) {
  means <- model[['fixedpopmeans']]
  if (!length(means)) return(NULL)
  stored <- model[['fixedpopmeansraw']]
  if (is.data.frame(stored) && identical(stored$param, names(means)) &&
      identical(stored$value, as.numeric(means))) {
    return(stored)
  }
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
    row <- .ctParamRows(model, name)
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
  rows <- match(means$param, table$param)
  old <- as.integer(table$parnumber[rows])
  if (anyNA(old) || any(old < 1L) || any(old > nfull)) {
    stop("Internal layout error: a stated population mean has no raw position. ",
      "This is a bug in ctsem rather than in the model; please report it.",
      call. = FALSE)
  }
  free <- setdiff(seq_len(nfull), old)
  map <- integer(nfull)
  map[free] <- seq_along(free)
  map[old] <- length(free) + seq_along(old)
  # A stated value takes no prior: it is not a parameter.
  if (!is.null(priors)) {
    keep <- !as.integer(priors$index) %in% old
    priors$index <- priors$index[keep]
    priors$scale <- priors$scale[keep]
  }
  moved <- .ctJuliaRenumberRaw(list(parameter_table = parameter_table,
    laplace = laplace, ti_effects = ti_effects, ti_fixed = ti_fixed,
    ti_missing = ti_missing, random_effects = random_effects, priors = priors),
    map)
  c(moved, list(anchors = list(param = means$param, value = means$value,
    raw = means$raw, index = as.integer(length(free) + seq_along(old)))))
}

# Every raw position `x` names, renumbered by `map` (`map[k]` the new position
# of position k): a specification, or a list of its pieces, with the same
# names. The parameter table's numbers and the `param[k]` its rendered
# transforms read, the Laplace levels' indices, the TI effects, fixed TI
# effects and sampled TI predictor values, the augmented route's random
# effects, the priors and the smooth terms' penalties. A piece `x` lacks is
# left out. The one place that knows which fields hold raw positions: a field
# added to the specification that holds one is renumbered here.
#' @keywords internal
.ctJuliaRenumberRaw <- function(x, map) {
  remap <- function(v) {
    if (is.null(v)) return(v)
    v <- as.integer(v)
    use <- !is.na(v) & v > 0L
    v[use] <- map[v[use]]
    v
  }
  # The rendered transforms read raw positions by number, `param[k]`.
  retext <- function(v) {
    if (is.null(v)) return(v)
    out <- as.character(v)
    for (i in which(!is.na(out) & grepl("param[", out, fixed = TRUE))) {
      at <- gregexpr("\\bparam\\[([0-9]+)\\]", out[i], perl = TRUE)
      hits <- regmatches(out[i], at)[[1]]
      if (!length(hits)) next
      k <- as.integer(sub("^param\\[([0-9]+)\\]$", "\\1", hits))
      regmatches(out[i], at) <- list(paste0("param[", map[k], "]"))
    }
    out
  }
  table <- x[["parameter_table"]]
  if (!is.null(table)) {
    table$parnumber <- remap(table$parnumber)
    for (column in c("transform", "predicttransform", "updatetransform",
        "tdtransform")) {
      if (!is.null(table[[column]])) table[[column]] <- retext(table[[column]])
    }
    x$parameter_table <- table
  }
  laplace <- x[["laplace"]]
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
    x$laplace <- laplace
  }
  if (NROW(x[["ti_effects"]])) {
    x$ti_effects$parameter <- remap(x$ti_effects$parameter)
    x$ti_effects$coefficient <- remap(x$ti_effects$coefficient)
  }
  if (NROW(x[["ti_fixed"]])) x$ti_fixed$parameter <- remap(x$ti_fixed$parameter)
  if (NROW(x[["ti_missing"]])) x$ti_missing$parameter <- remap(x$ti_missing$parameter)
  if (NROW(x[["random_effects"]])) {
    x$random_effects$parameter <- remap(x$random_effects$parameter)
  }
  if (!is.null(x[["priors"]])) x$priors$index <- remap(x$priors$index)
  if (NROW(x[["penalty"]])) {
    x$penalty$parameter <- remap(x$penalty$parameter)
    x$penalty$linear <- remap(x$penalty$linear)
    x$penalty$index <- lapply(x$penalty$index, remap)
  }
  x
}

# `spec` with its stated population means moved on to follow `nfree` free
# parameters, every position naming one renumbered (`.ctJuliaRenumberRaw()`).
# Raw parameters added to a prepared specification are free ones, and the
# engine reads whatever follows the free parameters as the stated values
# (`CTSEMAnchors`), so room is made for them before the stated means: numbered
# from `.ctBackendNpar(spec) + 1` without it, sampled TI predictor values took
# the stated means' positions. A specification without stated means, or with
# them already there, comes back as it is.
#' @keywords internal
.ctJuliaAnchorsAfter <- function(spec, nfree) {
  index <- as.integer(spec$anchors$index)
  if (!length(index)) return(spec)
  nfree <- as.integer(nfree)[1L]
  now <- .ctBackendNpar(spec)
  if (nfree == now) return(spec)
  if (is.na(nfree) || nfree < now || !identical(index, now + seq_along(index))) {
    stop("Internal layout error: stated population means at ",
      paste(index, collapse = ", "), " cannot follow ", nfree, " free ",
      "parameters. This is a bug in ctsem; please report it.", call. = FALSE)
  }
  map <- seq_len(now + length(index))
  map[index] <- nfree + seq_along(index)
  spec <- .ctJuliaRenumberRaw(spec, map)
  spec$anchors$index <- as.integer(nfree + seq_along(index))
  spec
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
