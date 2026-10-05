# Generating data from a model through the julia engine.
#
# This is generation from the *model*, not from any way of fitting it: each
# subject's parameters are drawn from the population the model describes --
# at every level of its hierarchy -- and its latent trajectory is simulated
# from the process at those parameters, each observation then drawn from its
# measurement model given the state at its row. Nothing is integrated and
# nothing is linearised, so a nonlinear drift, a categorical or count
# indicator and a random effect on the dynamics are all generated exactly.
# Which fitting route would suit the data is a separate question, and
# `ctGenerateFromFit()` is where the routes come back in.
#
# The engine pieces are the ones fitting and sampling already use: the
# Laplace objective (`intoverpop = 'laplace'`) describes the random effects as
# raw-parameter indices, `ctsem_laplace_population_deviations` turns standard
# normals into each block's deviation at the population covariance, and
# `ctsem_generate_states` simulates each subject at its shifted parameters.
#
# Every value the model leaves free is set at its prior's centre -- raw zero,
# under every transform -- unless `popmeans` states it, and every population
# sd is its `sdscale` on the raw scale unless RAWPOPVAR states it: the prior's
# scales used as fixed values. Free TI predictor effects are zero; a fixed one
# (`'mm||||TI1=0.5'`) shifts the raw parameter by that much per unit of the
# predictor, as fitting applies it. All of it is named in one message.
#
# The standard normals are drawn on the R side, in a fixed order -- the TI
# predictors, the random effects, the latent innovations, the observation
# deviates -- so `set.seed()` governs the whole dataset.

# The values generation supplies when a matrix is free, for the R generator
# (`ctModeltoNumeric()`), which integrates a linear system at fixed values and
# has no prior to centre on. Three diagonals cannot take a literal zero because
# zero is singular there: DRIFT (`fQinf()` cannot solve for the asymptotic
# covariance), DIFFUSION and T0VAR (a covariance that cannot be factorised).
#' @keywords internal
.ctGenerateDefaults <- function(continuoustime = TRUE) {
  driftdiagonal <- if(isTRUE(continuoustime)) -1e-6 else 1 - 1e-6
  list(
    DRIFT = list(diagonal = driftdiagonal, offdiagonal = 0),
    DIFFUSION = list(diagonal = 1e-6, offdiagonal = 0),
    T0VAR = list(diagonal = 1e-6, offdiagonal = 0),
    MANIFESTVAR = list(diagonal = 0, offdiagonal = 0),
    LAMBDA = list(diagonal = 0, offdiagonal = 0),
    T0MEANS = list(diagonal = 0, offdiagonal = 0),
    MANIFESTMEANS = list(diagonal = 0, offdiagonal = 0),
    CINT = list(diagonal = 0, offdiagonal = 0),
    TDPREDEFFECT = list(diagonal = 0, offdiagonal = 0),
    PARS = list(diagonal = 0, offdiagonal = 0))
}

# How many units at each level of the hierarchy, innermost first: subjects,
# then each grouping level in `c(model$subjectIDname, model$groupIDnames)`
# order. A level not given has one group.
#' @keywords internal
.ctGenerateLevelCounts <- function(model, n.subjects) {
  levels <- c(model$subjectIDname, model$groupIDnames)
  n <- suppressWarnings(as.integer(n.subjects))
  if (!length(n) || anyNA(n) || any(n < 1L) || any(n != n.subjects)) {
    stop("n must be positive whole numbers.", call. = FALSE)
  }
  if (length(n) > length(levels)) {
    stop("n has ", length(n), " counts for the model's ",
      length(levels), " id level", if (length(levels) > 1L) "s" else "", " (",
      paste(levels, collapse = ", "), ").", call. = FALSE)
  }
  n <- c(n, rep(1L, length(levels) - length(n)))
  if (is.unsorted(rev(n))) {
    stop("n must not grow from one level to the next: each ",
      "grouping level holds the units of the level inside it (",
      paste(levels, n, sep = " = ", collapse = ", "), ").", call. = FALSE)
  }
  stats::setNames(n, levels)
}

# A long-format skeleton with the requested units and times and no
# observations: the shape generation fills in.
#
# The id, time and group columns are named as the model names them:
# `.ctJuliaPrepare()` looks them up by name. Groups are nested and as even as
# the counts allow, each level's units split in order over the level above.
#' @keywords internal
.ctGenerateSkeleton <- function(model, n.subjects, times) {
  counts <- .ctGenerateLevelCounts(model, n.subjects)
  nsubjects <- counts[[1L]]
  idname <- model$subjectIDname
  timename <- model$timeName
  rows <- do.call(rbind, lapply(seq_len(nsubjects), function(i)
    stats::setNames(data.frame(i, times[[i]]), c(idname, timename))))
  unit <- seq_len(nsubjects)
  for (l in seq_along(model$groupIDnames)) {
    unit <- ceiling(unit * counts[[l + 1L]] / counts[[l]])
    rows[[model$groupIDnames[l]]] <- unit[rows[[idname]]]
  }
  # Zero, not NA. The engine generates only where an observation exists, so
  # the skeleton's missingness pattern is the input and its values are not.
  for (nm in model$manifestNames) rows[[nm]] <- 0
  # Time-dependent predictors at zero: the specification carries no
  # distribution for them.
  for (nm in model$TDpredNames) rows[[nm]] <- 0
  # Time-independent predictors are drawn, one standard normal per subject,
  # so an effect has variation to act on and reads as change per standard
  # deviation.
  for (nm in model$TIpredNames) {
    rows[[nm]] <- stats::rnorm(nsubjects)[rows[[idname]]]
  }
  rows
}

# A parameter's natural value at raw `raw`, through the transform the model
# states for it: an expression in `param`, or one of ctsem's numbered shapes.
#' @keywords internal
.ctGenerateNatural <- function(transform, raw) {
  code <- suppressWarnings(as.integer(transform))
  if (!is.na(code)) {
    return(vapply(raw, function(r) as.numeric(tform(r, code, 1, 1, 0, 0)),
      numeric(1)))
  }
  text <- parse(text = as.character(transform))
  vapply(raw, function(r) as.numeric(eval(text, list(param = r),
    environment(.ctGenerateNatural))), numeric(1))
}

# The raw value whose transform gives `value`. Transforms are monotone, so a
# bracket that changes sign holds the one root.
#' @keywords internal
.ctGenerateRawFor <- function(transform, value, name) {
  f <- function(r) .ctGenerateNatural(transform, r) - value
  for (width in 2^(0:10)) {
    ends <- f(c(-width, width))
    if (anyNA(ends)) break
    if (any(ends == 0)) return(c(-width, width)[ends == 0][1L])
    if (sign(ends[1L]) != sign(ends[2L])) {
      return(stats::uniroot(f, c(-width, width), tol = 1e-12)$root)
    }
  }
  stop("popmeans: ", name, " = ", format(value), " is not a value its ",
    "transform (", transform, ") can take.", call. = FALSE)
}

# The generation model, its prepared specification, and the raw vector that
# fixes every value generation needs, with a message naming them.
#
# TI effects are prepared as free coefficients whatever the model says --
# a fixed one is relabelled so it gets a coefficient of its own -- and the
# value is written into the raw vector. The specification then carries every
# effect the same way, and `.ctJuliaTIEffects()` is reused as it stands.
#' @keywords internal
.ctGeneratePrepare <- function(model, skeleton, popmeans = NULL,
  project = NULL, quiet = FALSE) {
  model <- .ctCovRefresh(model)
  fixedti <- list()
  for (k in seq_along(model$TIpredNames)) {
    column <- paste0(model$TIpredNames[k], "_effect")
    if (is.null(model$pars[[column]])) next
    value <- .ctTipredEffectValue(model$pars[[column]])
    fixed <- which(!is.na(value) & !is.na(model$pars$param))
    for (i in fixed) {
      fixedti[[paste(model$pars$param[i], k)]] <- value[i]
      model$pars[[column]][i] <- paste0("ctgen_", model$TIpredNames[k], "_",
        make.names(model$pars$param[i]))
    }
  }

  varying <- .ctAnyVarying(model)
  model <- .ctModelRawPopVarSync(model)
  # Messages from preparation are about fitting -- a T0VAR row the population
  # block covers, a level with fewer groups than its rank -- and generation
  # names its own values below.
  spec <- suppressMessages(.ctJuliaPrepare(skeleton, model, project = project,
    priors = FALSE, intoverpop = if (varying) "laplace" else "augmented"))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  raw <- numeric(.ctBackendNpar(spec))

  table <- spec$parameter_table
  free <- !is.na(table$parnumber)
  parnames <- unique(as.character(table$param[free]))
  transformOf <- function(name) {
    model$pars$transform[which(model$pars$param %in% name)[1L]]
  }

  # Population means: the prior's centre, raw zero, unless stated.
  if (length(popmeans)) {
    if (is.null(names(popmeans)) || any(!nzchar(names(popmeans)))) {
      stop("popmeans must be named by parameter.", call. = FALSE)
    }
    unknown <- setdiff(names(popmeans), parnames)
    if (length(unknown)) {
      stop("popmeans names ", paste(unknown, collapse = ", "), ", which ",
        if (length(unknown) > 1L) "are not free parameters" else
          "is not a free parameter", " of this model.", call. = FALSE)
    }
    for (nm in names(popmeans)) {
      index <- unique(table$parnumber[free & table$param %in% nm])
      raw[index] <- .ctGenerateRawFor(transformOf(nm), popmeans[[nm]], nm)
    }
  }
  centred <- setdiff(parnames, names(popmeans))
  natural <- vapply(centred, function(nm)
    .ctGenerateNatural(transformOf(nm), 0), numeric(1))

  # Population spread, per level. What RAWPOPVAR (RAWPOPVAR_<level>) states is
  # already fixed in the specification, as fitting holds it; a free sd is the
  # level's sdscale, a free correlation coordinate zero. All raw-scale.
  spread <- character()
  for (level in spec$laplace$levels) {
    if (length(level$sd_index)) {
      for (j in seq_along(level$sd_index)) {
        free <- level$sd_index[j] > 0L
        target <- if (free) level$sd_scale[j] else
          log1p(exp(2 * level$sd_fixed[j] - 1)) * level$sd_scale[j]
        if (free) raw[level$sd_index[j]] <- .ctJuliaRawPopSd(target, level$sd_scale[j])
        spread <- c(spread, sprintf("%s %s [%s]", level$param[j],
          format(signif(target, 4)), level$name))
      }
    } else if (length(level$load_index)) {
      # A reduced-rank level: each basis effect's own loading at raw 1, which
      # is one sdscale of spread, as a fit starts from.
      raw[.ctJuliaLoadingStart(spec)] <- 1
      spread <- c(spread, sprintf("rank %d loadings at their start [%s]",
        level$rank, level$name))
    }
  }

  # TI effects: a fixed value where the model states one, zero otherwise.
  zeroed <- character()
  ti <- spec$ti_effects
  for (r in seq_len(NROW(ti))) {
    name <- as.character(table$param[match(ti$parameter[r], table$parnumber)])
    value <- fixedti[[paste(name, ti$predictor[r])]]
    if (is.null(value)) {
      zeroed <- c(zeroed, paste0(model$TIpredNames[ti$predictor[r]], " on ", name))
    } else raw[ti$coefficient[r]] <- value
  }

  if (!quiet) {
    listed <- function(x) paste0(paste(utils::head(x, 8), collapse = ", "),
      if (length(x) > 8L) ", ..." else "")
    parts <- c(
      if (length(centred)) paste0("at prior centres (raw 0): ",
        listed(paste0(centred, "=", vapply(natural,
          function(x) format(signif(x, 4)), "")))),
      if (length(popmeans)) paste0("from popmeans: ",
        listed(names(popmeans))),
      if (length(spread)) paste0("population sd (raw): ", listed(spread)),
      if (length(zeroed)) paste0("TI effects at zero: ", listed(zeroed)))
    if (length(parts)) message("Generating values -- ",
      paste(parts, collapse = "; "), ".")
  }
  list(handle = handle, raw = raw, varying = varying)
}

# One draw of every random effect from the population, as natural deviations
# for `effects`. NULL when the specification has none.
#' @keywords internal
.ctGeneratePopulationEffects <- function(handle, raw) {
  spec <- .ctBackendSpec(handle)
  if (is.null(spec$laplace)) return(NULL)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(handle)
  n <- as.integer(.ctBackendJuliaValue(
    module$ctsem_laplace_coordinate_dimension(objective)))
  if (n < 1L) return(NULL)
  as.numeric(.ctBackendJuliaValue(module$ctsem_laplace_population_deviations(
    objective, .ctJuliaNumericVector(as.numeric(raw)),
    .ctJuliaNumericVector(stats::rnorm(n)))))
}

#' @keywords internal
.ctGenerateJulia <- function(model, n.subjects, times, popmeans = NULL,
  project = NULL, quiet = FALSE, intoverstates = FALSE) {

  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    stop("Generating with backend='julia' needs the JuliaConnectoR package.",
      call. = FALSE)
  }
  skeleton <- .ctGenerateSkeleton(model, n.subjects, times)
  prepared <- .ctGeneratePrepare(model, skeleton, popmeans = popmeans,
    project = project, quiet = quiet)
  handle <- prepared$handle
  # A zero-length vector deadlocks the JuliaConnectoR bridge, and a fully
  # fixed model has no free parameters, so one unread element is sent instead.
  raw <- if (length(prepared$raw)) prepared$raw else 0
  effects <- .ctGeneratePopulationEffects(handle, raw)

  nmanifest <- length(model$manifestNames)
  if (isTRUE(intoverstates)) {
    # The filter's one-step-ahead predictive, at the same subjects: the
    # diagnostic comparison, not the model.
    base <- matrix(stats::rnorm(nmanifest * nrow(skeleton)), nmanifest,
      nrow(skeleton))
    generated <- .ctBackendGenerate(handle, raw, base, effects = effects)
  } else {
    # The trajectory from the process, then each observation given the state
    # at its row. The engine says how many innovations the design needs.
    nz <- .ctBackendStateDimension(handle)
    z <- stats::rnorm(nz)
    base <- matrix(stats::rnorm(nmanifest * nrow(skeleton)), nmanifest,
      nrow(skeleton))
    generated <- .ctBackendGenerateStates(handle, raw, z, base,
      effects = effects)
  }
  drawn <- as.numeric(if (is.list(generated)) generated$Y else generated)
  if (length(drawn) != nmanifest * nrow(skeleton)) {
    stop("The engine returned ", length(drawn), " generated values where ",
      nmanifest * nrow(skeleton), " were expected.", call. = FALSE)
  }
  skeleton[, model$manifestNames] <- t(matrix(drawn, nrow = nmanifest))
  # A matrix, because that is what `ctGenerate`'s own path returns and callers
  # index it positionally.
  as.matrix(skeleton)
}
