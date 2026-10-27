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
# under every transform -- unless `values` states it, and every population
# sd is its `sdscale` on the raw scale unless RAWPOPVAR or `values` states
# it: the prior's scales used as fixed values. Free correlations are zero and
# free TI predictor effects zero unless `values` states them; a fixed effect
# (`'mm||||TI1=0.5'`) shifts the raw parameter by that much per unit of the
# predictor, as fitting applies it. `values` takes each value on the scale
# summary() reports it (`.ctGenerateValues`), so a model written for fitting
# generates from the values one would compare its fit with. All of it is named
# in one message.
#
# The design -- subjects, groups, times, missingness, predictor values -- is
# either built (`.ctGenerateSkeleton`) or the caller's (`.ctGenerateDesign`).
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

# A design the caller supplies (`ctGenerate(datastruct = )`), as the skeleton
# generation fills in: the subjects, groups, times and predictor values as
# given, ordered by subject and time as the engine reads them. A manifest cell
# that is NA stays missing; any other value, or a manifest column left out, is
# generated. A predictor the model names and the design lacks takes the
# balanced design's default -- zero for a TD predictor, a standard normal per
# subject for a TI predictor -- and is named in a message. `original` and
# `order` put the generated values back into the caller's own rows.
#' @keywords internal
.ctGenerateDesign <- function(model, datastruct, quiet = FALSE) {
  original <- as.data.frame(datastruct, stringsAsFactors = FALSE)
  idname <- model$subjectIDname
  timename <- model$timeName
  needed <- c(idname, timename, model$groupIDnames)
  absent <- setdiff(needed, names(original))
  if (length(absent)) stop("datastruct needs the column",
    if (length(absent) > 1L) "s", " ", paste(absent, collapse = ", "),
    ", which the model names for its subjects, times or groups.", call. = FALSE)
  if (!nrow(original)) stop("datastruct has no rows.", call. = FALSE)
  if (anyNA(original[[idname]]) || anyNA(original[[timename]])) {
    stop("datastruct has missing ids or times.", call. = FALSE)
  }
  # The order `.ctJuliaPrepare()` gives the data the engine reads, which is
  # the order generated values come back in: by id's own sort, then time.
  # First appearance put "b" before "a", and every subject's values back on
  # another subject's rows.
  order <- order(original[[idname]], original[[timename]])
  rows <- original[order, needed, drop = FALSE]
  rownames(rows) <- NULL
  observed <- matrix(TRUE, nrow(rows), length(model$manifestNames),
    dimnames = list(NULL, model$manifestNames))
  for (nm in model$manifestNames) {
    if (!is.null(original[[nm]])) observed[, nm] <- !is.na(original[[nm]][order])
    rows[[nm]] <- ifelse(observed[, nm], 0, NA_real_)
  }
  defaulted <- character()
  for (nm in model$TDpredNames) {
    if (is.null(original[[nm]])) {
      rows[[nm]] <- 0
      defaulted <- c(defaulted, paste0(nm, " at 0"))
    } else rows[[nm]] <- as.numeric(original[[nm]][order])
  }
  subject <- match(rows[[idname]], unique(rows[[idname]]))
  for (nm in model$TIpredNames) {
    if (is.null(original[[nm]])) {
      rows[[nm]] <- stats::rnorm(max(subject))[subject]
      defaulted <- c(defaulted, paste0(nm, " drawn standard normal"))
      next
    }
    value <- as.numeric(original[[nm]][order])
    if (anyNA(value)) stop("datastruct: TI predictor ", nm, " has missing ",
      "values, and generation needs each subject's.", call. = FALSE)
    if (any(tapply(value, subject, function(v) length(unique(v))) > 1L)) {
      stop("datastruct: TI predictor ", nm, " must be constant within each ",
        "subject.", call. = FALSE)
    }
    rows[[nm]] <- value
  }
  if (length(defaulted) && !quiet) message("Not in datastruct: ",
    paste(defaulted, collapse = ", "), ".")
  list(skeleton = rows, order = order, observed = observed, original = original)
}

# `values`: a named vector is the free parameters' population values, on
# their natural scale; a list holds summary()'s own sections -- `popmeans`
# (the same), `popsd` (or `popsd.<level>`) the random-effect sds,
# `rawpopcorr` (or `rawpopcorr.<level>`) their raw correlations, `tipreds` the
# linearised TI predictor effects -- each a named vector, named as that
# table's rows.
#' @keywords internal
.ctGenerateValues <- function(values) {
  if (!length(values)) return(list(means = NULL, sections = list()))
  named <- function(v) !is.null(names(v)) && all(nzchar(names(v)))
  if (!is.list(values)) {
    if (!is.numeric(values) || !named(values)) stop("values must be a named ",
      "numeric vector of population values, or a list of summary()'s ",
      "sections.", call. = FALSE)
    return(list(means = values, sections = list()))
  }
  section <- grepl("^(popmeans|tipreds)$|^(popsd|rawpopcorr)(\\..+)?$",
    names(values))
  if (!named(values) || !all(section)) stop("values, as a list, holds ",
    "summary()'s sections -- popmeans, popsd, rawpopcorr and tipreds, with ",
    ".<level> on popsd and rawpopcorr where there are levels",
    if (named(values)) paste0(" -- not ", paste(names(values)[!section],
      collapse = ", ")), ". A named vector gives population values alone.",
    call. = FALSE)
  for (nm in names(values)) {
    v <- values[[nm]]
    if (!is.numeric(v) || !named(v)) stop("values$", nm, " must be a named ",
      "numeric vector, named as summary() names that table's rows.",
      call. = FALSE)
  }
  list(means = values$popmeans,
    sections = values[setdiff(names(values), "popmeans")])
}

# The raw sd at which a random effect spreads its parameter by `target` on the
# natural scale, as summary() measures it: the Gauss-Hermite spread of the
# transform about the raw mean (`.ctBackendRandomEffectLevel`). The spread is
# zero at zero and grows with the raw sd for a monotone transform, so a
# bracket holds the one root.
#' @keywords internal
.ctGenerateRawSd <- function(transform, mean, target, name) {
  if (!is.finite(target) || target <= 0) stop("values$popsd: ", name,
    " must be a positive number; a parameter with no spread does not vary.",
    call. = FALSE)
  rule <- .ctBackendGaussHermite()
  spread <- function(s) {
    v <- .ctGenerateNatural(transform, mean + s * rule$node)
    sqrt(sum(rule$weight * (v - sum(rule$weight * v))^2))
  }
  upper <- 1
  while (!(spread(upper) >= target)) {
    upper <- upper * 2
    if (upper > 2^20) stop("values$popsd: ", name, " = ", format(target),
      " is not a spread its transform (", transform, ") can produce about ",
      "its mean.", call. = FALSE)
  }
  stats::uniroot(function(s) spread(s) - target, c(0, upper), tol = 1e-12)$root
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
  stop("values: ", name, " = ", format(value), " is not a value its ",
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
.ctGeneratePrepare <- function(model, skeleton, values = NULL,
  project = NULL, quiet = FALSE) {
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

  # The model set up as a fit sets it up (`.ctFitModelSetup()`, the same code
  # ctFit() runs), then prepared as a fit prepares it. Messages from both are
  # about fitting -- a T0VAR row the population block covers, a level with
  # fewer groups than its rank -- and generation names its own values below.
  # The setup draws from the random stream (`ctModelTransformsToNum()`'s grid),
  # which a fit keeps for its starting values; the generator's draws follow
  # from the seed alone, so the stream is put back.
  # The route is the one an optimising julia fit takes, chosen on the model
  # as set up: a 'stationary' T0MEANS cell is a varying intercept until the
  # setup fixes it.
  keep <- .ctStreamKeep()
  setup <- tryCatch(suppressMessages(.ctFitModelSetup(model, intoverpop = "auto",
    priors = FALSE, nsubjects = length(unique(skeleton[[model$subjectIDname]])),
    checks = FALSE)), finally = keep())
  route <- if (identical(setup$intoverpopmethod, "laplace")) "laplace" else "augmented"
  spec <- suppressMessages(.ctJuliaPrepare(skeleton, setup$ctm, project = project,
    priors = FALSE, intoverpop = route))
  # The transforms the reads below take, as the fit's are: covariances given
  # as ctCov() in the construction in force now.
  model <- .ctCovRefresh(model)
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  raw <- numeric(.ctBackendNpar(spec))

  table <- spec$parameter_table
  free <- !is.na(table$parnumber)
  parnames <- unique(as.character(table$param[free]))
  transformOf <- function(name) {
    model$pars$transform[which(model$pars$param %in% name)[1L]]
  }

  # Population means: the prior's centre, raw zero, unless stated.
  given <- .ctGenerateValues(values)
  means <- given$means
  sectionsused <- character()
  if (length(means)) {
    unknown <- setdiff(names(means), parnames)
    if (length(unknown)) {
      stop("values names ", paste(unknown, collapse = ", "), ", which ",
        if (length(unknown) > 1L) "are not free parameters" else
          "is not a free parameter", " of this model.", call. = FALSE)
    }
    stated <- intersect(names(means), spec$anchors$param)
    if (length(stated)) {
      stop("values gives a population mean for ", paste(stated, collapse = ", "),
        ", which the model states (fixedpopmeans); change it there.", call. = FALSE)
    }
    for (nm in names(means)) {
      index <- unique(table$parnumber[free & table$param %in% nm])
      raw[index] <- .ctGenerateRawFor(transformOf(nm), means[[nm]], nm)
    }
  }
  centred <- setdiff(parnames, c(names(means), spec$anchors$param))
  natural <- vapply(centred, function(nm)
    .ctGenerateNatural(transformOf(nm), 0), numeric(1))

  # Population spread, per level. What RAWPOPVAR (RAWPOPVAR_<level>) states is
  # already fixed in the specification, as fitting holds it; a free sd is the
  # level's sdscale, a free correlation coordinate zero. All raw-scale.
  spread <- character()
  # A level's section of `values`: `<name>.<level>`, or the bare name for
  # the subject level.
  sectionFor <- function(name, level) {
    keys <- c(paste0(name, ".", level$name),
      if (identical(level$name, model$subjectIDname)) name)
    hit <- intersect(keys, names(given$sections))
    if (length(hit) > 1L) stop("values gives both ", paste(hit, collapse = " and "),
      ", which name the same level.", call. = FALSE)
    if (!length(hit)) return(NULL)
    sectionsused <<- c(sectionsused, hit)
    given$sections[[hit]]
  }
  for (level in spec$laplace$levels) {
    sds <- sectionFor("popsd", level)
    cors <- sectionFor("rawpopcorr", level)
    if ((length(sds) || length(cors)) && length(level$load_index)) {
      stop("values: level ", level$name, " has a reduced-rank population ",
        "covariance (poprank), whose spread popsd and rawpopcorr cannot set.",
        call. = FALSE)
    }
    unknown <- setdiff(names(sds), level$param)
    if (length(unknown)) stop("values$popsd names ", paste(unknown,
      collapse = ", "), ", which ", if (length(unknown) > 1L) "do" else "does",
      " not vary at level ", level$name, ".", call. = FALSE)
    if (length(level$sd_index)) {
      for (j in seq_along(level$sd_index)) {
        free <- level$sd_index[j] > 0L
        name <- level$param[j]
        if (name %in% names(sds)) {
          if (!free) stop("values$popsd: ", name, "'s spread is stated in ",
            "the model's RAWPOPVAR.", call. = FALSE)
          # The mean may be stated in the model (R/ctFixedPopMeans.R), past
          # the free parameters.
          target <- .ctGenerateRawSd(transformOf(name),
            .ctBackendRawFull(spec, raw)[level$re_index[j]], sds[[name]], name)
        } else {
          target <- if (free) level$sd_scale[j] else
            log1p(exp(2 * level$sd_fixed[j] - 1)) * level$sd_scale[j]
        }
        if (free) raw[level$sd_index[j]] <- .ctJuliaRawPopSd(target, level$sd_scale[j])
        spread <- c(spread, sprintf("%s %s [%s]", name,
          format(signif(target, 4)), level$name))
      }
      if (length(cors)) {
        k <- length(level$param)
        R <- diag(k)
        for (nm in names(cors)) {
          pair <- match(strsplit(nm, "__", fixed = TRUE)[[1L]], level$param)
          if (length(pair) != 2L || anyNA(pair) || pair[1L] == pair[2L]) {
            stop("values$rawpopcorr: ", nm, " is not a pair of parameters ",
              "varying at level ", level$name, " (named a__b, as summary() ",
              "names them).", call. = FALSE)
          }
          R[pair[1L], pair[2L]] <- R[pair[2L], pair[1L]] <- cors[[nm]]
        }
        index <- as.integer(level$cor_index)
        if (any(index == 0L)) stop("values$rawpopcorr: level ", level$name,
          " states some of its correlations in RAWPOPVAR.", call. = FALSE)
        construction <- c("0" = "rawcorr", "2" = "z")[as.character(
          if (is.null(spec$covmatcode)) 0L else spec$covmatcode)]
        if (is.na(construction)) stop("values$rawpopcorr needs ",
          "covmattransform 'z' or 'rawcorr'.", call. = FALSE)
        cells <- .ctCovCells(R, construction,
          what = paste0("values$rawpopcorr at level ", level$name))
        lower <- which(lower.tri(R), arr.ind = TRUE)
        lower <- lower[order(lower[, 2L], lower[, 1L]), , drop = FALSE]
        raw[index] <- cells[lower]
      }
    } else if (length(level$load_index)) {
      # A reduced-rank level: each basis effect's own loading at raw 1, which
      # is one sdscale of spread, as a fit starts from.
      raw[.ctJuliaLoadingStart(spec)] <- 1
      spread <- c(spread, sprintf("rank %d loadings at their start [%s]",
        level$rank, level$name))
    }
    # Own variances (`popunique`) at one sdscale of spread, as a free
    # full-rank sd is; a rank-0 level has these and no loadings.
    own <- as.integer(level$unique_index)
    for (j in which(own > 0L)) {
      raw[own[j]] <- .ctJuliaRawPopSd(level$sd_scale[j], level$sd_scale[j])
      spread <- c(spread, sprintf("%s own %s [%s]", level$param[j],
        format(signif(level$sd_scale[j], 4)), level$name))
    }
  }

  # TI effects: a fixed value where the model states one, the raw coefficient
  # for a linearised effect `values$tipreds` states, zero otherwise. The
  # summary reports an effect as the raw coefficient times the transform's
  # slope at the raw mean (`.ctBackendTipredDraws`), so this divides by it.
  zeroed <- character()
  tipreds <- given$sections$tipreds
  if (!is.null(tipreds)) sectionsused <- c(sectionsused, "tipreds")
  stated <- character()
  ti <- spec$ti_effects
  for (r in seq_len(NROW(ti))) {
    name <- as.character(table$param[match(ti$parameter[r], table$parnumber)])
    value <- fixedti[[paste(name, ti$predictor[r])]]
    key <- paste0("tip_", model$TIpredNames[ti$predictor[r]], "_", name)
    if (key %in% names(tipreds)) {
      if (!is.null(value)) stop("values$tipreds: ", key, " is fixed in the ",
        "model.", call. = FALSE)
      # The target may be a stated population mean (R/ctFixedPopMeans.R).
      mean <- .ctBackendRawFull(spec, raw)[ti$parameter[r]]
      h <- 1e-5
      slope <- diff(.ctGenerateNatural(transformOf(name), mean + c(-h, h))) / (2 * h)
      if (!is.finite(slope) || abs(slope) < 1e-12) stop("values$tipreds: ",
        key, " cannot be set, because ", name, "'s transform is flat at its ",
        "mean.", call. = FALSE)
      raw[ti$coefficient[r]] <- tipreds[[key]] / slope
      stated <- c(stated, key)
    } else if (is.null(value)) {
      zeroed <- c(zeroed, paste0(model$TIpredNames[ti$predictor[r]], " on ", name))
    } else raw[ti$coefficient[r]] <- value
  }
  unknown <- setdiff(names(tipreds), stated)
  if (length(unknown)) stop("values$tipreds names ", paste(unknown,
    collapse = ", "), ", which ", if (length(unknown) > 1L) "are not" else
    "is not", " a free TI predictor effect of this model (named tip_<predictor>_",
    "<parameter>, as summary() names them).", call. = FALSE)
  unread <- setdiff(names(given$sections), sectionsused)
  if (length(unread)) stop("values: ", paste(unread, collapse = ", "),
    " name", if (length(unread) == 1L) "s", " no level of this model.",
    call. = FALSE)

  if (!quiet) {
    listed <- function(x) paste0(paste(utils::head(x, 8), collapse = ", "),
      if (length(x) > 8L) ", ..." else "")
    parts <- c(
      if (length(centred)) paste0("at prior centres (raw 0): ",
        listed(paste0(centred, "=", vapply(natural,
          function(x) format(signif(x, 4)), "")))),
      if (length(values)) paste0("from values: ",
        listed(names(values))),
      if (length(spread)) paste0("population sd (raw): ", listed(spread)),
      if (length(zeroed)) paste0("TI effects at zero: ", listed(zeroed)))
    if (length(parts)) message("Generating values -- ",
      paste(parts, collapse = "; "), ".")
  }
  list(handle = handle, raw = raw)
}

# How many standard normals one draw of every random effect takes; zero when the
# specification has none.
.ctGenerateEffectDimension <- function(handle) {
  spec <- .ctBackendSpec(handle)
  if (is.null(spec$laplace)) return(0L)
  module <- .ctJuliaModule(spec$project)
  as.integer(.ctBackendJuliaValue(
    module$ctsem_laplace_coordinate_dimension(.ctJuliaObjective(handle))))
}

# One draw of every random effect from the population, as natural deviations
# for `effects`. NULL when the specification has none.
#' @keywords internal
.ctGeneratePopulationEffects <- function(handle, raw) {
  n <- .ctGenerateEffectDimension(handle)
  if (n < 1L) return(NULL)
  spec <- .ctBackendSpec(handle)
  module <- .ctJuliaModule(spec$project)
  as.numeric(.ctBackendJuliaValue(module$ctsem_laplace_population_deviations(
    .ctJuliaObjective(handle), .ctJuliaNumericVector(as.numeric(raw)),
    .ctJuliaNumericVector(stats::rnorm(n)))))
}

#' @keywords internal
.ctGenerateJulia <- function(model, n.subjects, times, values = NULL,
  project = NULL, quiet = FALSE, intoverstates = FALSE, design = NULL) {

  if (!.ctJuliaConnectoROk()) stop(.ctJuliaConnectoRNeeded(), call. = FALSE)
  designed <- if (is.null(design)) NULL else
    .ctGenerateDesign(model, design, quiet = quiet)
  skeleton <- if (is.null(designed)) .ctGenerateSkeleton(model, n.subjects, times) else
    designed$skeleton
  prepared <- .ctGeneratePrepare(model, skeleton, values = values,
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
  if (!is.null(designed)) {
    for (nm in model$manifestNames) skeleton[[nm]][!designed$observed[, nm]] <- NA
  }
  # A matrix, because that is what `ctGenerate`'s own path returns and callers
  # index it positionally. What generated it rides along as attributes, keyed
  # by id (and time) so they survive the burnin trim: each subject's
  # parameters, and on the state-explicit route the latent states.
  out <- if (is.null(designed)) as.matrix(skeleton) else {
    # The caller's own data frame and rows, with the generated values (and any
    # predictor it lacked) filled in.
    filled <- designed$original
    for (nm in c(model$manifestNames, model$TDpredNames, model$TIpredNames)) {
      if (is.null(filled[[nm]])) filled[[nm]] <- NA_real_
      filled[designed$order, nm] <- skeleton[[nm]]
    }
    filled
  }
  attr(out, "subjectpars") <- .ctGenerateSubjectPars(handle, raw, effects,
    model, skeleton)
  grouppars <- .ctGenerateGroupPars(handle, raw, effects, model, skeleton)
  if (length(grouppars)) attr(out, "grouppars") <- grouppars
  if (!isTRUE(intoverstates) && is.list(generated) && !is.null(generated$states)) {
    states <- t(as.matrix(generated$states))
    nl <- length(model$latentNames)
    colnames(states) <- c(model$latentNames,
      if (ncol(states) > nl) paste0("state", seq_len(ncol(states) - nl)))
    attr(out, "states") <- cbind(skeleton[, c(model$subjectIDname,
      model$timeName)], states[, seq_len(min(nl, ncol(states))), drop = FALSE])
  }
  out
}

# Each grouping level's own parameters, as generated: the population values
# shifted by that level's random effects and those of the levels above it,
# with every level inside it at zero. One data frame per grouping level, one
# row per group. TI predictor effects are a subject's and are left out.
#' @keywords internal
.ctGenerateGroupPars <- function(handle, raw, effects, model, skeleton) {
  spec <- .ctBackendSpec(handle)
  if (is.null(spec$laplace) || !length(effects) || !length(model$groupIDnames)) {
    return(NULL)
  }
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(handle)
  position <- as.integer(.ctBackendJuliaValue(
    module$ctsem_laplace_deviation_layout(objective))$level)
  out <- list()
  for (l in seq_along(model$groupIDnames)) {
    kept <- effects
    kept[position <= l] <- 0
    pars <- .ctGenerateSubjectPars(handle, raw, kept, model, skeleton,
      ti_effects = FALSE)
    group <- model$groupIDnames[l]
    pars <- pars[!duplicated(pars[[group]]), , drop = FALSE]
    pars <- pars[, setdiff(names(pars), c(model$subjectIDname,
      model$groupIDnames[seq_len(l - 1L)])), drop = FALSE]
    rownames(pars) <- NULL
    out[[group]] <- pars
  }
  out
}

# Each subject's parameters as generated, on their natural scale: the
# population values shifted by that subject's drawn random effects at every
# level and by its TI predictor effects. One row per subject, one column per
# free parameter, with the subject's groups beside its id.
#' @keywords internal
.ctGenerateSubjectPars <- function(handle, raw, effects, model, skeleton,
  ti_effects = TRUE) {
  spec <- .ctBackendSpec(handle)
  # Read by position below, a stated population mean's among them
  # (R/ctFixedPopMeans.R).
  raw <- .ctBackendRawFull(spec, as.numeric(raw))
  table <- spec$parameter_table
  free <- !is.na(table$parnumber)
  if (!any(free)) return(NULL)
  parnames <- unique(as.character(table$param[free]))
  parnumber <- table$parnumber[free][match(parnames, as.character(table$param[free]))]
  first <- !duplicated(skeleton[[model$subjectIDname]])
  nsubjects <- sum(first)
  values <- if (!is.null(spec$laplace) && length(effects)) {
    module <- .ctJuliaModule(spec$project)
    matrix(as.numeric(.ctBackendJuliaValue(module$ctsem_laplace_subject_values(
      .ctJuliaObjective(handle), .ctJuliaNumericVector(as.numeric(raw)),
      .ctJuliaNumericVector(as.numeric(effects)), ti_effects = ti_effects))),
      nrow = nsubjects)
  } else {
    # No random effects: the population vector, shifted by TI effects alone.
    v <- matrix(as.numeric(raw), nsubjects, length(raw), byrow = TRUE)
    ti <- spec$ti_effects
    for (r in seq_len(NROW(ti))) {
      predictor <- skeleton[[model$TIpredNames[ti$predictor[r]]]][first]
      v[, ti$parameter[r]] <- v[, ti$parameter[r]] + raw[ti$coefficient[r]] * predictor
    }
    v
  }
  transformOf <- function(name) {
    model$pars$transform[which(model$pars$param %in% name)[1L]]
  }
  natural <- vapply(seq_along(parnames), function(k)
    .ctGenerateNatural(transformOf(parnames[k]), values[, parnumber[k]]),
    numeric(nsubjects))
  natural <- matrix(natural, nrow = nsubjects, dimnames = list(NULL, parnames))
  keys <- skeleton[first, c(model$subjectIDname, model$groupIDnames), drop = FALSE]
  rownames(keys) <- NULL
  cbind(keys, as.data.frame(natural))
}
