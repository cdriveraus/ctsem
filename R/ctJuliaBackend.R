# Julia likelihood backend -------------------------------------------------
#
# The Julia process is deliberately process-local.  Fit objects retain only a
# serializable specification and rebuild their proxy on demand after saveRDS.

.ct_julia_cache <- new.env(parent = emptyenv())
.ct_julia_cache$objectives <- new.env(parent = emptyenv())
.ct_julia_engine_file <- function() system.file("julia", "engine.json", package = "ctsem")
.ctJuliaOr <- function(x, default) if (is.null(x)) default else x
.ctJuliaString <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || grepl('"', value, fixed = TRUE) ||
      grepl("\r", value, fixed = TRUE) || grepl("\n", value, fixed = TRUE)) {
    stop("Julia string literal must be one non-empty-line value.", call. = FALSE)
  }
  paste0('"', value, '"')
}

.ctJuliaBin <- function(julia_bin = NULL) {
  if (!is.null(julia_bin)) {
    return(normalizePath(julia_bin, winslash = "/", mustWork = TRUE))
  }
  configured <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (nzchar(configured) && file.exists(file.path(configured, "julia.exe"))) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  if (.Platform$OS.type == "windows") {
    profiles <- unique(c(Sys.getenv("USERPROFILE", unset = ""), path.expand("~")))
    for (profile in profiles[nzchar(profiles)]) {
      juliaup <- file.path(profile, ".julia", "juliaup")
      installs <- list.dirs(juliaup, recursive = FALSE, full.names = TRUE)
      bins <- file.path(installs, "bin")
      bins <- bins[file.exists(file.path(bins, "julia.exe"))]
      if (length(bins)) return(normalizePath(bins[[length(bins)]], winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

.ctJuliaRequire <- function() {
  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    stop("backend='julia' requires the suggested package JuliaConnectoR. Install it, install Julia, then call ctJuliaSetup().", call. = FALSE)
  }
  invisible(TRUE)
}

.ctJuliaEngineLock <- function() {
  path <- .ct_julia_engine_file()
  if (!nzchar(path) || !file.exists(path)) {
    stop("ctsem Julia engine lock file is unavailable; reinstall ctsem.", call. = FALSE)
  }
  json <- paste(readLines(path, warn = FALSE), collapse = "\n")
  fields <- c("url", "revision", "subdir", "julia")
  out <- lapply(fields, function(field) {
    hit <- regmatches(json, regexec(paste0('"', field, '"\\s*:\\s*"([^"]+)"'), json))[[1]]
    if (length(hit) < 2L) stop("Malformed ctsem Julia engine lock file.", call. = FALSE)
    hit[[2]]
  })
  names(out) <- fields
  out
}

#' Configure the Julia engine used by ctsem
#'
#' Installs or instantiates the version of ContinuousTimeSEM.jl pinned by this
#' ctsem source tree. This is explicit: ctFit never downloads Julia packages.
#' @param project Optional local ContinuousTimeSEM.jl project for development.
#' @param revision Engine revision, normally \code{"locked"}.
#' @param julia_bin Optional Julia binary directory.
#' @param force Reconfigure an existing Julia session.
#' @return A Julia-engine status list, invisibly.
#' @export
ctJuliaSetup <- function(project = NULL, revision = "locked", julia_bin = NULL,
  force = FALSE) {
  .ctJuliaRequire()
  if (isTRUE(force)) .ctJuliaClearSession()
  julia_bin <- .ctJuliaBin(julia_bin)
  if (!is.null(julia_bin)) Sys.setenv(JULIA_BINDIR = julia_bin)
  lock <- .ctJuliaEngineLock()
  if (!is.null(project)) project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  if (identical(revision, "locked")) revision <- lock$revision

  # JuliaConnectoR discovers Julia from JULIA_BINDIR/PATH. Its setup is
  # intentionally left to the user rather than triggering a download here.
  JuliaConnectoR::juliaEval("using Pkg")
  if (!is.null(project)) {
    cmd <- sprintf("Pkg.activate(%s); Pkg.resolve(); Pkg.instantiate()", .ctJuliaString(project))
  } else {
    spec <- sprintf("Pkg.PackageSpec(url=%s, rev=%s, subdir=%s)",
      .ctJuliaString(lock$url), .ctJuliaString(revision), .ctJuliaString(lock$subdir))
    cmd <- sprintf("Pkg.add(%s); Pkg.resolve(); Pkg.instantiate()", spec)
  }
  JuliaConnectoR::juliaEval(cmd)
  JuliaConnectoR::juliaEval("using ContinuousTimeSEM")
  .ct_julia_cache$project <- project
  .ct_julia_cache$revision <- revision
  .ct_julia_cache$module <- JuliaConnectoR::juliaImport("ContinuousTimeSEM")
  invisible(ctJuliaStatus())
}

#' Report Julia backend availability
#' @inheritParams ctJuliaSetup
#' @return A list describing the selected Julia engine.
#' @export
ctJuliaStatus <- function(project = NULL, julia_bin = NULL) {
  .ctJuliaRequire()
  julia_bin <- .ctJuliaBin(julia_bin)
  if (!is.null(julia_bin)) Sys.setenv(JULIA_BINDIR = julia_bin)
  available <- tryCatch({
    JuliaConnectoR::juliaEval("VERSION")
    TRUE
  }, error = function(e) FALSE)
  lock <- .ctJuliaEngineLock()
  list(available = available, project = .ctJuliaOr(project, .ct_julia_cache$project),
    revision = .ctJuliaOr(.ct_julia_cache$revision, lock$revision),
    lock = lock)
}

.ctJuliaModule <- function(project = NULL) {
  if (!is.null(.ct_julia_cache$module) && identical(project, .ct_julia_cache$project)) return(.ct_julia_cache$module)
  ctJuliaSetup(project = project)
  .ct_julia_cache$module
}

.ctJuliaClearSession <- function() {
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    try(JuliaConnectoR::stopJulia(), silent = TRUE)
  }
  .ct_julia_cache$module <- NULL
  .ct_julia_cache$project <- NULL
  .ct_julia_cache$revision <- NULL
  .ct_julia_cache$objectives <- new.env(parent = emptyenv())
  invisible(NULL)
}

.ctJuliaObjectiveKey <- function(spec) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Julia objective caching requires the suggested package digest.", call. = FALSE)
  }
  digest::digest(list(spec$parameter_table, spec$subject_starts, spec$times,
    spec$manifest_data, spec$tdpred_data, spec$tipred_data,
    spec$ti_effects, spec$max_timestep, spec$project, spec$engine), algo = "sha256")
}

.ctJuliaSubjectStarts <- function(ids) {
  if (!length(ids)) stop("The Julia backend requires at least one observation.", call. = FALSE)
  if (length(ids) == 1L) return(1L)
  c(1L, which(ids[-1L] != ids[-length(ids)]) + 1L)
}

.ctJuliaNumericVector <- function(values) {
  values <- as.numeric(values)
  if (!length(values)) return(JuliaConnectoR::juliaEval("Float64[]"))
  # JuliaConnectoR maps an R length-one numeric to a scalar.  A list keeps the
  # intended vector shape for Julia functions that require AbstractVector.
  JuliaConnectoR::juliaPut(as.list(values))
}

.ctJuliaInitialValues <- function(npar, inits = NULL) {
  # Match stanoptimis(): absent initial values are small, R-seeded draws in
  # unconstrained space rather than an exact all-zero vector.
  if (is.null(inits) || identical(inits, "random")) return(stats::rnorm(npar, 0, .01))
  values <- as.numeric(inits)
  if (length(values) != npar) {
    stop("Julia initial values must have one entry per free parameter.", call. = FALSE)
  }
  if (anyNA(values)) stop("Julia initial values must be numeric.", call. = FALSE)
  values
}

.ctJuliaUnsupported <- function(model, optimize, priors, intoverpop, vb, gendata,
  stanmodeltext, compileArgs, forcerecompile) {
  failures <- character()
  if (!isTRUE(optimize)) failures <- c(failures, "optimize=FALSE (HMC)")
  if (isTRUE(priors)) failures <- c(failures, "priors=TRUE")
  if (!isTRUE(model$continuoustime)) failures <- c(failures, "discrete-time model")
  if (any(model$manifesttype > 0)) failures <- c(failures, "non-Gaussian manifest variables")
  if (isTRUE(vb)) failures <- c(failures, "variational Bayes")
  if (isTRUE(gendata)) failures <- c(failures, "generation")
  if (!is.na(stanmodeltext)[1] || length(compileArgs) > 0L || isTRUE(forcerecompile)) failures <- c(failures, "Stan compilation controls")
  if (length(failures)) stop("Julia backend v1 does not support: ", paste(failures, collapse = ", "), ".", call. = FALSE)
}

.ctJuliaTDExpression <- function(expression) {
  # ctModelStatesAndPARS emits Stan-style row references. Julia receives the
  # current row as a vector in CTSEMRowContext instead.
  expression <- gsub("tdpreds\\[rowi\\s*,\\s*([0-9]+)\\]", "ctx.tdpreds[\\1]",
    expression, perl = TRUE)
  expression
}

.ctJuliaCanonicalModel <- function(model) {
  # ctFit has already established the canonical, augmented layout before it
  # dispatches to a backend. Reusing it keeps the state and parameter ordering
  # identical to Stan's matsetup/matvalues contract.
  if (!is.null(model$modelmats) && !is.null(model$intoverpop)) return(model)
  ctm <- ctsem:::ctModel0DRIFT(model, model$continuoustime)
  ctm$pars <- ctsem:::ctModelStatesAndPARS(
    ctm$pars, statenames = ctm$latentNames, tdprednames = ctm$TDpredNames
  )
  jacobian <- try(ctsem:::ctJacobian(ctm), silent = TRUE)
  if (inherits(jacobian, "try-error")) jacobian <- ctsem:::ctJacobian(ctm, simplify = FALSE)
  jacobian <- jacobian[names(ctsem:::ctStanMatricesList()$jacobian)]
  jacobian_rows <- ctsem:::ctModelUnlist(jacobian, names(jacobian))
  jacobian_rows <- jacobian_rows[apply(jacobian_rows, 1L, function(x) any(!is.na(x))), , drop = FALSE]
  if (nrow(jacobian_rows)) {
    template <- ctm$pars[rep(1L, nrow(jacobian_rows)), , drop = FALSE]
    for (name in intersect(c("matrix", "row", "col", "param", "value"), names(jacobian_rows))) {
      template[[name]] <- jacobian_rows[[name]]
    }
    # Jacobian entries derived from a named ctsem parameter inherit its
    # transform and free-coordinate mapping; state expressions are rewritten
    # below after the Jacobian rows have been appended.
    for (i in seq_len(nrow(template))) {
      source <- if (is.na(template$param[i])) NA_integer_ else match(template$param[i], ctm$pars$param)
      if (!is.na(source)) {
        for (name in setdiff(names(ctm$pars), c("matrix", "row", "col"))) {
          template[[name]][i] <- ctm$pars[[name]][source]
        }
      }
    }
    existing_key <- paste(ctm$pars$matrix, ctm$pars$row, ctm$pars$col, sep = "\r")
    template_key <- paste(template$matrix, template$row, template$col, sep = "\r")
    existing <- match(template_key, existing_key)
    for (i in seq_len(nrow(template))) {
      if (is.na(existing[i])) {
        ctm$pars <- rbind(ctm$pars, template[i, , drop = FALSE])
      } else {
        ctm$pars[existing[i], names(template)] <- template[i, names(template)]
      }
    }
  }
  ctm$pars <- ctsem:::ctModelStatesAndPARS(
    ctm$pars, statenames = ctm$latentNames, tdprednames = ctm$TDpredNames
  )
  ctsem:::T0VARredundancies(ctm)
}

.ctJuliaParameterTable <- function(model) {
  ctm <- .ctJuliaCanonicalModel(model)
  p <- as.data.frame(ctm$pars, stringsAsFactors = FALSE)
  p$matrix <- as.character(p$matrix)
  render_transform <- function(transform, multiplier, meanscale, offset, inneroffset) {
    transform <- as.integer(transform)
    if (transform == 0L) {
      return(sprintf("%.17g + %.17g * (param * %.17g + %.17g)",
        offset, multiplier, meanscale, inneroffset))
    }
    ctsem:::tform("param", transform, multiplier, meanscale, offset,
      inneroffset, singletext = TRUE)
  }
  numeric_transform <- !is.na(suppressWarnings(as.integer(p$transform)))
  for (i in which(numeric_transform)) {
    rendered <- render_transform(p$transform[i], p$multiplier[i],
      p$meanscale[i], p$offset[i], p$inneroffset[i])
    if (length(rendered) != 1L) {
      stop("Julia backend cannot render scalar transform ", p$transform[i],
        " for ", p$matrix[i], "[", p$row[i], ",", p$col[i], "].", call. = FALSE)
    }
    p$transform[[i]] <- rendered[[1L]]
  }
  p$parnumber <- NA_integer_
  p$predicttransform <- NA_character_
  p$updatetransform <- NA_character_
  free <- !is.na(p$param) & !grepl("[", p$param, fixed = TRUE) &
    is.na(suppressWarnings(as.numeric(p$value)))
  parnames <- unique(as.character(p$param[free]))
  p$parnumber[free] <- match(as.character(p$param[free]), parnames)
  if (!is.null(ctm$modelmats$matsetup)) {
    matrix_codes <- ctsem:::ctStanMatricesList()$all
    setup <- as.data.frame(ctm$modelmats$matsetup)
    values <- as.data.frame(ctm$modelmats$matvalues)
    setup$matrix_name <- names(matrix_codes)[match(setup$matrix, matrix_codes)]
    setup_key <- paste(setup$matrix_name, setup$row, setup$col, sep = "\r")
    parameter_key <- paste(p$matrix, p$row, p$col, sep = "\r")
    setup_row <- match(parameter_key, setup_key)
    mapped <- setup$param[setup_row]
    use_mapped <- free & !is.na(mapped) & mapped > 0L
    p$parnumber[use_mapped] <- as.integer(mapped[use_mapped])

    # Stan transforms each raw population coordinate once, using the first
    # matching base/PARS setup row. Reusing a parameter in another matrix must
    # not reapply that cell's local transform.
    setup_parameter <- as.integer(setup$param)
    setup_when <- as.integer(setup$when)
    candidate <- which(setup_parameter > 0L & setup_when %in% c(0L, 100L))
    canonical <- rep(NA_character_, max(c(0L, setup_parameter), na.rm = TRUE))
    for (row in candidate) {
      parameter <- setup_parameter[row]
      if (!is.na(canonical[parameter])) next
      canonical[parameter] <- render_transform(setup$transform[row],
        values$multiplier[row], values$meanscale[row], values$offset[row],
        values$inneroffset[row])
    }
    use_canonical <- free & !is.na(p$parnumber) & p$parnumber <= length(canonical) &
      !is.na(canonical[p$parnumber])
    p$transform[use_canonical] <- canonical[p$parnumber[use_canonical]]
  }
  dynamic <- !is.na(p$param) & grepl("[", p$param, fixed = TRUE)
  predict_matrices <- c("PARS", "DRIFT", "CINT", "DIFFUSION", "JAx")
  update_matrices <- c("LAMBDA", "MANIFESTMEANS", "MANIFESTVAR", "Jy")
  td_matrices <- c("TDPREDEFFECT", "Jtd")
  p$predicttransform[dynamic & p$matrix %in% predict_matrices] <- .ctJuliaTDExpression(p$param[dynamic & p$matrix %in% predict_matrices])
  p$updatetransform[dynamic & p$matrix %in% update_matrices] <- .ctJuliaTDExpression(p$param[dynamic & p$matrix %in% update_matrices])
  p$tdtransform[dynamic & p$matrix %in% td_matrices] <- .ctJuliaTDExpression(p$param[dynamic & p$matrix %in% td_matrices])
  unsupported_dynamic <- dynamic & !(p$matrix %in% c(predict_matrices, update_matrices, td_matrices))
  if (any(unsupported_dynamic)) {
    bad <- p[which(unsupported_dynamic)[1L], c("matrix", "row", "col"), drop = FALSE]
    stop("Julia backend does not support a state-dependent expression in ",
      bad$matrix, "[", bad$row, ",", bad$col, "].", call. = FALSE)
  }
  p$param[dynamic] <- NA_character_
  bare <- free & is.na(p$transform)
  p$transform[bare] <- paste0("param[", p$parnumber[bare], "]")
  for (i in which(free & !is.na(p$transform))) {
    p$transform[i] <- gsub("\\bparam\\b", paste0("param[", p$parnumber[i], "]"),
      as.character(p$transform[i]), perl = TRUE)
  }
  # Julia consumes these canonical columns directly.  Preserve absent optional
  # transforms as explicit NA columns: omitting them changes the DataFrame
  # schema and prevents a serialised specification from being reconstructed.
  keep <- c("matrix", "row", "col", "param", "parnumber", "value", "transform",
    "predicttransform", "updatetransform", "tdtransform", "indvarying",
    grep("_effect$", names(p), value = TRUE))
  for (name in setdiff(keep, names(p))) p[[name]] <- NA
  out <- p[, keep, drop = FALSE]
  # The buffered Julia kernels always carry a PARS component, even for models
  # without user-declared state parameters.  Give such models one fixed dummy
  # entry so the component axis is stable and complex-transform calls remain
  # well-defined.
  if (!"PARS" %in% out$matrix) {
    dummy <- out[1L, , drop = FALSE]
    dummy[] <- NA
    dummy$matrix <- "PARS"
    dummy$row <- 1L
    dummy$col <- 1L
    dummy$value <- 0
    dummy$indvarying <- FALSE
    out <- rbind(out, dummy)
  }
  # Complex transforms can use PARS values. Keep that component first so its
  # state-derived entries are refreshed before dependent DRIFT/JAc transforms.
  matrix_order <- c("PARS", unique(as.character(out$matrix[out$matrix != "PARS"])))
  out <- out[order(match(out$matrix, matrix_order)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.ctJuliaPadMatrix <- function(table, matrix, nrow, ncol) {
  present <- table$matrix == matrix
  if (!any(present)) return(table)
  template <- table[which(present)[1L], , drop = FALSE]
  for (row in seq_len(nrow)) for (col in seq_len(ncol)) {
    if (any(table$matrix == matrix & table$row == row & table$col == col)) next
    entry <- template
    entry[] <- NA
    entry$matrix <- matrix
    entry$row <- row
    entry$col <- col
    entry$value <- 0
    entry$indvarying <- FALSE
    effect_columns <- grep("_effect$", names(entry), value = TRUE)
    if (length(effect_columns)) entry[effect_columns] <- FALSE
    table <- rbind(table, entry)
  }
  table
}

.ctJuliaAugmentRandomEffects <- function(model) {
  original_nlatent <- model$n.latent
  prepared_augmentation <- !is.null(model$intoverpopindvaryingindex)
  has_random_effects <- prepared_augmentation || any(model$pars$indvarying %in% TRUE &
    is.na(suppressWarnings(as.numeric(model$pars$value))), na.rm = TRUE)
  if (!has_random_effects) {
    return(list(parameter_table = .ctJuliaParameterTable(model),
      nlatent = original_nlatent, nlatent_augmented = original_nlatent,
      dynamic_state_indices = seq_len(original_nlatent),
      random_effects = data.frame(), rewritten_cells = data.frame()))
  }

  expanded <- if (prepared_augmentation) model else {
    prepared <- ctsem:::ctModel0DRIFT(model, model$continuoustime)
    prepared$pars <- ctsem:::ctModelStatesAndPARS(prepared$pars,
      statenames = prepared$latentNames, tdprednames = prepared$TDpredNames)
    ctsem:::ctStanModelIntOverPop(prepared)
  }
  augmented_indices <- as.integer(expanded$intoverpopindvaryingindex)
  nlatent_augmented <- max(expanded$pars$row[expanded$pars$matrix == "T0MEANS"])
  table <- .ctJuliaParameterTable(expanded)

  for (matrix in c("DRIFT", "DIFFUSION", "JAx", "Jtd")) {
    table <- .ctJuliaPadMatrix(table, matrix, nlatent_augmented, nlatent_augmented)
  }
  for (matrix in c("CINT", "T0MEANS")) {
    table <- .ctJuliaPadMatrix(table, matrix, nlatent_augmented, 1L)
  }
  for (matrix in c("LAMBDA", "Jy")) {
    if (matrix %in% table$matrix) {
      table <- .ctJuliaPadMatrix(table, matrix,
        max(table$row[table$matrix == matrix]), nlatent_augmented)
    }
  }
  if ("TDPREDEFFECT" %in% table$matrix) {
    table <- .ctJuliaPadMatrix(table, "TDPREDEFFECT", nlatent_augmented,
      max(table$col[table$matrix == "TDPREDEFFECT"]))
  }
  table <- .ctJuliaPadMatrix(table, "T0VAR", nlatent_augmented, nlatent_augmented)

  next_parameter <- max(table$parnumber, na.rm = TRUE)
  random_sd_scale <- rep(1, length(augmented_indices))
  if (!is.null(expanded$modelmats$matsetup)) {
    setup <- as.data.frame(expanded$modelmats$matsetup)
    values <- as.data.frame(expanded$modelmats$matvalues)
    varying_parameters <- unique(setup$param[setup$indvarying > 0L])
    varying_parameters <- varying_parameters[varying_parameters > 0L]
    random_sd_scale <- values$sdscale[match(varying_parameters, setup$param)]
    random_sd_scale[is.na(random_sd_scale)] <- 1
  }
  if (length(random_sd_scale) != length(augmented_indices)) {
    stop("Prepared random-effect covariance metadata does not match the augmented state layout.", call. = FALSE)
  }
  covariance_rows <- list()
  # Match Stan's unconstrained parameter order exactly: all population scales,
  # then lower-triangular correlation coordinates column by column.
  for (position in seq_along(augmented_indices)) {
    row <- augmented_indices[position]
    col <- row
    index <- which(table$matrix == "T0VAR" & table$row == row & table$col == col)
    length(index) == 1L || stop("Internal Julia augmentation error: missing T0VAR entry.", call. = FALSE)
    next_parameter <- next_parameter + 1L
    table$param[index] <- sprintf("julia_popcov_%d_%d", row, col)
    table$parnumber[index] <- next_parameter
    table$value[index] <- NA_real_
    table$transform[index] <- sprintf("1e-10 + %.17g * log1p_exp(2 * param[%d] - 1)",
      random_sd_scale[position], next_parameter)
    covariance_rows[[length(covariance_rows) + 1L]] <- data.frame(
      row = row, col = col, parameter = next_parameter,
      type = "sd"
    )
  }
  if (length(augmented_indices) > 1L) for (column_position in seq_len(length(augmented_indices) - 1L)) {
    for (row_position in (column_position + 1L):length(augmented_indices)) {
      row <- augmented_indices[row_position]
      col <- augmented_indices[column_position]
      index <- which(table$matrix == "T0VAR" & table$row == row & table$col == col)
      length(index) == 1L || stop("Internal Julia augmentation error: missing T0VAR entry.", call. = FALSE)
      next_parameter <- next_parameter + 1L
      table$param[index] <- sprintf("julia_popcov_%d_%d", row, col)
      table$parnumber[index] <- next_parameter
      table$value[index] <- NA_real_
      table$transform[index] <- sprintf("2 / (1 + exp(-param[%d])) - 1", next_parameter)
      covariance_rows[[length(covariance_rows) + 1L]] <- data.frame(
        row = row, col = col, parameter = next_parameter, type = "correlation"
      )
    }
  }
  predict_transform <- replace(table$predicttransform, is.na(table$predicttransform), "")
  update_transform <- replace(table$updatetransform, is.na(table$updatetransform), "")
  td_transform <- replace(table$tdtransform, is.na(table$tdtransform), "")
  rewritten <- table[grepl("state\\[", predict_transform) |
    grepl("state\\[", update_transform) |
    grepl("state\\[", td_transform),
    c("matrix", "row", "col"), drop = FALSE]
  list(parameter_table = table, nlatent = original_nlatent,
    nlatent_augmented = nlatent_augmented,
    dynamic_state_indices = setdiff(seq_len(nlatent_augmented), augmented_indices),
    random_effects = do.call(rbind, covariance_rows), rewritten_cells = rewritten)
}

.ctJuliaTIData <- function(dat, model) {
  subject_ids <- unique(dat[[model$subjectIDname]])
  if (!model$n.TIpred) return(matrix(numeric(), nrow = length(subject_ids), ncol = 0L))
  values <- matrix(NA_real_, nrow = length(subject_ids), ncol = model$n.TIpred,
    dimnames = list(NULL, model$TIpredNames))
  for (i in seq_along(subject_ids)) {
    rows <- dat[[model$subjectIDname]] == subject_ids[i]
    subject_values <- as.matrix(dat[rows, model$TIpredNames, drop = FALSE])
    if (anyNA(subject_values)) {
      stop("Julia backend currently requires complete TI predictors; missing values must be handled before fitting.", call. = FALSE)
    }
    first <- subject_values[1L, ]
    if (any(vapply(seq_len(ncol(subject_values)), function(j) {
      any(subject_values[, j] != first[j])
    }, logical(1)))) {
      stop("Julia backend requires TI predictors to be constant within subject.", call. = FALSE)
    }
    values[i, ] <- first
  }
  values
}

.ctJuliaValidateTIConstancy <- function(dat, model) {
  if (!model$n.TIpred) return(invisible(NULL))
  for (subject in unique(dat[[model$subjectIDname]])) {
    rows <- dat[[model$subjectIDname]] == subject
    values <- dat[rows, model$TIpredNames, drop = FALSE]
    for (name in model$TIpredNames) {
      observed <- unique(values[[name]][!is.na(values[[name]])])
      if (length(observed) > 1L) {
        stop("Julia backend requires TI predictors to be constant within subject.", call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

.ctJuliaTIEffects <- function(table, model) {
  if (!model$n.TIpred) {
    return(data.frame(parameter = integer(), predictor = integer(), coefficient = integer()))
  }
  effect_columns <- paste0(model$TIpredNames, "_effect")
  available <- intersect(effect_columns, names(table))
  if (!length(available)) {
    return(data.frame(parameter = integer(), predictor = integer(), coefficient = integer()))
  }
  direct <- !is.na(table$parnumber) & !grepl("[", table$param, fixed = TRUE)
  entries <- list()
  coefficient <- max(table$parnumber, na.rm = TRUE)
  for (predictor in seq_along(model$TIpredNames)) {
    column <- effect_columns[predictor]
    if (!column %in% available) next
    parameters <- sort(unique(table$parnumber[direct & (table[[column]] %in% TRUE)]))
    for (parameter in parameters) {
      coefficient <- coefficient + 1L
      entries[[length(entries) + 1L]] <- data.frame(
        parameter = as.integer(parameter), predictor = as.integer(predictor),
        coefficient = as.integer(coefficient)
      )
    }
  }
  if (!length(entries)) return(data.frame(parameter = integer(), predictor = integer(), coefficient = integer()))
  do.call(rbind, entries)
}

.ctJuliaPrepare <- function(datalong, model, prepared_data = NULL, project = NULL) {
  dat <- data.frame(datalong)
  dat <- dat[order(dat[[model$subjectIDname]], dat[[model$timeName]]), , drop = FALSE]
  .ctJuliaValidateTIConstancy(dat, model)
  subject_starts <- .ctJuliaSubjectStarts(dat[[model$subjectIDname]])
  tdpred_data <- if (is.null(prepared_data)) {
    if (model$n.TDpred) as.matrix(dat[, model$TDpredNames, drop = FALSE]) else matrix(numeric(), nrow(dat), 0L)
  } else prepared_data$tdpreds
  tipred_data <- if (is.null(prepared_data)) .ctJuliaTIData(dat, model) else prepared_data$tipredsdata
  if (is.null(tdpred_data)) tdpred_data <- matrix(numeric(), nrow(dat), 0L)
  if (is.null(tipred_data)) tipred_data <- matrix(numeric(), nrow = length(unique(dat[[model$subjectIDname]])), ncol = 0L)
  tipred_data <- as.matrix(tipred_data)
  if (!ncol(tipred_data)) tipred_data <- matrix(numeric(), nrow = length(subject_starts), ncol = 0L)
  if (nrow(tdpred_data) != nrow(dat)) stop("Prepared TD predictor rows do not match the fitted data.", call. = FALSE)
  if (nrow(tipred_data) != length(subject_starts)) stop("Prepared TI predictor rows do not match the fitted subjects.", call. = FALSE)
  max_timestep <- if (!is.null(prepared_data$maxtimestep)) {
    as.numeric(prepared_data$maxtimestep)[1L]
  } else if (!is.null(model$nlcontrol$maxtimestep)) {
    as.numeric(model$nlcontrol$maxtimestep)[1L]
  } else 999999
  if (!is.finite(max_timestep) || max_timestep <= 0) stop("Julia maxtimestep must be a positive finite number.", call. = FALSE)
  augmented <- .ctJuliaAugmentRandomEffects(model)
  parameter_table <- augmented$parameter_table
  list(
    class = "ctJuliaModel",
    model = model,
    data = dat,
    parameter_table = parameter_table,
    subject_starts = as.integer(subject_starts),
    times = as.numeric(dat[[model$timeName]]),
    manifest_data = t(as.matrix(dat[, model$manifestNames, drop = FALSE])),
    tdpred_data = t(tdpred_data),
    tipred_data = as.matrix(tipred_data),
    ti_effects = .ctJuliaTIEffects(parameter_table, model),
    max_timestep = max_timestep,
    TDpredNames = model$TDpredNames,
    TIpredNames = model$TIpredNames,
    nlatent = augmented$nlatent,
    nlatent_augmented = augmented$nlatent_augmented,
    dynamic_state_indices = augmented$dynamic_state_indices,
    random_effects = augmented$random_effects,
    rewritten_cells = augmented$rewritten_cells,
    project = project,
    engine = .ctJuliaEngineLock()
  )
}

.ctJuliaObjective <- function(object) {
  stopifnot(inherits(object, "ctJuliaModel") || inherits(object, "ctJuliaFit"))
  spec <- if (inherits(object, "ctJuliaFit")) object$model_spec else object
  key <- .ctJuliaObjectiveKey(spec)
  if (exists(key, envir = .ct_julia_cache$objectives, inherits = FALSE)) {
    return(get(key, envir = .ct_julia_cache$objectives, inherits = FALSE))
  }
  module <- .ctJuliaModule(spec$project)
  # DataFrames is a dependency of ContinuousTimeSEM, not necessarily of
  # Julia's active global environment. Construct through the engine module.
  df <- module$DataFrame(spec$parameter_table)
  effects <- module$DataFrame(spec$ti_effects)
  params <- module$ekf_from_data_frame(df, effects,
    .ctJuliaNumericVector(spec$dynamic_state_indices))
  objective <- module$ctsem_objective(params, JuliaConnectoR::juliaPut(spec$subject_starts),
    JuliaConnectoR::juliaPut(spec$times), JuliaConnectoR::juliaPut(spec$manifest_data),
    JuliaConnectoR::juliaPut(spec$tdpred_data), JuliaConnectoR::juliaPut(spec$tipred_data),
    spec$max_timestep)
  assign(key, objective, envir = .ct_julia_cache$objectives)
  objective
}

#' Evaluate a prepared Julia ctsem likelihood
#' @param object A \code{ctJuliaModel} or \code{ctJuliaFit}.
#' @param pars Unconstrained parameters; defaults to the fitted estimate.
#' @param gradient Return an analytic/automatic-differentiation gradient.
#' @param contributions Return subject and row contribution details when available.
#' @return A list containing log likelihood and, when requested, gradient.
#' @export
ctJuliaEvaluate <- function(object, pars = NULL, gradient = TRUE, contributions = FALSE) {
  if (!inherits(object, c("ctJuliaModel", "ctJuliaFit"))) stop("object must be a ctJuliaModel or ctJuliaFit", call. = FALSE)
  if (is.null(pars)) {
    if (inherits(object, "ctJuliaFit")) pars <- object$estimate$raw else stop("pars must be supplied for a prepared ctJuliaModel", call. = FALSE)
  }
  module <- .ctJuliaModule(if (inherits(object, "ctJuliaFit")) object$model_spec$project else object$project)
  result <- module$ctsem_evaluate(.ctJuliaObjective(object), .ctJuliaNumericVector(pars),
    gradient = isTRUE(gradient), contributions = isTRUE(contributions))
  JuliaConnectoR::juliaGet(result)
}

#' @export
ctExtract.ctJuliaFit <- function(object, subjectMatrices = FALSE, cores = 2,
  nsamples = "all", subjects = "all", ...) {
  if (isTRUE(subjectMatrices)) {
    stop("Subject matrices are not yet available for ctJuliaFit objects.", call. = FALSE)
  }
  if (!identical(nsamples, "all") || !identical(subjects, "all")) {
    stop("ctJuliaFit contains a point estimate, not posterior samples.", call. = FALSE)
  }
  list(rawpars = object$estimate$raw, loglik = object$estimate$loglik,
    gradient = object$estimate$gradient,
    subject_loglik = object$estimate$subject_loglik)
}

#' @export
ctSummaryMatrices.ctJuliaFit <- function(fit, ...) {
  stop("ctSummaryMatrices() for ctJuliaFit requires the Julia parameter-matrix reconstruction API, which is not yet implemented.", call. = FALSE)
}

ctFitJuliaBackend <- function(datalong, model, prepared_data = NULL, inits = NULL, cores = 1L,
  backendcontrol = list(), verbose = 0L, fit = TRUE) {
  if (isTRUE(backendcontrol$restart_session)) .ctJuliaClearSession()
  project <- .ctJuliaOr(backendcontrol$julia_project, NULL)
  gradient <- .ctJuliaOr(backendcontrol$gradient, "forward")
  if (!gradient %in% c("forward", "adjoint")) stop("backendcontrol$gradient must be 'forward' or 'adjoint'", call. = FALSE)
  if (identical(gradient, "adjoint")) stop("The main-based Julia engine has not yet validated its adjoint implementation; use backendcontrol=list(gradient='forward').", call. = FALSE)
  model_spec <- .ctJuliaPrepare(datalong, model, prepared_data = prepared_data, project = project)
  if (!fit) return(structure(model_spec, class = c("ctJuliaModel", "ctFitModel")))

  objective <- .ctJuliaObjective(structure(model_spec, class = c("ctJuliaModel", "ctFitModel")))
  npar <- max(c(model_spec$parameter_table$parnumber, model_spec$ti_effects$coefficient), na.rm = TRUE)
  start <- .ctJuliaInitialValues(npar, inits)
  module <- .ctJuliaModule(project)
  result <- JuliaConnectoR::juliaGet(module$ctsem_optimize(objective, .ctJuliaNumericVector(start),
    maxiter = as.integer(.ctJuliaOr(backendcontrol$maxiter, 1000L)),
    g_tol = .ctJuliaOr(backendcontrol$g_tol, 1e-8),
    f_tol = .ctJuliaOr(backendcontrol$f_tol, 0),
    x_tol = .ctJuliaOr(backendcontrol$x_tol, 0),
    verbose = verbose > 0L))
  out <- list(backend = "julia", model = model, model_spec = model_spec,
    data = datalong, estimate = list(raw = as.numeric(result$minimizer),
      loglik = as.numeric(result$maximum_loglik), gradient = as.numeric(result$gradient),
      subject_loglik = result$subject_loglik, converged = isTRUE(result$converged),
      iterations = as.integer(result$iterations)), engine = model_spec$engine,
    args = list(backend = "julia", backendcontrol = backendcontrol, cores = cores))
  class(out) <- c("ctJuliaFit", "ctFit")
  out
}

#' @export
print.ctJuliaFit <- function(x, ...) {
  cat("ctsem Julia fit\n")
  cat("  log likelihood:", format(x$estimate$loglik), "\n")
  cat("  converged:", x$estimate$converged, " iterations:", x$estimate$iterations, "\n")
  invisible(x)
}

#' @export
coef.ctJuliaFit <- function(object, ...) object$estimate$raw

#' @export
logLik.ctJuliaFit <- function(object, ...) {
  structure(object$estimate$loglik, df = length(object$estimate$raw), nobs = nrow(object$data), class = "logLik")
}

#' @export
summary.ctJuliaFit <- function(object, ...) {
  parameter_table <- object$model_spec$parameter_table
  free_rows <- parameter_table[!is.na(parameter_table$parnumber), , drop = FALSE]
  free_rows <- free_rows[match(seq_len(length(object$estimate$raw)), free_rows$parnumber), , drop = FALSE]
  list(backend = "julia", loglik = object$estimate$loglik,
    coefficients = stats::setNames(object$estimate$raw, free_rows$param),
    gradient = object$estimate$gradient, converged = object$estimate$converged,
    iterations = object$estimate$iterations, note = "Point estimates only; uncertainty methods are not available for Julia fits.")
}
