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
  # Backslashes must be escaped, not merely tolerated: a Windows path pasted
  # into a Julia string literal otherwise becomes an escape sequence, and
  # a path like ...\Users\... fails to parse as "invalid unicode escape"
  # rather than as anything that points at the real problem.
  value <- gsub("\\", "\\\\", value, fixed = TRUE)
  paste0('"', value, '"')
}

# Where Julia is, in order of authority: what the caller named, what the user
# configured, what ctsem installed itself, what is on the PATH, and finally
# juliaup's own installs -- which R started from a launcher rather than a shell
# will not have inherited a PATH entry for. NULL means "nowhere we can see",
# which is what makes ctJuliaInstall() offer to fetch one.
.ctJuliaBin <- function(julia_bin = NULL) {
  if (!is.null(julia_bin)) {
    return(normalizePath(julia_bin, winslash = "/", mustWork = TRUE))
  }
  configured <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (.ctJuliaIsBinDir(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  managed <- .ctJuliaManagedBin()
  if (!is.null(managed)) return(managed)
  onpath <- Sys.which("julia")[[1]]
  if (nzchar(onpath)) {
    return(normalizePath(dirname(onpath), winslash = "/", mustWork = TRUE))
  }
  .ctJuliaJuliaupBin()
}

.ctJuliaRequire <- function() {
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) return(invisible(TRUE))
  # Offered rather than demanded: this is one CRAN package away from working,
  # and sending the user off to install it by hand and come back is the friction
  # ctJuliaInstall() exists to remove.
  if (isTRUE(tryCatch(.ctJuliaInstallConnectoR(), error = function(e) FALSE))) {
    return(invisible(TRUE))
  }
  stop("backend='julia' needs the JuliaConnectoR package, which is not installed.\n",
    .ctJuliaDeclined(), call. = FALSE)
}

.ctJuliaEngineLock <- function() {
  path <- .ct_julia_engine_file()
  if (!nzchar(path) || !file.exists(path)) {
    stop("ctsem Julia engine lock file is unavailable; reinstall ctsem.", call. = FALSE)
  }
  json <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # Provenance for the vendored engine in inst/julia/, not an install spec: the
  # url/branch/revision say which upstream commit the vendored copy was taken
  # from, so the two can be compared. tools/sync-julia-engine.sh writes it.
  fields <- c("url", "branch", "revision", "subdir", "julia")
  out <- lapply(fields, function(field) {
    hit <- regmatches(json, regexec(paste0('"', field, '"\\s*:\\s*"([^"]+)"'), json))[[1]]
    if (length(hit) < 2L) stop("Malformed ctsem Julia engine lock file.", call. = FALSE)
    hit[[2]]
  })
  names(out) <- fields
  out
}

# Path to the copy of ContinuousTimeSEM.jl shipped inside this ctsem install.
.ctJuliaEnginePath <- function() {
  path <- system.file("julia", "ContinuousTimeSEM", package = "ctsem")
  if (!nzchar(path) || !file.exists(file.path(path, "Project.toml"))) {
    stop("The Julia engine is missing from this ctsem installation; reinstall ctsem.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# A writable project directory for the engine, keyed by the vendored revision so
# a ctsem upgrade gets a fresh environment instead of reusing a stale manifest.
#
# The vendored tree is *copied* here rather than activated in place: activating a
# project writes to its Manifest.toml, and an R library directory is frequently
# read-only. The copy is under 400 KB.
.ctJuliaEnvDir <- function(lock) {
  base <- tools::R_user_dir("ctsem", which = "cache")
  # Forward slashes throughout: this path is used by R and also embedded in
  # Julia source, and Julia accepts them on every platform.
  gsub("\\", "/", file.path(base, "julia", paste0("engine-", substr(lock$revision, 1, 12))),
    fixed = TRUE)
}

# Whether a Julia session already exists. JuliaConnectoR starts one lazily on
# the first call, so "has anything talked to Julia yet" is the question.
.ctJuliaSessionRunning <- function() {
  !is.null(.ct_julia_cache$module) ||
    isTRUE(tryCatch(JuliaConnectoR::juliaEval("true"), error = function(e) FALSE))
}

.ctJuliaCheckAvailable <- function() {
  ok <- tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)
  if (isTRUE(ok)) return(invisible(TRUE))
  # Same reasoning as .ctJuliaRequire(): ask, rather than end the session's work
  # with an instruction to install something and start again.
  if (isTRUE(tryCatch(.ctJuliaOfferJulia(), error = function(e) FALSE)) &&
      isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE))) {
    return(invisible(TRUE))
  }
  stop("Julia was not found, and backend='julia' needs it (", .ct_julia_minimum, " or newer).\n",
    .ctJuliaDeclined(), "\n",
    "  ctJuliaInstall() downloads the official build into ", .ctJuliaInstallRoot(), ".\n",
    "  To use a Julia you already have, set JULIA_BINDIR, e.g.\n",
    "    Sys.setenv(JULIA_BINDIR = \"/path/to/julia/bin\")",
    call. = FALSE)
}

#' Configure the Julia engine used by ctsem
#'
#' Prepares the copy of ContinuousTimeSEM.jl that ships inside this ctsem
#' installation. No network access and no repository credentials are involved:
#' the engine source is vendored in \code{inst/julia/}, and this only creates a
#' Julia project for it and instantiates its dependencies.
#'
#' The first call downloads and precompiles those dependencies (roughly 120 MB
#' and a minute or two); later calls in new sessions reuse them.
#'
#' If \pkg{JuliaConnectoR} or Julia itself is missing, this offers to install it
#' rather than failing -- the same thing \code{\link{ctJuliaInstall}} does, which
#' is the function to reach for when setting the backend up deliberately, or
#' from a script.
#' @param project Optional local ContinuousTimeSEM.jl checkout to use instead of
#'   the vendored copy, for engine development.
#' @param revision Ignored; retained for backward compatibility. The engine
#'   revision is whatever is vendored, and is reported by \code{ctJuliaStatus()}.
#' @param julia_bin Optional Julia binary directory.
#' @param threads Number of Julia threads. The engine splits its subject loop
#'   across them. Julia fixes its thread count at process start, so this only
#'   takes effect if no Julia session is running yet -- pass \code{force = TRUE}
#'   to restart one. \code{NULL} leaves it to Julia's own default (one thread
#'   unless \code{JULIA_NUM_THREADS} is already set).
#' @param force Reconfigure an existing Julia session.
#' @return A Julia-engine status list, invisibly.
#' @seealso \code{\link{ctJuliaInstall}}, \code{\link{ctJuliaStatus}}
#' @export
ctJuliaSetup <- function(project = NULL, revision = "locked", julia_bin = NULL,
  threads = NULL, force = FALSE) {
  .ctJuliaRequire()
  if (isTRUE(force)) .ctJuliaClearSession()
  julia_bin <- .ctJuliaBin(julia_bin)
  if (!is.null(julia_bin)) Sys.setenv(JULIA_BINDIR = julia_bin)
  if (!is.null(threads)) {
    threads <- max(1L, as.integer(threads)[1L])
    # Julia fixes Threads.nthreads() at process start and JuliaConnectoR's
    # subprocess inherits this environment, so it has to be set before the
    # session exists. Warning rather than silently doing nothing matters here:
    # a user asking for 8 threads and getting 1 would otherwise just see a
    # disappointing benchmark.
    if (.ctJuliaSessionRunning()) {
      if (!identical(Sys.getenv("JULIA_NUM_THREADS", unset = ""), as.character(threads))) {
        warning("A Julia session is already running, so threads=", threads,
          " has no effect. Use ctJuliaSetup(threads=", threads,
          ", force=TRUE) to restart it.", call. = FALSE)
      }
    } else {
      Sys.setenv(JULIA_NUM_THREADS = as.character(threads))
    }
  }
  .ctJuliaCheckAvailable()
  lock <- .ctJuliaEngineLock()

  if (!is.null(project)) {
    # Developer override: use the checkout as its own project, in place.
    project <- normalizePath(project, winslash = "/", mustWork = TRUE)
    env_dir <- project
  } else {
    env_dir <- .ctJuliaEnvDir(lock)
    if (!file.exists(file.path(env_dir, "Project.toml"))) {
      dir.create(dirname(env_dir), recursive = TRUE, showWarnings = FALSE)
      unlink(env_dir, recursive = TRUE)
      # copy the vendored tree, then rename to the target so an interrupted
      # copy cannot leave a half-populated environment behind
      staging <- paste0(env_dir, "-partial")
      unlink(staging, recursive = TRUE)
      dir.create(staging, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(.ctJuliaEnginePath(), full.names = TRUE), staging,
        recursive = TRUE)
      file.rename(staging, env_dir)
    }
  }

  JuliaConnectoR::juliaEval("using Pkg")
  activate <- sprintf("Pkg.activate(%s)", .ctJuliaString(env_dir))
  instantiated <- tryCatch({
    JuliaConnectoR::juliaEval(paste0(activate, "; Pkg.instantiate()"))
    TRUE
  }, error = function(e) FALSE)
  if (!instantiated) {
    # The vendored manifest pins the versions this ctsem release was tested
    # against, but it can be unsatisfiable on a different Julia version. Falling
    # back to a fresh resolve is better than refusing to run; the compat bounds
    # in Project.toml still apply.
    unlink(file.path(env_dir, "Manifest.toml"))
    JuliaConnectoR::juliaEval(paste0(activate, "; Pkg.resolve(); Pkg.instantiate()"))
  }
  JuliaConnectoR::juliaEval("using ContinuousTimeSEM")
  .ct_julia_cache$project <- project
  .ct_julia_cache$revision <- lock$revision
  .ct_julia_cache$module <- JuliaConnectoR::juliaImport("ContinuousTimeSEM")
  invisible(ctJuliaStatus())
}

#' Report Julia backend availability
#'
#' Reports what \code{backend='julia'} would find, and installs nothing: it is
#' safe to call on a machine with no Julia and no \pkg{JuliaConnectoR}, where it
#' says so rather than erroring. \code{\link{ctJuliaInstall}} supplies whatever
#' it reports as missing.
#'
#' @inheritParams ctJuliaSetup
#' @return A list describing the selected Julia engine: whether it is
#'   \code{available}, whether the \code{connectoR} bridge package is installed,
#'   the \code{julia_bin} directory in use and the \code{julia} version there,
#'   the engine \code{revision}, the number of \code{threads} in a running
#'   session, and the engine provenance \code{lock}.
#' @export
ctJuliaStatus <- function(project = NULL, julia_bin = NULL) {
  connectoR <- requireNamespace("JuliaConnectoR", quietly = TRUE)
  julia_bin <- .ctJuliaBin(julia_bin)
  if (!is.null(julia_bin)) Sys.setenv(JULIA_BINDIR = julia_bin)
  # A running session can be asked its own version; only a Julia that has never
  # started needs a subprocess spawned to find out.
  version <- if (connectoR) {
    tryCatch(as.character(JuliaConnectoR::juliaEval("string(VERSION)")),
      error = function(e) NA_character_)
  } else NA_character_
  available <- !is.na(version)
  if (!available) version <- .ctJuliaBinVersion(julia_bin)
  lock <- .ctJuliaEngineLock()
  threads <- if (available) {
    tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
      error = function(e) NA_integer_)
  } else NA_integer_
  list(available = available, connectoR = connectoR,
    julia_bin = .ctJuliaOr(julia_bin, NA_character_), julia = version,
    project = .ctJuliaOr(project, .ct_julia_cache$project),
    revision = .ctJuliaOr(.ct_julia_cache$revision, lock$revision),
    threads = threads, lock = lock)
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
    spec$ti_effects, spec$priors, spec$max_timestep, spec$project, spec$engine),
    algo = "sha256")
}

.ctJuliaSubjectStarts <- function(ids) {
  if (!length(ids)) stop("The Julia backend requires at least one observation.", call. = FALSE)
  if (length(ids) == 1L) return(1L)
  c(1L, which(ids[-1L] != ids[-length(ids)]) + 1L)
}

.ctJuliaNumericVector <- function(values) {
  values <- as.numeric(values)
  if (!length(values)) return(JuliaConnectoR::juliaEval("Float64[]"))
  # JuliaConnectoR maps an R length-one numeric to a scalar, so a length-one
  # vector has to go across as a list to arrive as an AbstractVector. Every
  # other length must NOT: `juliaPut` marshals a plain numeric vector as a
  # binary block and a list element by element, and both arrive as the same
  # `Vector{Float64}`. Measured on this machine, marshalling the parameter
  # vector of a 1490-parameter model took 0.069 s as a list against 0.0002 s
  # as a vector -- 344x, and *half the total wall time of a gradient call*,
  # which made every published Julia backend timing mostly JuliaConnectoR
  # rather than Julia. This function is called once per objective evaluation,
  # so it sits directly in the optimizer's inner loop.
  if (length(values) == 1L) return(JuliaConnectoR::juliaPut(list(values)))
  JuliaConnectoR::juliaPut(values)
}

# Replace NA with the sentinel the engine's column API expects. Done here rather
# than in Julia so every column arrives concretely typed, with no
# `Union{Missing,T}` in the model constructor.
.ctJuliaNoNA <- function(values, sentinel) {
  values[is.na(values)] <- sentinel
  values
}

# Marshal an atomic vector, preserving vector-ness. JuliaConnectoR maps a
# length-one R vector to a Julia scalar, so that case has to go across as a
# list; every other length must not, because a list is marshalled element by
# element (see .ctJuliaNumericVector).
.ctJuliaVector <- function(values) {
  if (!length(values)) {
    stop("Internal error: an empty vector cannot be marshalled to Julia; ",
      "omit the argument instead.", call. = FALSE)
  }
  if (length(values) == 1L) return(JuliaConnectoR::juliaPut(list(values)))
  JuliaConnectoR::juliaPut(values)
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
  # Stan's population covariance is built entirely in raw-parameter units
  # (`rawpopcovbase`/`rawpopsd`, via `sdscale`), then explicitly rescaled to
  # state units: for every indvarying T0MEANS row, `ctModelWriter.R` multiplies
  # T0cov's corresponding row *and* column by that row's `multiplier*meanscale`
  # (see `T0cov[matsetup[ri,1], ] *= matvalues[ri,2] * matvalues[ri,3]` and the
  # matching column update, ctModelWriter.R:903-910). That doubles-up to a
  # `k_i*k_j` scaling of every population-covariance entry, where
  # `k_i = multiplier_i*meanscale_i` comes from state i's own T0MEANS
  # transform. Without it, a state whose T0MEANS uses a non-unit
  # multiplier/meanscale (e.g. the default `10*param` for T0MEANS/CINT-type
  # customs pars) gets a population SD that's wrong by a factor of `k_i`
  # relative to Stan -- exactly the T0VAR population-SD mismatch seen for
  # ctsemTutorial.qmd's individual-differences model. Folding `k_i` into the
  # diagonal (SD) transform here reproduces Stan's row+column rescaling
  # exactly, since sdcovsqrt2cov later builds T0cov[i,j] = sd_i*sd_j*corr_ij:
  # scaling sd_i by k_i and sd_j by k_j automatically scales every entry
  # (including off-diagonals) by k_i*k_j, with no separate adjustment needed
  # for the correlation parameters themselves (they stay dimensionless).
  t0means_state_scale <- rep(1, length(augmented_indices))
  if (!is.null(expanded$modelmats$matsetup)) {
    setup <- as.data.frame(expanded$modelmats$matsetup)
    values <- as.data.frame(expanded$modelmats$matvalues)
    varying_parameters <- unique(setup$param[setup$indvarying > 0L])
    varying_parameters <- varying_parameters[varying_parameters > 0L]
    random_sd_scale <- values$sdscale[match(varying_parameters, setup$param)]
    random_sd_scale[is.na(random_sd_scale)] <- 1
    t0means_code <- ctsem:::ctStanMatricesList()$all[["T0MEANS"]]
    t0means_rows <- setup$matrix == t0means_code & setup$col == 1L
    match_position <- match(augmented_indices, setup$row[t0means_rows])
    t0means_setup_rows <- which(t0means_rows)[match_position]
    t0means_state_scale <- values$multiplier[t0means_setup_rows] * values$meanscale[t0means_setup_rows]
    t0means_state_scale[is.na(t0means_state_scale)] <- 1
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
    table$transform[index] <- sprintf("%.17g * (1e-10 + %.17g * log1p_exp(2 * param[%d] - 1))",
      t0means_state_scale[position], random_sd_scale[position], next_parameter)
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
  # `augmented_indices` (= Stan's `intoverpopindvaryingindex`) is every state
  # with population-varying T0VAR: both the newly-created carrier states for
  # non-T0MEANS random effects (DRIFT/CINT/etc., always appended contiguously
  # after the original states -- see `extralatents` in
  # ctModelWriter.R::ctStanModelIntOverPop) AND any *original* state whose
  # own T0MEANS is directly indvarying (e.g. a random initial value). Only the
  # former are static/no-own-dynamics carriers; the latter are still genuine
  # dynamic states that need their own diffusion/Lyapunov treatment, exactly
  # like Stan's own `derrind` (ctData.R) explicitly excludes only indices
  # beyond `standata$nlatent` ("stable individual differences"), not every
  # index with population-varying T0VAR. `setdiff(seq_len(nlatent_augmented),
  # augmented_indices)` conflated these two groups, silently emptying
  # dynamic_state_indices (and forcing the Lyapunov solve over the *entire*
  # augmented state space, including static carriers with structurally zero
  # diffusion) for any model combining T0MEANS random effects with other
  # (DRIFT/CINT/etc.) random effects -- exactly the combination in
  # ctsemTutorial.qmd's individual-differences example, which failed with a
  # LAPACKException from the Schur-based Lyapunov solver once the augmented
  # dimension exceeded 4. This does not yet replicate Stan's further
  # optimization of also excluding original states with structurally zero,
  # uncoupled diffusion (`derrind`'s first two steps) -- it conservatively
  # includes every original state, which is correct but not maximally
  # reduced.
  list(parameter_table = table, nlatent = original_nlatent,
    nlatent_augmented = nlatent_augmented,
    dynamic_state_indices = seq_len(original_nlatent),
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

.ctJuliaPrepare <- function(datalong, model, prepared_data = NULL, project = NULL,
  priors = FALSE) {
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
  ti_effects <- .ctJuliaTIEffects(parameter_table, model)
  npar <- max(c(parameter_table$parnumber, ti_effects$coefficient), na.rm = TRUE)
  prior_spec <- if (isTRUE(priors)) .ctBackendPriorSpec(prepared_data, npar) else NULL
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
    ti_effects = ti_effects,
    priors = prior_spec,
    max_timestep = max_timestep,
    # A discrete-time model is the same filter with a different discretization:
    # DRIFT, CINT and DIFFUSION are already the one-step quantities, so the
    # exponential, the Lyapunov solve and the intercept solve all collapse.
    continuoustime = isTRUE(model$continuoustime),
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
  # Plain column vectors, not a DataFrame. The engine dropped DataFrames as a
  # dependency (it was its most expensive one and was used only as a row
  # container here), so absent entries arrive as sentinels -- 0 for parnumber,
  # NaN for value, "" for a transform -- rather than as NA/missing. This also
  # removes two RPC round trips per objective build.
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) {
    effects <- data.frame(parameter = integer(), predictor = integer(), coefficient = integer())
  }
  # The optional columns are Julia keyword arguments and are passed only when
  # they have entries: JuliaConnectoR hangs marshalling an empty vector, so a
  # model with no TI predictors must not send one at all.
  arguments <- list(
    .ctJuliaVector(as.character(table$matrix)),
    .ctJuliaVector(as.integer(table$row)),
    .ctJuliaVector(as.integer(table$col)),
    .ctJuliaVector(.ctJuliaNoNA(as.integer(table$parnumber), 0L)),
    .ctJuliaVector(.ctJuliaNoNA(as.numeric(table$value), NaN)),
    .ctJuliaVector(.ctJuliaNoNA(as.character(table$transform), "")),
    .ctJuliaVector(.ctJuliaNoNA(as.character(table$predicttransform), "")),
    .ctJuliaVector(.ctJuliaNoNA(as.character(table$updatetransform), "")),
    .ctJuliaVector(.ctJuliaNoNA(as.character(table$tdtransform), "")))
  if (nrow(effects)) {
    arguments$ti_parameter <- .ctJuliaVector(as.integer(effects$parameter))
    arguments$ti_predictor <- .ctJuliaVector(as.integer(effects$predictor))
    arguments$ti_coefficient <- .ctJuliaVector(as.integer(effects$coefficient))
  }
  if (length(spec$dynamic_state_indices)) {
    arguments$diffusion_state_indices <-
      .ctJuliaVector(as.integer(spec$dynamic_state_indices))
  }
  arguments$continuous_time <- isTRUE(spec$continuoustime)
  params <- do.call(module$ekf_from_columns, arguments)
  # .ctJuliaVector, not juliaPut, for the vectors: JuliaConnectoR marshals a
  # length-one R vector as a *scalar*, so a single-subject model would hand the
  # objective constructor an Int where it wants an AbstractVector. That never
  # showed up while every caller was a whole fitted dataset; ctPredict() on one
  # subject is a caller where it does.
  objective_args <- list(params, .ctJuliaVector(spec$subject_starts),
    .ctJuliaVector(spec$times), JuliaConnectoR::juliaPut(spec$manifest_data),
    JuliaConnectoR::juliaPut(spec$tdpred_data), JuliaConnectoR::juliaPut(spec$tipred_data),
    spec$max_timestep)
  if (!is.null(spec$priors) && length(spec$priors$index)) {
    objective_args$prior_index <- .ctJuliaVector(spec$priors$index)
    objective_args$prior_scale <- .ctJuliaVector(spec$priors$scale)
    objective_args$prior_weight <- spec$priors$weight
  }
  objective <- do.call(module$ctsem_objective, objective_args)
  assign(key, objective, envir = .ct_julia_cache$objectives)
  objective
}

#' Evaluate a prepared Julia ctsem likelihood
#' @param object A \code{ctJuliaModel} or \code{ctJuliaFit}.
#' @param pars Unconstrained parameters; defaults to the fitted estimate.
#' @param gradient Return an analytic/automatic-differentiation gradient.
#' @param contributions Return subject and row contribution details when available.
#' @param gradient_method Either "adjoint" (reverse mode, the default) or
#'   "forward" (ForwardDiff). Both compute the same gradient. "adjoint" costs
#'   the same regardless of how many free parameters a model has, and is
#'   faster than "forward" at every model size measured, by a margin that
#'   grows with the parameter count. See \code{optimcontrol$gradient} in
#'   \code{ctFit}.
#' @return A list containing log likelihood and, when requested, gradient.
#' @export
ctJuliaEvaluate <- function(object, pars = NULL, gradient = TRUE, contributions = FALSE,
  gradient_method = c("adjoint", "forward")) {
  gradient_method <- match.arg(gradient_method)
  if (!inherits(object, c("ctJuliaModel", "ctJuliaFit"))) stop("object must be a ctJuliaModel or ctJuliaFit", call. = FALSE)
  if (is.null(pars)) {
    if (inherits(object, "ctJuliaFit")) pars <- object$estimate$raw else stop("pars must be supplied for a prepared ctJuliaModel", call. = FALSE)
  }
  module <- .ctJuliaModule(if (inherits(object, "ctJuliaFit")) object$model_spec$project else object$project)
  result <- module$ctsem_evaluate(.ctJuliaObjective(object), .ctJuliaNumericVector(pars),
    gradient = isTRUE(gradient), contributions = isTRUE(contributions),
    gradient_method = gradient_method)
  JuliaConnectoR::juliaGet(result)
}

#' @export
summary.ctJuliaFit <- function(object, timeinterval = 1, digits = 3, parmatrices = TRUE, ...) {
  .ctBackendSummary(object, timeinterval = timeinterval, digits = digits,
    parmatrices = parmatrices, ...)
}

#' @export
ctExtract.ctJuliaFit <- function(object, subjectMatrices = FALSE, cores = 2,
  nsamples = "all", subjects = "all", ...) {
  .ctBackendExtract(object, subjectMatrices = subjectMatrices, nsamples = nsamples,
    subjects = subjects, ...)
}

#' @export
ctSummaryMatrices.ctJuliaFit <- function(fit, calcfunc = quantile,
  calcfuncargs = list(probs = 0.5), timeinterval = 1, ...) {
  .ctBackendSummaryMatrices(fit, calcfunc = calcfunc, calcfuncargs = calcfuncargs,
    timeinterval = timeinterval, ...)
}

ctFitJuliaBackend <- function(datalong, model, prepared_data = NULL, inits = NULL, cores = 1L,
  backendcontrol = list(), optimcontrol = list(), verbose = 0L, fit = TRUE,
  priors = FALSE) {
  if (isTRUE(backendcontrol$restart_session)) .ctJuliaClearSession()
  project <- .ctJuliaOr(backendcontrol$julia_project, NULL)
  # `cores` splits the engine's subject loop. It is requested as a Julia thread
  # count before the session starts (which is the only time that can be set),
  # and capped per fit afterwards, so a session started with more threads is not
  # forced to use them all.
  cores <- max(1L, suppressWarnings(as.integer(cores)[1L]))
  if (is.na(cores)) cores <- 1L
  if (cores > 1L && !.ctJuliaSessionRunning() &&
      !nzchar(Sys.getenv("JULIA_NUM_THREADS", unset = ""))) {
    Sys.setenv(JULIA_NUM_THREADS = as.character(cores))
  }
  # `optimcontrol$gradient` is the documented control; `backendcontrol$gradient`
  # is still honoured because it was the only way to set this before, and
  # silently ignoring it would change results for anyone already passing it.
  gradient <- .ctJuliaOr(optimcontrol$gradient, .ctJuliaOr(backendcontrol$gradient, "adjoint"))
  if (!gradient %in% c("forward", "adjoint")) stop("gradient must be 'forward' or 'adjoint'", call. = FALSE)
  # 'adjoint' selects the Julia engine's reverse-mode gradient. Its cost is
  # independent of the free-parameter count (one traced forward sweep plus one
  # reverse sweep per subject), where 'forward' (ForwardDiff) costs one dual
  # pass per chunk of parameters. Measured on one gradient evaluation, 20
  # subjects x 5 waves, default ctsem parameterisation:
  #
  #   latents  free pars   forward     adjoint
  #        2         23    0.009 s     0.006 s
  #        6        153    0.062 s     0.019 s
  #       20       1490   39.1   s     0.168 s
  #
  # 'adjoint' is faster at every size measured, and the margin grows without
  # bound with the parameter count. 'forward' remains the default because it
  # is the longer-tested path, not because it is faster; there is no silent
  # fallback between them in either direction.
  model_spec <- .ctJuliaPrepare(datalong, model, prepared_data = prepared_data,
    project = project, priors = priors)
  if (!fit) return(structure(model_spec, class = c("ctJuliaModel", "ctFitModel")))

  objective <- .ctJuliaObjective(structure(model_spec, class = c("ctJuliaModel", "ctFitModel")))
  npar <- max(c(model_spec$parameter_table$parnumber, model_spec$ti_effects$coefficient), na.rm = TRUE)
  start <- .ctJuliaInitialValues(npar, inits)
  module <- .ctJuliaModule(project)
  # Called by name rather than through the imported module: the Julia function
  # ends in `!`, which is not a syntactic R name.
  JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!", cores)
  result <- JuliaConnectoR::juliaGet(module$ctsem_optimize(objective, .ctJuliaNumericVector(start),
    maxiter = as.integer(.ctJuliaOr(backendcontrol$maxiter, 1000L)),
    g_tol = .ctJuliaOr(backendcontrol$g_tol, 1e-8),
    f_tol = .ctJuliaOr(backendcontrol$f_tol, 0),
    x_tol = .ctJuliaOr(backendcontrol$x_tol, 0),
    verbose = verbose > 0L,
    gradient_method = gradient))
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
