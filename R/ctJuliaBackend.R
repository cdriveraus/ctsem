# Julia likelihood backend -------------------------------------------------
#
# The Julia process is deliberately process-local.  Fit objects retain only a
# serializable specification and rebuild their proxy on demand after saveRDS.

.ct_julia_cache <- new.env(parent = emptyenv())
.ct_julia_cache$objectives <- new.env(parent = emptyenv())
# Keyed with the same objective key: a layout describes the objective's model,
# so it is stale in exactly the cases the objective is.
.ct_julia_cache$layouts <- new.env(parent = emptyenv())
.ctJuliaOr <- function(x, default) if (is.null(x)) default else x

# Can a carriage return move the cursor here?
#
# The progress reporter overwrites one line in place, which requires that `\r`
# reaches something that treats it as a cursor movement. A console does; a log
# file, a knitr chunk and a captured stream do not, and there the same updates
# have to be rarer and on their own lines or the output becomes one long line of
# accumulated garbage.
#
# `interactive()` was the test and is the wrong one: it asks whether a human is
# at a prompt, not whether output reaches a console. Under Shiny it is TRUE
# while stdout is being captured for a log pane -- the case that motivated
# replacing it -- and it is FALSE in a plain `R -f` run whose output is going
# straight to a terminal that handles `\r` perfectly well.
#
# So: detect the cases that are actually known to capture, and let
# `options(ctsem.progress.overwrite = )` settle it outright for anything this
# does not know about. A front end that captures output should set it to FALSE
# once at startup rather than passing an argument through every fitting call.
#' @keywords internal
.ctProgressConsole <- function() {
  # `options(ctsem.progress.overwrite = )` is *not* read here, and used to be.
  # This function answers "is anyone watching", and its answer decides whether
  # a fit reports at all; the option says whether a carriage return moves the
  # cursor, which `.ctProgressOverwrite()` reads. Conflating them meant that
  # setting the option FALSE -- to ask for plain lines rather than an
  # overwritten one -- silently turned reporting off altogether, which is not
  # what it says and not what anyone setting it wants.
  # Inside a Shiny session: a non-NULL reactive domain is the reliable signal,
  # and reaching it through the namespace keeps shiny a suggestion rather than
  # a dependency.
  if ("shiny" %in% loadedNamespaces()) {
    domain <- try(get("getDefaultReactiveDomain",
      envir = asNamespace("shiny"))(), silent = TRUE)
    if (!inherits(domain, "try-error") && !is.null(domain)) return(FALSE)
  }
  # sink() and capture.output() divert stdout to a connection; knitr collects
  # chunk output. Both keep the carriage return as a character.
  if (sink.number() > 0L) return(FALSE)
  if (isTRUE(getOption("knitr.in.progress"))) return(FALSE)
  interactive()
}

# Whether to overwrite, which is the console question and the history question
# together. `verbose >= 2` asks to keep every update, and at that point the
# history is the reason it was turned on.
#
# One line that is rewritten is the whole design, and it is not given up on for
# a console that might not honour a carriage return: the alternative is a stream
# of updates down the screen, which is what this exists to avoid. A slower
# cadence only makes that stream shorter. So overwriting stays on wherever
# anything is being watched, and the option below settles the cases this cannot
# know about -- in both directions, which it previously did not.
#
# The RStudio case is why that matters. Its console renders each arriving chunk
# as its own block, and while the engine printed to *stdout* while R messaged on
# *stderr*, the two interleaved and the line could not survive: reported as
# "extra line breaks". Both now arrive as messages on one stream, which removes
# that cause -- so in-place is the default there too, and
# `options(ctsem.progress.overwrite = FALSE)` is the way out if a particular
# console still cannot do it.
#' @keywords internal
.ctProgressOverwrite <- function(verbose = 0) {
  # `verbose` is a level on the fitting paths and a flag on `ctSample()`; a
  # flag means level one, which still overwrites.
  level <- if (is.numeric(verbose) && length(verbose) == 1L && !is.na(verbose)) {
    verbose
  } else if (isTRUE(verbose)) 1 else 0
  option <- getOption("ctsem.progress.overwrite")
  if (is.logical(option) && length(option) == 1L && !is.na(option)) {
    return(option && level < 2)
  }
  .ctProgressConsole() && level < 2
}

#' Deliver an engine progress line as an R message
#'
#' Returns a function of `(text, kind)` -- the protocol `_deliver!()` in
#' `progress.jl` calls, with `kind` one of `"update"`, `"done"` or `"break"`.
#'
#' It exists because a console draws a *message* differently from text that
#' merely arrives on the same stream. In RStudio, `message()` output is a
#' shaded block and a raw write to `stderr()` is plain body text, so an engine
#' that printed its own line -- however correct the stream -- produced a fit
#' whose progress did not look like anything else that fit said. The engine
#' still owns the line's content: only its delivery is here, so there is no
#' second spelling of the fields in a second language.
#'
#' The carriage return and the padding are this side's business too, for the
#' same reason: whether a cursor can be moved at all is a property of the
#' console, not of the engine.
#'
#' Wrapped in `tryCatch` throughout. This is called from inside a running Julia
#' call, and an error thrown back across that boundary does not merely lose a
#' progress line -- it desynchronises the bridge, which costs the fit. The
#' engine drops the sink on its first failure for the same reason.
#'
#' @param overwrite Update one line in place rather than emitting each report
#'   on its own line.
#' @return A function of `(text, kind)`, invisibly stateful.
#' @keywords internal
.ctProgressSink <- function(overwrite = .ctProgressOverwrite(1)) {
  pad <- 0L
  open <- FALSE
  function(text, kind = "update") {
    tryCatch({
      text <- as.character(text)[1L]
      if (is.na(text)) text <- ""
      if (identical(kind, "break")) {
        # Only ends the line, so that something else can print on its own.
        if (open) {
          message("")
          open <<- FALSE
          pad <<- 0L
        }
        return(invisible(NULL))
      }
      final <- identical(kind, "done")
      if (!overwrite) {
        message(text)
        open <<- FALSE
        return(invisible(NULL))
      }
      pad <<- max(pad, nchar(text))
      # `appendLF = FALSE` for an update, so the next carriage return lands on
      # this line; the closing line ends itself. The padding is what stops a
      # shorter update leaving the tail of a longer one behind it.
      message("\r", formatC(text, width = -pad), appendLF = final)
      open <<- !final
      if (final) pad <<- 0L
      invisible(NULL)
    }, error = function(e) invisible(NULL))
  }
}

# Is progress reporting asked for at all -- `verbose` read the way it is
# everywhere else, a level as well as a flag.
#
# `isTRUE(verbose)` is the wrong test once `verbose` can be `2`: `isTRUE(2)` is
# `FALSE`, so a caller who raised the level for more detail would silently get
# no progress output at all. This is the same reading `.ctProgressOverwrite`
# already gives `verbose`, pulled out because it is also needed everywhere a
# default of `progress = isTRUE(verbose)` has to become `progress =
# .ctVerboseOn(verbose)` instead.
#' @keywords internal
.ctVerboseOn <- function(verbose) {
  isTRUE(verbose) ||
    (is.numeric(verbose) && length(verbose) == 1L && !is.na(verbose) && verbose > 0)
}

# `1m 04s`, `12.4s`, `2h 05m` -- the R-side twin of `_duration()` in
# `progress.jl`. Needed on this side too because `.ctBackendReportProcesses`
# formats an ETA from what a worker process wrote to a file, not from
# anything the engine printed, so there is no Julia-side line to reuse.
#' @keywords internal
.ctDuration <- function(seconds) {
  seconds <- max(0, as.numeric(seconds)[1L])
  if (!is.finite(seconds)) return("--")
  if (seconds < 60) return(sprintf("%.1fs", seconds))
  if (seconds < 3600) {
    return(sprintf("%dm %02ds", as.integer(seconds %/% 60),
      as.integer(round(seconds %% 60))))
  }
  sprintf("%dh %02dm", as.integer(seconds %/% 3600),
    as.integer(round((seconds %% 3600) / 60)))
}
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

# The engine's identity: a hash of its own source.
#
# This used to be a git commit recorded in `inst/julia/engine.json` by a script
# that copied the engine in from a separate repository. The engine is now simply
# part of ctsem -- `inst/julia/ContinuousTimeSEM/` is the source, edited in place
# like any other file here -- so there is no upstream commit to record, and no
# bookkeeping step that can be forgotten.
#
# Hashing the content rather than maintaining a version is what makes the cached
# project self-invalidating. That cache is keyed on this and is only populated
# when empty, so any identifier that has to be *updated by hand* strands a user
# on a stale engine the moment someone edits the source without bumping it: they
# upgrade ctsem, get the new R code, and run it against the old engine, silently.
# A content hash cannot be forgotten.
# Deliberately not memoized. Hashing 45 files takes 5.6 ms against 1.1 ms for a
# cached lookup, and the saving is invisible next to anything that asks for it
# -- but a cache keyed on the path alone would report a stale version after an
# edit within the same session, which is precisely the property this exists to
# provide. A guarantee with an exception is not a guarantee.
.ctJuliaEngineVersion <- function(path = .ctJuliaEnginePath()) {
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE))
  manifest <- paste(substring(files, nchar(path) + 2L), tools::md5sum(files),
    collapse = "
")
  handle <- tempfile(fileext = ".txt")
  on.exit(unlink(handle), add = TRUE)
  writeLines(manifest, handle)
  substr(unname(tools::md5sum(handle)), 1L, 12L)
}

# Path to the copy of ContinuousTimeSEM.jl shipped inside this ctsem install.
.ctJuliaEnginePath <- function() {
  path <- system.file("julia", "ContinuousTimeSEM", package = "ctsem")
  if (!nzchar(path) || !file.exists(file.path(path, "Project.toml"))) {
    stop("The Julia engine is missing from this ctsem installation; reinstall ctsem.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# A writable project directory for the engine, keyed by the engine's content
# hash so that any change to it -- a ctsem upgrade, or an edit made while
# working on the engine -- gets a fresh environment instead of reusing a stale
# one.
#
# The engine tree is *copied* here rather than activated in place: activating a
# project writes to its Manifest.toml, and an R library directory is frequently
# read-only. The copy is under 400 KB.
.ctJuliaEnvDir <- function(version = .ctJuliaEngineVersion()) {
  base <- tools::R_user_dir("ctsem", which = "cache")
  # Forward slashes throughout: this path is used by R and also embedded in
  # Julia source, and Julia accepts them on every platform.
  gsub("\\", "/", file.path(base, "julia", paste0("engine-", version)), fixed = TRUE)
}

# Is a Julia session live? *Without* starting one.
#
# This used to fall back to `JuliaConnectoR::juliaEval("true")`, which starts a
# session rather than reporting on one -- so the predicate created the situation
# it was asked about. Two things followed, both silent:
#
#   * `ctJuliaSetup(threads = 10, force = TRUE)` cleared the session, then asked
#     this question, then warned that a session was already running and refused
#     to apply `threads`. `force` could not work.
#   * `ctFit(backend = 'julia', cores = n)` sets `JULIA_NUM_THREADS` only when no
#     session is running. The check started one first, so the variable was never
#     set and every fit ran single-threaded whatever `cores` said.
#
# JuliaConnectoR exports no non-starting predicate, so its connection is read
# directly and defensively: if that internal ever moves, the answer degrades to
# "not running", which is the safe direction -- we then set the thread count,
# which is a no-op when it was already right.
.ctJuliaSessionRunning <- function() {
  if (!is.null(.ct_julia_cache$module)) return(TRUE)
  isTRUE(tryCatch({
    connection <- get("pkgLocal", envir = asNamespace("JuliaConnectoR"))$con
    !is.null(connection) && isOpen(connection)
  }, error = function(e) FALSE))
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
#' installation. The engine source itself is vendored in \code{inst/julia/} --
#' no repository credentials and no network access are involved in getting
#' that -- but its Julia package dependencies still have to be instantiated,
#' which reaches the network on a machine where they are not already cached.
#' That step is consent-gated exactly as in \code{\link{ctJuliaInstall}}: asked
#' for interactively, skipped with a message when the environment is already
#' instantiated, and declined by default in a non-interactive session unless
#' \code{agree = TRUE} or \code{CTSEM_JULIA_AGREE=yes} said otherwise in
#' advance.
#'
#' The first call downloads and precompiles those dependencies (roughly 120 MB
#' and a minute or two); later calls in new sessions reuse them.
#'
#' If \pkg{JuliaConnectoR} or Julia itself is missing, this offers to install it
#' rather than failing -- the same thing \code{\link{ctJuliaInstall}} does, which
#' is the function to reach for when setting the backend up deliberately, or
#' from a script.
#' @param project Optional local ContinuousTimeSEM.jl checkout to use instead of
#'   the copy that ships with ctsem. Rarely needed: the engine is part of ctsem
#'   and is edited in place under \code{inst/julia/}.
#' @param revision Ignored; retained for backward compatibility. The engine is
#'   part of ctsem, and the version \code{ctJuliaStatus()} reports is a hash of
#'   its source rather than something selectable.
#' @param julia_bin Optional Julia binary directory.
#' @param threads Number of Julia threads, fixed for the life of the session.
#'   The engine splits its subject loop across them, and \code{ctFit}'s
#'   \code{cores} cannot exceed this count: Julia fixes it at process start, so
#'   setting it here takes effect only if no Julia session is running yet --
#'   pass \code{force = TRUE} to restart one. \code{NULL} leaves it to Julia's
#'   own default (one thread unless \code{JULIA_NUM_THREADS} is already set),
#'   which is why starting the engine before the first fit and then asking for
#'   \code{cores = 4} gives four chunks' worth of nothing. A fit that asks for
#'   more cores than the session has threads says so; see \code{\link{ctFit}}.
#' @param force Reconfigure an existing Julia session.
#' @param agree \code{TRUE} to consent to instantiating the engine's Julia
#'   package dependencies without being asked, \code{FALSE} to refuse.
#'   \code{NULL}, the default, asks in an interactive session and refuses
#'   otherwise; see \code{\link{ctJuliaInstall}}.
#' @return A Julia-engine status list, invisibly.
#' @seealso \code{\link{ctJuliaInstall}}, \code{\link{ctJuliaStatus}}
#' @export
ctJuliaSetup <- function(project = NULL, revision = "locked", julia_bin = NULL,
  threads = NULL, force = FALSE, agree = NULL) {
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
  engineversion <- .ctJuliaEngineVersion()

  if (!is.null(project)) {
    # Developer override: use the checkout as its own project, in place.
    project <- normalizePath(project, winslash = "/", mustWork = TRUE)
    env_dir <- project
  } else {
    env_dir <- .ctJuliaEnvDir(engineversion)
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

  JuliaConnectoR::juliaEval("using Pkg, Logging")
  # Quietly. `Pkg.activate()` announces itself, and the environment is ctsem's
  # own vendored one -- the first thing a user saw was
  # `Activating project at C:\Users\...\engine-4f87c9793c9f`, naming a cache
  # path, followed by a three-line Pkg warning recommending `Pkg.resolve()` on a
  # manifest they did not write. Alarming, and about a situation this code
  # already handles: the resolve fallback below is exactly that recommendation,
  # taken automatically.
  #
  # `io=devnull` silences the announcements and a NullLogger the warnings.
  # Failures are exceptions rather than logs, so they still propagate to the
  # fallback and to the user.
  activate <- sprintf("Pkg.activate(%s; io=devnull)", .ctJuliaString(env_dir))
  silently <- function(code) paste0(
    "Logging.with_logger(Logging.NullLogger()) do; ", code, "; end")
  JuliaConnectoR::juliaEval(silently(activate))

  # `Pkg.instantiate()` is the network-touching step: it fetches whatever the
  # manifest lists that is not already in the user's Julia depot. Consent-gated
  # exactly like the other two installs `ctJuliaInstall()` offers (see
  # `.ctJuliaAgreed()`), and skipped entirely -- no prompt, no network -- when
  # the environment already loads, which is the ordinary case after the first
  # session on a machine.
  ready <- tryCatch({
    JuliaConnectoR::juliaEval("using ContinuousTimeSEM")
    TRUE
  }, error = function(e) FALSE)

  if (!ready) {
    agreed <- .ctJuliaAgreed(agree, paste0(
      "The julia backend needs to install the vendored engine's Julia package ",
      "dependencies, which are not yet present.\n",
      "  environment: ", env_dir, "\n",
      "  This downloads packages from the Julia package registry."))
    if (!agreed) {
      stop("The julia backend needs its Julia package dependencies installed, ",
        "and consent was not given.\n", .ctJuliaDeclined(), call. = FALSE)
    }

    # `Pkg.instantiate()` precompiles the environment, and then `using` precompiles
    # the engine again -- two full builds of the same package on every fresh
    # session with a cold cache. That was cheap until the engine started
    # precompiling model shapes; now it is a hundred seconds paid twice.
    # `JULIA_PKG_PRECOMPILE_AUTO` turns off only Pkg's automatic pass, so `using`
    # still builds whatever is stale, and `withenv` puts it back rather than
    # leaving the session's Pkg quietly reconfigured.
    quiet_instantiate <- paste0('withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do; ',
      'Pkg.instantiate(io=devnull); end')
    instantiated <- tryCatch({
      JuliaConnectoR::juliaEval(silently(paste0(activate, "; ", quiet_instantiate)))
      TRUE
    }, error = function(e) FALSE)
    if (!instantiated) {
      # The vendored manifest pins the versions this ctsem release was tested
      # against, but it can be unsatisfiable on a different Julia version. Falling
      # back to a fresh resolve is better than refusing to run; the compat bounds
      # in Project.toml still apply.
      unlink(file.path(env_dir, "Manifest.toml"))
      JuliaConnectoR::juliaEval(silently(paste0(activate,
        "; Pkg.resolve(io=devnull); ", quiet_instantiate)))
    }
    JuliaConnectoR::juliaEval("using ContinuousTimeSEM")
  }
  .ctJuliaTuneBridge()
  .ct_julia_cache$project <- project
  .ct_julia_cache$engine <- engineversion
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
#'   the \code{engine} version (a hash of the engine source shipped with this
#'   ctsem, which is also what its cached project directory is keyed on), and
#'   the number of \code{threads} in a running session.
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
  threads <- if (available) {
    tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
      error = function(e) NA_integer_)
  } else NA_integer_
  list(available = available, connectoR = connectoR,
    julia_bin = .ctJuliaOr(julia_bin, NA_character_), julia = version,
    project = .ctJuliaOr(project, .ct_julia_cache$project),
    engine = .ctJuliaOr(.ct_julia_cache$engine,
      tryCatch(.ctJuliaEngineVersion(), error = function(e) NA_character_)),
    threads = threads)
}

# The socket the bridge runs on ----------------------------------------------
#
# JuliaConnectoR sends every message as a run of small writes -- an indicator
# byte, a length, a name, a count, then each argument -- and reads the reply the
# same way. On Linux that write-write-read pattern meets Nagle's algorithm on
# the sending side and a delayed ACK on the receiving side, and the second write
# waits up to 40 ms for an acknowledgement the kernel is deliberately holding
# back. Measured on a bare socket with neither Julia nor R packages in the way
# (23-core Linux box): one write per request 0.033 ms, two writes 40.7 ms, eight
# writes 40.8 ms. One stall per message, independent of payload.
#
# That is the whole of what looked like a per-round-trip latency. A bridge
# message cost 41 ms when one direction split into several writes and 82 ms when
# both did, which made `ctJuliaEvaluate` ~370 ms and read as "about nine round
# trips of 42 ms". It sends six messages, not nine, and the unit is a TCP stall
# rather than a round trip. (Six: marshal the parameter vector, release the
# proxies R collected since the last call, evaluate, then three for `juliaGet`,
# which sets a flag, fetches, and unsets it.)
#
# `ctsem_tune_bridge!` in the engine turns both halves off. Same box, per
# message: 82 ms untuned, 41 ms with Nagle disabled, 2.5 ms with the quickack
# task as well. Windows and macOS are unaffected either way -- Windows measures
# ~0.3 ms a message untuned, and TCP_QUICKACK does not exist there.
#
# What is left, quantified rather than guessed at: three of the six messages are
# `juliaGet`, which exists only because the engine returns a NamedTuple and
# JuliaConnectoR sends those by reference. A flat `Array{Float64}` comes back by
# value in the call's own message, so a flattened R-facing return would drop
# those three. That was worth ~220 ms a call before the socket fix and is worth
# ~5 ms after it, which is why it has not been done -- and why it should not be
# done for speed alone. The cleaner fix belongs upstream: if JuliaConnectoR
# buffered each message into one `send()`, neither the stall nor the task below
# would exist.
#
# Everything here is advisory: a failure costs the tuning, never the session.
.ctJuliaCommunicator <- function() {
  # The Julia object that owns the socket. JuliaConnectoR keeps it in a
  # package-private environment and passes it to Julia itself when a call needs
  # it -- `juliaGet` does exactly this. `asNamespace` rather than `:::` because
  # the latter is an R CMD check note, and the tryCatch is for the day the name
  # changes.
  tryCatch({
    local_env <- get("pkgLocal", envir = asNamespace("JuliaConnectoR"))
    communicator <- local_env$communicator
    if (inherits(communicator, "JuliaProxy")) communicator else NULL
  }, error = function(e) NULL)
}

.ctJuliaTuneBridge <- function() {
  # `ctsem.julia.tunebridge = FALSE` leaves the socket exactly as JuliaConnectoR
  # opened it. It exists as an escape hatch and as the control arm for measuring
  # what the tuning is worth: the options are set on the connection, so an A/B
  # needs two sessions rather than a switch inside one.
  if (!isTRUE(getOption("ctsem.julia.tunebridge", TRUE))) return(invisible(NA_integer_))
  # `ctsem.julia.quickack` is how often the Julia task re-arms TCP_QUICKACK, in
  # seconds; zero or less asks for no task at all, which leaves the Nagle half of
  # the fix in place and the R->Julia half not. Floored at 0.1 ms rather than
  # taken literally: the task's `sleep` becomes a spin below that, and the stall
  # it is removing is 40 ms, so there is nothing to win by asking for less.
  interval <- getOption("ctsem.julia.quickack", 0.001)
  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval)) {
    interval <- 0.001
  } else if (interval > 0) {
    interval <- max(interval, 1e-4)
  }
  communicator <- .ctJuliaCommunicator()
  if (is.null(communicator)) return(invisible(NA_integer_))
  out <- tryCatch(as.integer(.ctBackendJuliaValue(JuliaConnectoR::juliaCall(
    "ContinuousTimeSEM.ctsem_tune_bridge!", communicator,
    quickack_interval = as.numeric(interval)))),
    error = function(e) NA_integer_)
  invisible(out)
}

.ctJuliaModule <- function(project = NULL) {
  if (!is.null(.ct_julia_cache$module) && identical(project, .ct_julia_cache$project)) return(.ct_julia_cache$module)
  ctJuliaSetup(project = project)
  .ct_julia_cache$module
}

# Drop every Julia proxy *before* stopping Julia, and collect them while the
# session that issued them is still there to hear it.
#
# JuliaConnectoR hands R proxy objects that hold a reference counted on the
# Julia side, and releases them from an R finalizer. Stopping Julia first and
# clearing the cache afterwards leaves those proxies alive with nothing to talk
# to: their finalizers run whenever R next collects, which may be after a new
# session has started, and the decrement then names a reference that session
# never issued. That surfaces as
#
#     KeyError: key 0x0000... not found
#       decrefcount!(communicator, ref) at sharing.jl:134
#
# from inside an unrelated later call. Order is the whole fix -- release, then
# collect, then stop.
.ctJuliaClearSession <- function() {
  .ct_julia_cache$module <- NULL
  .ct_julia_cache$project <- NULL
  .ct_julia_cache$engine <- NULL
  .ct_julia_cache$objectives <- new.env(parent = emptyenv())
  .ct_julia_cache$layouts <- new.env(parent = emptyenv())
  # Two passes: the first frees the proxies, the second any proxy a finalizer
  # from the first pass happened to drop.
  gc(verbose = FALSE); gc(verbose = FALSE)
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    try(JuliaConnectoR::stopJulia(), silent = TRUE)
  }
  invisible(NULL)
}

# Julia's bridge is a request/response socket, and an interrupt does not reach
# it. Pressing Escape between R sending a request and reading the reply leaves
# that reply sitting unread: the connection is not merely interrupted, it is
# *desynchronised*, and every later call reads the answer to the call before it.
# Reported from a real session as "can't do any more fits" -- the whole R
# session is unusable for ctsem until it is restarted.
#
# Dropping the connection here costs the Julia session, and with it the engine's
# compiled shapes, so the next fit pays a recompile. That is a far better trade
# than an R session that silently returns the wrong answers, or refuses to fit
# at all, until someone thinks to restart it.
#
# The worker pool has to go the same way, and for the same reason one layer up.
# A sampled fit runs its chains in `future::multisession` workers, and Escape
# during one leaves two kinds of wreckage: the chains keep running, so every
# worker is still busy when the next fit asks for one, and if the interrupt
# landed while the parent was reading a worker's result then that socket is
# desynchronised exactly as the Julia bridge is. Reported from a real session as
# workers that "seem broken", after which sampling silently fell back to one
# process for the rest of the session -- which is what the fallback is *for*,
# but nobody asked for it and nothing said the pool was the reason.
#
# Dropping the pool costs each worker's compiled shapes, as dropping the Julia
# session costs the parent's. Same trade, and it is not close: a recompile is
# tens of seconds and an unusable pool lasts as long as the session.
#
# `withCallingHandlers` rather than `tryCatch`: the handler runs and the
# interrupt then carries on unwinding, so Escape still aborts the fit. It only
# stops leaving wreckage behind.
#' @keywords internal
.ctJuliaInterruptSafe <- function(expr) {
  withCallingHandlers(expr, interrupt = function(cnd) {
    message("Interrupted. Restarting the Julia session and any sampling ",
      "workers, because a half-finished request would leave every later fit ",
      "in this session reading the wrong reply. The next fit will recompile ",
      "the engine for its model shape.")
    try(.ctBackendWarmStop(NULL), silent = TRUE)
    try(.ctJuliaClearSession(), silent = TRUE)
  })
}

.ctJuliaObjectiveKey <- function(spec) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Julia objective caching requires the suggested package digest.", call. = FALSE)
  }
  digest::digest(list(spec$parameter_table, spec$subject_starts, spec$times,
    spec$manifest_data, spec$tdpred_data, spec$tipred_data,
    spec$ti_effects, spec$priors, spec$max_timestep, spec$project, spec$engine,
    # Two fits differing only in how random effects are integrated share every
    # field above and are not the same objective. Nor do two differing only in
    # the covariance construction: without `covmatcode` here, a 'z' model built
    # after a 'rawcorr' one was handed the earlier objective and silently
    # returned the earlier likelihood.
    spec$intoverpop, spec$laplace, spec$covmatcode),
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

.ctJuliaInitialValues <- function(npar, inits = NULL, initsd = .01) {
  # Match stanoptimis(): absent initial values are small, R-seeded draws in
  # unconstrained space rather than an exact all-zero vector, and `initsd` is
  # its name for the scale, with its default.
  if (is.null(inits) || identical(inits, "random")) {
    return(stats::rnorm(npar, 0, as.numeric(initsd)[1L]))
  }
  values <- as.numeric(inits)
  if (length(values) != npar) {
    stop("Julia initial values must have one entry per free parameter.", call. = FALSE)
  }
  if (anyNA(values)) stop("Julia initial values must be numeric.", call. = FALSE)
  values
}

.ctJuliaUnsupported <- function(model, optimize, priors, intoverpop, gendata,
  stanmodeltext, compileArgs, forcerecompile, intoverstates = TRUE,
  optimcontrol = list()) {
  failures <- character()
  # `optimcontrol$is` selects Stan's optimization-plus-importance-sampling
  # route (see R/ctFit.R's estimation-method message); the julia path's only
  # uncertainty methods are the Hessian and ctOptimUncertainty()'s sampling,
  # neither of which is importance sampling, so silently ignoring it would
  # leave the user told one thing happened when another did.
  if (isTRUE(optimcontrol$is)) {
    failures <- c(failures, "optimcontrol$is (importance sampling)")
  }
  # `intoverstates=FALSE` fits over the joint density of parameters and states
  # (see `state_sampling.jl`), which composes with the default `intoverpop` --
  # augmented random effects are extra states, so they are sampled along with
  # every other one -- and does not compose with the Laplace route, which wraps
  # the likelihood in an inner problem the state path replaces.
  if (!isTRUE(intoverstates) && identical(as.character(intoverpop)[1L], "laplace")) {
    failures <- c(failures,
      "intoverstates=FALSE together with intoverpop='laplace'")
  }
  # `optimize=FALSE` is supported now: the engine has its own No-U-Turn sampler,
  # and which target it samples is decided by `intoverpop`. See
  # `.ctJuliaSampleFit`.
  # Binary and ordinal manifest variables are supported now: the filter
  # integrates the observation rather than linearising it, which is why it is
  # worth having here at all. See
  # inst/julia/ContinuousTimeSEM/src/binary_measurement.jl.
  if (any(!model$manifesttype %in% 0:4)) {
    failures <- c(failures, "manifest types beyond censored (manifesttype > 4)")
  }
  # `covmattransform` reaches the engine nowhere: every call to
  # `sdcovsqrt2cov!` passes a literal 0, and the branch the other values would
  # select is commented out. R computes `standata$choleskymats` for both
  # backends one dispatch above this, so the setting was validated and then
  # ignored, which reads as support -- a 'cholesky' model measured 3.76 log
  # units away from the same model on stan, with no error and no warning, and
  # its summary matrices would be wrong too. Refused for every non-default
  # value rather than only 'cholesky': 'rawcorr_indep' happens to be inert on
  # the stan path as well, so nothing a user can rely on is removed, and the
  # julia contract does not inherit that accident.
  if (!is.null(model$covmattransform) &&
      !as.character(model$covmattransform)[1L] %in% c("rawcorr", "z")) {
    failures <- c(failures, paste0("covmattransform='",
      as.character(model$covmattransform)[1L],
      "' (only 'rawcorr' and 'z')"))
  }
  # A TI effect fixed to a value ('TI1=4.3') is honoured by generation and not
  # by fitting: the coefficient it occupies is an ordinary free parameter to
  # the optimiser, so a fit would estimate it and silently ignore the value.
  fixedeffects <- .ctTipredFixedEffects(model)
  if (length(fixedeffects)) {
    failures <- c(failures, paste0("time independent predictor effects fixed ",
      "to a value (", paste(utils::head(fixedeffects, 4), collapse = ", "),
      ") -- these are for generation; use TRUE for an effect to estimate"))
  }
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
  ctm <- ctModel0DRIFT(model, model$continuoustime)
  ctm$pars <- ctModelStatesAndPARS(
    ctm$pars, statenames = ctm$latentNames, tdprednames = ctm$TDpredNames
  )
  jacobian <- try(ctJacobian(ctm), silent = TRUE)
  if (inherits(jacobian, "try-error")) jacobian <- ctJacobian(ctm, simplify = FALSE)
  jacobian <- jacobian[names(.ctMatricesList()$jacobian)]
  jacobian_rows <- ctModelUnlist(jacobian, names(jacobian))
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
  ctm$pars <- ctModelStatesAndPARS(
    ctm$pars, statenames = ctm$latentNames, tdprednames = ctm$TDpredNames
  )
  T0VARredundancies(ctm)
}

# A reference inside a matrix the engine evaluates once, from the parameters
# alone -- T0MEANS and T0VAR are the two -- rather than per row or per state.
#
# `PARS[k,1]` there is not dynamics, it is aliasing. Declaring a parameter in
# PARS and naming it in several cells is how ctsem shares one parameter
# between them, and `.ctModelIntOverPop()` rewrites every such use into a
# reference. The value is known as soon as the parameter vector is; there is
# nothing to evaluate per row. This used to be refused wholesale, by a check
# that fired on any `[` in the cell and reported every case as "does not
# support a state-dependent expression in ..." -- including, as here, cells
# with no state in them at all. What it actually stopped was a
# behaviour-genetics model whose T0MEANS is four cells sharing two PARS
# parameters.
#
# So resolve it the way the T0MEANS *state* references above are resolved:
# compose the transform text at model-build time, leaving an ordinary
# transformed free parameter that needs no engine support. Refuse only what
# genuinely cannot be resolved, and say which of the four reasons it is.
#
# ONE distinct reference per cell, and that is a hard constraint rather than
# a convenience. `parameter_transforms.jl` states that each regular transform
# reads exactly one entry of the parameter vector, `sp.parnumber[tf_idx]`,
# and `adjoint_parameters.jl` verifies it when the adjoint workspace is
# built. A cell over two raw parameters has no single `parnumber` to give, so
# it is refused rather than mis-indexed -- which would be a wrong gradient
# rather than an error.
#
# The referencing cell's own `transform` is dropped, not composed on top of
# the referenced value. Measured against stan, which is the reference: for
# `PARS = 'a'` with `T0MEANS = c('a','a')`, whose cells carry ctsem's default
# `10 * param`, stan reports T0MEANS[1,1] equal to PARS[1,1] and not ten
# times it. The comment above `use_mapped` says the same thing from the other
# side -- a parameter reused in another matrix must not reapply that cell's
# local transform.
#
# Circularity is caught here by `visited`, and again, earlier and with a
# better message, by `ctModelCycleCheck.R` when the model is built. Both are
# worth having: this one also covers a chain that only exists after
# `.ctModelIntOverPop()` has rewritten the cells.
.CT_JULIA_PARSREF <- "\\bPARS\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]"

.ctJuliaResolveStaticRefs <- function(p, rows) {

  refs_of <- function(x) unique(gsub("[[:space:]]", "",
    regmatches(x, gregexpr(.CT_JULIA_PARSREF, x, perl = TRUE))[[1L]]))
  index_of <- function(ref) as.integer(regmatches(ref, gregexpr("\\d+", ref))[[1L]])
  cellname <- function(i) paste0(p$matrix[i], "[", p$row[i], ",", p$col[i], "]")

  # `resolve` and `compose` are mutually recursive -- a PARS cell may itself
  # be an expression over another PARS cell -- which works because both are
  # looked up when they are called, not when they are defined.
  resolve <- function(key, visited) {
    if (key %in% visited) return(list(kind = "no", why = paste0(
      "It is a circular reference: ", paste(c(visited, key), collapse = " -> "),
      ". Break the chain by making one of those cells a fixed value")))
    kc <- index_of(key)
    row_k <- which(p$matrix == "PARS" & p$row == kc[1L] & p$col == kc[2L])
    if (!length(row_k)) return(list(kind = "no", why = paste0(
      "It references ", key, ", which the model does not define")))
    row_k <- row_k[1L]
    text_k <- p$param[row_k]
    if (is.na(text_k)) {
      value <- suppressWarnings(as.numeric(p$value[row_k]))
      if (is.na(value)) return(list(kind = "no", why = paste0(
        key, " is neither a free parameter nor a fixed value")))
      return(list(kind = "fixed", value = value))
    }
    if (grepl("[", text_k, fixed = TRUE)) return(compose(text_k, c(visited, key), key))
    parnumber <- p$parnumber[row_k]
    if (is.na(parnumber)) return(list(kind = "no", why = paste0(
      key, " has no parameter number, so there is nothing to reference")))
    inner <- if (is.na(p$transform[row_k])) "param" else as.character(p$transform[row_k])
    list(kind = "free", parnumber = parnumber,
      text = gsub("\\bparam\\b", paste0("param[", parnumber, "]"), inner, perl = TRUE))
  }

  compose <- function(text, visited, what) {
    if (grepl("tdpreds[", text, fixed = TRUE)) return(list(kind = "no", why = paste0(
      what, " depends on the time-dependent predictor data, which varies by row, ",
      "and this matrix is evaluated once from the parameters")))
    if (grepl("\\bstate\\s*\\[", text, perl = TRUE)) return(list(kind = "no", why = paste0(
      what, " depends on a latent state")))
    refs <- refs_of(text)
    if (!length(refs)) return(list(kind = "no", why = paste0(
      what, " contains an index the julia backend does not recognise: ", text)))
    if (length(refs) > 1L) return(list(kind = "no", why = paste0(
      what, " combines ", length(refs), " PARS cells (", paste(refs, collapse = ", "),
      "). The engine's transform layer reads exactly one raw parameter per cell, ",
      "so this one cannot be resolved to a single parameter")))
    inner <- resolve(refs[1L], visited)
    if (identical(inner$kind, "no")) return(inner)
    outer <- gsub(.CT_JULIA_PARSREF, "param", text, perl = TRUE)
    if (identical(inner$kind, "fixed")) {
      folded <- try(eval(parse(text = gsub("\\bparam\\b",
        format(inner$value, digits = 17, scientific = FALSE), outer, perl = TRUE))),
        silent = TRUE)
      if (inherits(folded, "try-error") || !is.numeric(folded) || length(folded) != 1L ||
          !is.finite(folded)) {
        return(list(kind = "no", why = paste0(what,
          " resolves to a fixed value the expression could not be evaluated at: ", text)))
      }
      return(list(kind = "fixed", value = as.numeric(folded)))
    }
    list(kind = "free", parnumber = inner$parnumber,
      text = gsub("\\bparam\\b", paste0("(", inner$text, ")"), outer, perl = TRUE))
  }

  for (i in rows) {
    r <- compose(as.character(p$param[i]), character(0), cellname(i))
    if (identical(r$kind, "no")) {
      stop("Julia backend cannot resolve ", cellname(i), " = ", p$param[i], ". ",
        r$why, ".", call. = FALSE)
    }
    p$param[i] <- NA_character_
    if (identical(r$kind, "fixed")) {
      p$value[i] <- r$value
      p$transform[i] <- NA_character_
      p$parnumber[i] <- NA_integer_
    } else {
      # Already in final `param[k]` form. This row is not in `free` -- that
      # was computed with `!grepl("[", param)` and this cell's param had one
      # -- so the substitution loop below leaves the text alone, which is
      # what keeps `param[3]` from becoming `param[3][3]`.
      p$value[i] <- NA_real_
      p$parnumber[i] <- r$parnumber
      p$transform[i] <- r$text
    }
  }
  p
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
    tform("param", transform, multiplier, meanscale, offset,
      inneroffset, singletext = TRUE)
  }
  # The model's own transform text, matched by cell rather than by row order.
  # `ctFit` records it before `ctModelTransformsToNum` reduces each transform
  # to four fitted numbers; see the comment there for what that reduction
  # loses and why it matters.
  exact <- rep(NA_character_, nrow(p))
  if (!is.null(ctm$transformtext)) {
    recorded <- ctm$transformtext
    exact <- recorded$text[match(paste(p$matrix, p$row, p$col, sep = "
"),
      paste(recorded$matrix, recorded$row, recorded$col, sep = "
"))]
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

  # The exact transform text wherever it survived, for the free parameters only.
  #
  # `ctModelTransformsToNum` turns each transform string into four numbers by a
  # grid search and keeps the string beside them as `transformtext`;
  # re-rendering from the numbers loses any constant the search could not see,
  # which is every variance diagonal's `1e-10` floor -- see the comment there
  # for what a floorless variance does to the filter.
  #
  # Free parameters only, because a cell that is not one has no `param[k]` to
  # write and its text still says `param`. Under `intoverpop` the CINT cells
  # become references to augmented states and are exactly that: substituting
  # there left `10 * param` in the table with no index on it, which the engine
  # compiled into a closure over the whole parameter vector.
  substitute_exact <- free & !is.na(exact)
  p$transform[substitute_exact] <- exact[substitute_exact]
  if (!is.null(ctm$modelmats$matsetup)) {
    matrix_codes <- .ctMatricesList()$all
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
    #
    # Which row is canonical comes from `matsetup`, but the transform *text*
    # comes from the model, and re-rendering it from `matvalues` instead is
    # lossy. `ctModelTransformsToNum` recovers `matvalues` from the model's
    # transform string by fitting the four numbers to it by least squares and
    # then rounding to six decimal places, so any constant below 5e-7 is gone:
    # every variance diagonal carries a `1e-10` floor -- `1e-10 + 5 *
    # log1p_exp(2 * param)` -- and re-rendering gives `0 + 5 * (...)` back.
    # DRIFT's `1e-06` floor survives the rounding, which is why this showed up
    # only in the variances.
    #
    # A variance with no floor reaches *exactly* zero, because `log1p_exp(x)`
    # is exactly zero once `1 + exp(x)` rounds to one, at about `x = -37` --
    # which one subject's TI-predictor offset is quite enough to reach. What
    # follows is division by that zero and derivatives of `sqrt` at it, and
    # since `isfinite` on a dual number tests only its value, the NaN travels
    # in the partials with nothing to stop it. It surfaced as trial points
    # rejected for a non-finite gradient -- 64 of 229 in one Laplace fit --
    # and a fit that stopped 0.6 log units short.
    setup_parameter <- as.integer(setup$param)
    setup_when <- as.integer(setup$when)
    # when == 0 only. A `when == 100` row is a PARS cell whose text references
    # a state, and its `param` column holds that *state* index rather than a
    # parameter number, so including it let a carrier cell claim a real
    # parameter and replace its transform. This mirrors the same correction in
    # stan's parvectform; see the comment there for why 100 was ever accepted.
    candidate <- which(setup_parameter > 0L & setup_when %in% 0L)
    canonical <- rep(NA_character_, max(c(0L, setup_parameter), na.rm = TRUE))
    own_row <- match(setup_key, parameter_key)
    for (row in candidate) {
      parameter <- setup_parameter[row]
      if (!is.na(canonical[parameter])) next
      own <- own_row[row]
      text <- if (!is.na(own) && !is.na(p$transform[own]))
        as.character(p$transform[own]) else NA_character_
      canonical[parameter] <- if (!is.na(text)) text else
        render_transform(setup$transform[row], values$multiplier[row],
          values$meanscale[row], values$offset[row], values$inneroffset[row])
    }
    use_canonical <- free & !is.na(p$parnumber) & p$parnumber <= length(canonical) &
      !is.na(canonical[p$parnumber])
    p$transform[use_canonical] <- canonical[p$parnumber[use_canonical]]
  }
  # A T0MEANS cell that references exactly one other latent state -- the
  # standard idiom for a stable latent intercept, ctsem's state-space
  # replacement for MANIFESTTRAITVAR -- is not state-dependent in the sense
  # the predict/update/td groups below handle. The state it reads is fully
  # determined at t0 by other free parameters and fixed values, with no
  # runtime dependency on the Kalman filter's evolving state: `ctStanModel`'s
  # cycle check (see `ctModelCycleCheck.R`) has already refused any model
  # where that is not true, so every reference here resolves, in a finite
  # number of steps, to a free parameter or a fixed constant.
  #
  # Rather than adding a fourth state-dependent transform group to the engine
  # -- mirrored through the adjoint tape, the sampler and the summary
  # routines, for predict/update/td -- that chain is resolved here, at
  # model-build time, by composing transform text: T0MEANS[1,1] =
  # tform_outer(state[3]) becomes T0MEANS[1,1] = tform_outer(tform_inner(param[pn])),
  # the same free-parameter transform machinery every other T0MEANS cell
  # already uses. This is numerically identical to Stan's own resolution,
  # which evaluates state 3 first and then applies tform_outer to that number
  # (see the two-pass t0 block in `ctModelWriter.R`) -- composing the
  # transforms symbolically is the same computation done in one pass instead
  # of two.
  #
  # Individual variation on the *referenced* state (an indvarying T0MEANS
  # cell feeding a state-reference chain) is not propagated to the composed
  # cell here -- out of scope for the idiom this resolves, which by
  # construction carries no individual variation of its own.
  .ctJuliaStateIndex <- function(x) {
    m <- regmatches(x, regexpr("\\bstate\\[\\d+\\]", x, perl = TRUE))
    if (!length(m) || !nzchar(m)) return(NA_integer_)
    as.integer(gsub("\\D", "", m))
  }
  t0_resolve_cache <- new.env(parent = emptyenv())
  resolve_t0_state <- function(state_idx, visited) {
    key <- as.character(state_idx)
    cached <- t0_resolve_cache[[key]]
    if (!is.null(cached)) return(cached)
    if (state_idx %in% visited) {
      stop("Circular t0 state reference at state ", state_idx,
        " -- this should already have been refused when the model was built.",
        call. = FALSE)
    }
    row_k <- which(p$matrix == "T0MEANS" & p$row == state_idx & p$col == 1L)
    if (!length(row_k)) {
      stop("T0MEANS[", state_idx, ",1] not found while resolving a state reference.",
        call. = FALSE)
    }
    row_k <- row_k[1L]
    text_k <- p$param[row_k]
    result <- if (is.na(text_k)) {
      # A fixed constant: nothing to compose, this state's t0 value is a number.
      list(kind = "fixed", value = as.numeric(p$value[row_k]))
    } else if (simpleStateCheck(text_k)) {
      # This state's own t0 value is itself a reference to another state --
      # resolve that one first, then apply this cell's own transform to it.
      inner <- resolve_t0_state(.ctJuliaStateIndex(text_k), c(visited, state_idx))
      outer_text <- gsub("\\bstate\\[\\d+\\]", "param", text_k, perl = TRUE)
      if (identical(inner$kind, "fixed")) {
        list(kind = "fixed", value = eval(parse(text = gsub("\\bparam\\b",
          format(inner$value, digits = 17, scientific = FALSE), outer_text, perl = TRUE))))
      } else {
        inner_text <- if (is.na(inner$text)) "param" else inner$text
        list(kind = "free", parnumber = inner$parnumber,
          text = gsub("\\bparam\\b", paste0("(", inner_text, ")"), outer_text, perl = TRUE))
      }
    } else if (grepl("[", text_k, fixed = TRUE)) {
      stop("T0MEANS[", state_idx, ",1] is a state-dependent expression the ",
        "julia backend cannot resolve at model-build time: ", text_k, ".", call. = FALSE)
    } else {
      # An ordinary free parameter.
      pn <- p$parnumber[row_k]
      if (is.na(pn)) {
        stop("T0MEANS[", state_idx, ",1] must be a free parameter or a fixed ",
          "value to be referenced by another state's t0 cell.", call. = FALSE)
      }
      list(kind = "free", parnumber = pn,
        text = if (is.na(p$transform[row_k])) NA_character_ else as.character(p$transform[row_k]))
    }
    t0_resolve_cache[[key]] <- result
    result
  }
  t0state_rows <- which(p$matrix == "T0MEANS" & !is.na(p$param) &
    vapply(p$param, function(x) !is.na(x) && simpleStateCheck(x), logical(1)))
  if (length(t0state_rows)) {
    # Resolve every reference first, entirely from the table's original
    # (pre-resolution) contents -- then apply the results in a second pass.
    # Mutating `p` between resolutions would let one cell's resolved,
    # already-rewritten row be misread as a fixed value by a later lookup of
    # the same state.
    #
    # Resolved on the state *this row defines* (`p$row[i]`), not the state its
    # text references: `resolve_t0_state` reads that row's own text and
    # composes its own outer transform on top of what the reference resolves
    # to. Resolving the referenced state directly and using it unchanged would
    # silently drop this row's own transform -- exactly right for a bare
    # reference, silently wrong for anything else, e.g. `0.5 * eta2`.
    resolved <- lapply(t0state_rows, function(i)
      resolve_t0_state(p$row[i], integer(0)))
    for (j in seq_along(t0state_rows)) {
      i <- t0state_rows[j]
      r <- resolved[[j]]
      p$param[i] <- NA_character_
      if (identical(r$kind, "fixed")) {
        p$value[i] <- r$value
        p$transform[i] <- NA_character_
        p$parnumber[i] <- NA_integer_
      } else {
        inner_text <- if (is.na(r$text)) "param" else r$text
        p$value[i] <- NA_real_
        p$parnumber[i] <- r$parnumber
        # Substituted to its final `param[pn]` form directly, rather than left
        # for the `free`-only substitution loop below: `free` was computed
        # before this row was resolved, so it does not include this row.
        p$transform[i] <- gsub("\\bparam\\b", paste0("param[", r$parnumber, "]"),
          inner_text, perl = TRUE)
      }
    }
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
    p <- .ctJuliaResolveStaticRefs(p, which(unsupported_dynamic))
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
    if (length(effect_columns)) entry[effect_columns] <- "FALSE"
    table <- rbind(table, entry)
  }
  table
}

# Which parameter cells are written by a state-dependent expression. Shared by
# both random-effect routes, which need the same answer for the same reason:
# a cell whose value depends on the current state cannot be treated as a
# constant of the filter.
.ctJuliaRewrittenCells <- function(table) {
  predict_transform <- replace(table$predicttransform, is.na(table$predicttransform), "")
  update_transform <- replace(table$updatetransform, is.na(table$updatetransform), "")
  td_transform <- replace(table$tdtransform, is.na(table$tdtransform), "")
  table[grepl("state\\[", predict_transform) |
    grepl("state\\[", update_transform) |
    grepl("state\\[", td_transform),
    c("matrix", "row", "col"), drop = FALSE]
}

# Which raw parameters the engine flagged as saturated, by name.
#
# The engine reports indices, and never an empty vector -- a zero-length
# vector deadlocks the JuliaConnectoR bridge, so "nothing saturated" arrives
# as the single index 0. "Raw parameter 10" is not something a user can act
# on; "rawcor_mm2__mm1" is, and it is usually the whole diagnosis.
.ctJuliaSaturatedNames <- function(result, model_spec, npar) {
  index <- suppressWarnings(as.integer(result$saturated_parameters))
  index <- index[!is.na(index) & index >= 1L & index <= npar]
  if (!length(index)) return(character())
  names <- try(.ctBackendRawParameterNames(list(model_spec = model_spec), npar),
    silent = TRUE)
  if (inherits(names, "try-error") || length(names) != npar) {
    return(paste0("raw[", index, "]"))
  }
  as.character(names[index])
}

# Which random-effect correlations, if any, ended the fit on their cap.
#
# Reported as a data frame naming the parameter pair and the level, because
# "a correlation hit the boundary" is only actionable if you know which one.
.ctJuliaBoundaryReport <- function(fit, model_spec, result) {
  empty <- data.frame(level = character(), param = character(),
    correlation = numeric(), stringsAsFactors = FALSE)
  laplace <- model_spec$laplace
  if (is.null(laplace)) return(empty)
  module <- .ctJuliaModule(model_spec$project)
  hit <- try(.ctBackendJuliaValue(module$ctsem_laplace_boundary(
    .ctJuliaObjective(fit), .ctJuliaNumericVector(fit$estimate$raw))), silent = TRUE)
  if (inherits(hit, "try-error") || !length(hit$level)) return(empty)
  rows <- lapply(seq_along(hit$level), function(i) {
    level <- laplace$levels[[hit$level[i]]]
    names <- level$param
    lower <- which(lower.tri(diag(max(1L, length(names)))), arr.ind = TRUE)
    pair <- if (hit$position[i] <= nrow(lower))
      paste0(names[lower[hit$position[i], 1L]], "__", names[lower[hit$position[i], 2L]])
      else paste0("correlation", hit$position[i])
    # The capped value, which is what the model used -- transforming the raw
    # coordinate as it stands would report the correlation the optimizer was
    # heading for rather than the one it was held at.
    capped <- max(-.ctJuliaCorrelationCap(), min(.ctJuliaCorrelationCap(), hit$value[i]))
    data.frame(level = level$name, param = pair,
      correlation = 2 / (1 + exp(-capped)) - 1, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# The engine's correlation cap, read rather than duplicated.
.ctJuliaCorrelationCap <- function() {
  value <- try(JuliaConnectoR::juliaEval("ContinuousTimeSEM._LAPLACE_COR_CAP[]"),
    silent = TRUE)
  if (inherits(value, "try-error")) 5.2933 else as.numeric(value)[1L]
}

# Every index the specification references must lie inside the raw parameter
# vector the fit will actually allocate.
#
# This invariant is here because breaking it has been the single recurring bug
# of this feature, three times over: a level's population scales, and then the
# TI-predictor coefficients, each indexed past the end of a vector sized from
# only part of the layout. Neither failed usefully. The first left a level
# pinned at its starting value while the summary printed the transform of that
# value as though it were an estimate -- the same number for every dataset. The
# second threw a bounds error that the optimizer's invalid-point guard
# swallowed, which then surfaced as a line search crashing on an unrelated
# assertion. Both were found by recovery studies, days later.
#
# Checking it costs nothing and turns that whole class into an immediate,
# specific error.
.ctJuliaCheckLayout <- function(table, laplace, ti_effects, npar, ti_missing = NULL) {
  used <- list(
    `model parameters` = as.integer(table$parnumber),
    `TI-predictor coefficients` = as.integer(ti_effects$coefficient))
  if (!is.null(ti_missing) && nrow(ti_missing)) {
    used[["sampled TI-predictor values"]] <- as.integer(ti_missing$parameter)
  }
  if (!is.null(laplace)) {
    for (level in laplace$levels) {
      used[[paste0("'", level$name, "' population scales")]] <- as.integer(level$sd_index)
      used[[paste0("'", level$name, "' correlations")]] <- as.integer(level$cor_index)
      used[[paste0("'", level$name, "' varying parameters")]] <- as.integer(level$re_index)
    }
  }
  for (what in names(used)) {
    index <- used[[what]]
    index <- index[!is.na(index)]
    if (length(index) && max(index) > npar) {
      stop("Internal layout error: the ", what, " reach raw parameter ",
        max(index), " but the raw vector holds ", npar,
        ". This is a bug in ctsem rather than in the model; please report it.",
        call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Which parameters vary at a given level, and the column that says so.
#
# The subject level keeps `indvarying`, unchanged and backward compatible. Each
# grouping level above it gets `indvarying_<idname>` -- explicit rather than
# positional, so a model carrying three levels reads as three named columns
# rather than as a matrix nobody can check by eye.
.ctJuliaLevelColumn <- function(model, level) {
  if (level == 1L) return("indvarying")
  paste0("indvarying_", model$groupIDnames[level - 1L])
}

# The matching sdscale column. One value per id element, so a level scales its
# own population sd rather than borrowing the subject level's.
.ctJuliaLevelScaleColumn <- function(model, level) {
  if (level == 1L) return("sdscale")
  paste0("sdscale_", model$groupIDnames[level - 1L])
}

# Each subject's group at each level, and a check that the nesting is strict.
#
# Strict nesting is what makes the hierarchy a tree, and a tree is what lets the
# integral factorise into one unit per outermost group. A subject appearing in
# two studies is not a hierarchy this method can integrate, so it is refused
# here rather than producing a quietly wrong answer later.
.ctJuliaHierarchy <- function(dat, model) {
  subjects <- unique(dat[[model$subjectIDname]])
  levels <- list(list(name = model$subjectIDname,
    group = seq_along(subjects), ngroups = length(subjects),
    labels = subjects))
  for (nm in model$groupIDnames) {
    if (!nm %in% colnames(dat)) {
      stop("Grouping id column '", nm, "' not found in the data.", call. = FALSE)
    }
    rows <- match(subjects, dat[[model$subjectIDname]])
    labels <- dat[[nm]][rows]
    # One group per subject, and the same one on every row of that subject.
    perrow <- split(dat[[nm]], dat[[model$subjectIDname]])
    ambiguous <- names(perrow)[vapply(perrow, function(x) length(unique(x)) > 1L, logical(1))]
    if (length(ambiguous)) {
      stop("Subject(s) ", paste(utils::head(ambiguous, 5), collapse = ", "),
        " have more than one value of '", nm, "'. The hierarchy must be strictly ",
        "nested: each subject belongs to exactly one group at every level.",
        call. = FALSE)
    }
    unique_labels <- unique(labels)
    levels[[length(levels) + 1L]] <- list(name = nm,
      group = match(labels, unique_labels), ngroups = length(unique_labels),
      labels = unique_labels)
  }
  .ctJuliaCheckNesting(levels)
  levels
}

# Every level must nest inside the next one out.
#
# The per-level check above catches a subject that spans two groups at one
# level. It does not catch two levels that cross *each other* -- pupils in
# classes and in neighbourhoods, where a class draws from several
# neighbourhoods -- and that is the case the engine cannot represent at all.
#
# Nesting is what makes a unit's curvature a tree: a block couples to another
# only when a member depends on both, and with nesting that means one contains
# the other, so eliminating innermost-first produces no fill-in and the cost is
# linear rather than cubic in a group's members. Crossing puts a cycle in that
# graph, and the block factorisation and the nested quadrature both stop being
# the right decomposition. Refusing is honest; silently treating one grouping
# as nested inside the other would not be.
.ctJuliaCheckNesting <- function(levels) {
  if (length(levels) < 2L) return(invisible(NULL))
  for (inner in seq_len(length(levels) - 1L)) {
    for (outer in seq(inner + 1L, length(levels))) {
      split_outer <- split(levels[[outer]]$group, levels[[inner]]$group)
      crossed <- vapply(split_outer, function(x) length(unique(x)) > 1L, logical(1))
      if (any(crossed)) {
        offending <- levels[[inner]]$labels[as.integer(names(which(crossed)))]
        stop("Grouping '", levels[[inner]]$name, "' is not nested within '",
          levels[[outer]]$name, "': ", sum(crossed), " group(s) of '",
          levels[[inner]]$name, "' span more than one '", levels[[outer]]$name,
          "' (for example ", paste(utils::head(offending, 3), collapse = ", "),
          "). The julia backend represents strictly nested hierarchies only -- ",
          "that is what makes each unit's curvature a tree and its cost linear ",
          "rather than cubic in group size. Crossed designs need a different ",
          "factorisation and are not supported.", call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

# The Laplace route's counterpart to `.ctJuliaAugmentRandomEffects`.
#
# Where the augmented route grows the latent state, so that the ordinary filter
# integrates the random effects out alongside the dynamic states, this route
# leaves the model's state space alone and describes the random effects to the
# engine as *indices*: which raw parameters vary between subjects, and which
# raw parameters say how much. The engine then integrates them out per subject.
#
# The raw vector's layout is Stan's, deliberately: population means, then
# population scales, then lower-triangular correlation coordinates, then TI
# predictor effects. `.ctBackendPriorSpec` maps ctsem's priors onto exactly
# that layout, and `.ctJuliaAugmentRandomEffects` already reproduces it, so
# reusing it means priors, uncertainty and the summary need no Laplace-specific
# case -- and the two routes' raw vectors are directly comparable, which is
# what makes an augmented-versus-Laplace check meaningful at all.
.ctJuliaLaplaceSpec <- function(model, table, prepared_data = NULL, dat) {
  base_npar <- suppressWarnings(max(c(0L, as.integer(table$parnumber)), na.rm = TRUE))
  varying <- integer()
  sdscale <- numeric()
  if (length(prepared_data$indvaryingindex)) {
    # Stan's own `indvaryingindex`, built by .ctPrepareData from the same matsetup.
    # Preferred when present so both backends agree on the set and its order.
    varying <- as.integer(prepared_data$indvaryingindex)
    sdscale <- as.numeric(prepared_data$sdscale)
  } else if (!is.null(model$modelmats$matsetup)) {
    setup <- as.data.frame(model$modelmats$matsetup)
    values <- as.data.frame(model$modelmats$matvalues)
    rows <- which(as.integer(setup$indvarying) > 0L & as.integer(setup$param) > 0L)
    rows <- rows[!duplicated(as.integer(setup$param)[rows])]
    varying <- as.integer(setup$param)[rows]
    sdscale <- as.numeric(values$sdscale)[rows]
  } else {
    # Neither, which is a model straight from `ctModel()`. `ctFit()` always
    # supplies both -- it builds `modelmats` and hands `.ctPrepareData()`'s
    # output over as `prepared_data` -- so this branch is for the callers that
    # prepare a specification without fitting one, `ctIdentify()` among them.
    # Without it `nrandom` came out zero on a model that plainly has random
    # effects and `intoverpop='laplace'` was refused as "no parameters are
    # marked indvarying", which is the opposite of what the model said.
    #
    # `table` is the canonical parameter table, so its `indvarying` is the
    # flag after `.ctJuliaCanonicalModel()` has had its say, and its
    # `parnumber` is the raw index. `sdscale` is not carried on the table, so
    # it is read from the model's own `pars` by parameter name.
    hit <- which(table$indvarying %in% TRUE & !is.na(table$parnumber) &
      table$parnumber > 0L)
    varying <- sort(unique(as.integer(table$parnumber[hit])))
    sdscale <- rep(1, length(varying))
    if (!is.null(model$pars$sdscale)) {
      names_at_index <- .ctJuliaLaplaceNames(table, varying)
      matched <- match(names_at_index, as.character(model$pars$param))
      supplied <- suppressWarnings(as.numeric(model$pars$sdscale[matched]))
      sdscale[is.finite(supplied)] <- supplied[is.finite(supplied)]
    }
  }
  if (length(sdscale) != length(varying)) sdscale <- rep(1, length(varying))
  usable <- !is.na(varying) & varying > 0L & varying <= base_npar
  varying <- varying[usable]
  sdscale <- sdscale[usable]
  sdscale[!is.finite(sdscale)] <- 1

  hierarchy <- .ctJuliaHierarchy(dat, model)

  # Level 1 is the subject level and uses the `indvarying` set found above.
  # Outer levels read their own column.
  setup <- if (!is.null(model$modelmats$matsetup)) as.data.frame(model$modelmats$matsetup) else NULL
  cursor <- base_npar
  levels <- list()
  for (l in seq_along(hierarchy)) {
    if (l == 1L) {
      lv_varying <- varying; lv_scale <- sdscale
    } else {
      column <- .ctJuliaLevelColumn(model, l)
      lv_varying <- integer(); lv_scale <- numeric()
      if (column %in% names(model$pars)) {
        names_at_level <- unique(model$pars$param[model$pars[[column]] %in% TRUE &
          !is.na(model$pars$param) & is.na(model$pars$value)])
        hit <- which(!is.na(table$parnumber) & table$param %in% names_at_level)
        lv_varying <- sort(unique(as.integer(table$parnumber[hit])))
        lv_varying <- lv_varying[lv_varying <= base_npar]
        # This level's own sdscale, not the subject level's: `sdscale` is one
        # value per id element and each level uses its own.
        scalecolumn <- .ctJuliaLevelScaleColumn(model, l)
        lv_scale <- rep(1, length(lv_varying))
        if (scalecolumn %in% names(model$pars)) {
          for (position in seq_along(lv_varying)) {
            row <- which(!is.na(table$parnumber) & table$parnumber == lv_varying[position])
            name <- if (length(row)) as.character(table$param[row[1L]]) else NA_character_
            match_row <- which(model$pars$param %in% name)
            if (length(match_row)) {
              value <- as.numeric(model$pars[[scalecolumn]][match_row[1L]])
              if (is.finite(value)) lv_scale[position] <- value
            }
          }
        }
      }
    }
    k <- length(lv_varying)
    noff <- as.integer(k * (k - 1L) / 2L)
    levels[[l]] <- list(
      name = hierarchy[[l]]$name,
      re_index = as.integer(lv_varying),
      sd_index = if (k) as.integer(cursor + seq_len(k)) else integer(),
      cor_index = if (noff) as.integer(cursor + k + seq_len(noff)) else integer(),
      sd_scale = as.numeric(lv_scale),
      param = .ctJuliaLaplaceNames(table, lv_varying),
      nrandom = as.integer(k),
      group = as.integer(hierarchy[[l]]$group),
      ngroups = as.integer(hierarchy[[l]]$ngroups),
      labels = hierarchy[[l]]$labels)
    cursor <- cursor + k + noff
  }

  total <- sum(vapply(levels, function(x) x$nrandom, integer(1)))
  list(
    levels = levels,
    # Level-one fields kept at the top for every existing single-level caller.
    re_index = levels[[1]]$re_index,
    sd_index = levels[[1]]$sd_index,
    cor_index = levels[[1]]$cor_index,
    sd_scale = levels[[1]]$sd_scale,
    param = levels[[1]]$param,
    nrandom = as.integer(total),
    nlevels = length(levels),
    base_npar = as.integer(base_npar),
    npar = as.integer(cursor)
  )
}

# The model's own name for each varying parameter, so every later report can
# say which parameter a random effect belongs to without re-deriving it.
.ctJuliaLaplaceNames <- function(table, varying) {
  vapply(as.integer(varying), function(parameter) {
    hit <- which(!is.na(table$parnumber) & table$parnumber == parameter &
      !is.na(table$param))
    if (!length(hit)) NA_character_ else as.character(table$param[hit[1L]])
  }, character(1L))
}

# Stan's `derrind` (ctData.R:310-315) reduction, ported: among the states in
# `1:span` (the augmented layout when there is one, else just the original
# states), which ones need their own row/column in the diffusion covariance
# and Lyapunov solve. Two steps, matching Stan exactly:
#
#  1. A state has "structurally zero" diffusion only when every DIFFUSION cell
#     in its row and column is a *fixed* value of exactly zero -- `table$value`
#     is NA for every free parameter and for every state-dependent transform
#     (see `.ctJuliaParameterTable`/`.ctJuliaCanonicalModel`, which leave
#     `value` NA whenever the cell is anything but a literal constant), so
#     testing `!is.na(value) && value == 0` never misclassifies a parameter
#     that merely *evaluates* to zero at the current draw. This mirrors
#     `listOfMatrices(ctm$pars)$DIFFUSION` where a free cell holds a parameter
#     name (never equal to the string "0") and only a literal fixed zero
#     prints as "0".
#  2. Zero own diffusion is not enough to drop a state: if its drift/Jacobian
#     (`JAx`, or `DRIFT` when no JAx rows exist -- both padded to `span` by
#     the caller) structurally couples it to a diffusing state, covariance
#     still reaches it through the transition, exactly as a matrix
#     exponential of a coupled-but-noise-free block is not itself block
#     diagonal. This is computed as the graph-connected-component closure
#     Stan gets via `JAxsubsets` (`expmGetSubsets`, ctData.R:1-23): coupling
#     is treated as undirected (`JAx[i,j]` nonzero links i and j both ways),
#     because that is what `expmGetSubsets` does -- a row-only edge still
#     grows the same connected component there.
#
# Indices beyond `nlatent` (augmented random-effect carriers) are dropped at
# the end, same as Stan's `derrind[derrind <= standata$nlatent]`; carrier
# states are never diffusing themselves, so they only ever enter as
# intermediate coupling nodes on the way to another original state, and this
# still finds them because the closure runs over the full `span`.
.ctJuliaDerrind <- function(table, nlatent, span) {
  if (span <= 0L) return(integer())
  structurally_nonzero <- function(matrix_name) {
    m <- matrix(FALSE, span, span)
    rows <- which(table$matrix == matrix_name & !is.na(table$row) & !is.na(table$col) &
      table$row >= 1L & table$row <= span & table$col >= 1L & table$col <= span)
    for (i in rows) {
      zero <- !is.na(table$value[i]) && isTRUE(as.numeric(table$value[i]) == 0)
      if (!zero) m[table$row[i], table$col[i]] <- TRUE
    }
    m
  }
  diffusion <- structurally_nonzero("DIFFUSION")
  coupling <- if (any(table$matrix == "JAx")) structurally_nonzero("JAx") else structurally_nonzero("DRIFT")
  reach <- apply(diffusion, 1L, any) | apply(diffusion, 2L, any)
  if (!any(reach)) return(integer())
  repeat {
    grown <- reach
    for (state in which(reach)) grown <- grown | coupling[state, ] | coupling[, state]
    if (!any(grown & !reach)) break
    reach <- grown
  }
  sort(which(reach)[which(reach) <= nlatent])
}

.ctJuliaAugmentRandomEffects <- function(model) {
  original_nlatent <- model$n.latent
  prepared_augmentation <- !is.null(model$intoverpopindvaryingindex)
  has_random_effects <- prepared_augmentation || any(model$pars$indvarying %in% TRUE &
    is.na(suppressWarnings(as.numeric(model$pars$value))), na.rm = TRUE)
  if (!has_random_effects) {
    no_re_table <- .ctJuliaParameterTable(model)
    return(list(parameter_table = no_re_table,
      nlatent = original_nlatent, nlatent_augmented = original_nlatent,
      dynamic_state_indices = .ctJuliaDerrind(no_re_table, original_nlatent, original_nlatent),
      random_effects = data.frame(), rewritten_cells = data.frame()))
  }

  expanded <- if (prepared_augmentation) model else {
    prepared <- ctModel0DRIFT(model, model$continuoustime)
    prepared$pars <- ctModelStatesAndPARS(prepared$pars,
      statenames = prepared$latentNames, tdprednames = prepared$TDpredNames)
    .ctModelIntOverPop(prepared)
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
    t0means_code <- .ctMatricesList()$all[["T0MEANS"]]
    t0means_rows <- setup$matrix == t0means_code & setup$col == 1L
    match_position <- match(augmented_indices, setup$row[t0means_rows])
    t0means_setup_rows <- which(t0means_rows)[match_position]
    t0means_state_scale <- values$multiplier[t0means_setup_rows] * values$meanscale[t0means_setup_rows]
    t0means_state_scale[is.na(t0means_state_scale)] <- 1
  } else if (!is.null(expanded$pars$sdscale)) {
    # `modelmats` is built by ctFit on its way to a backend, so the branch above
    # covers every fit. It is *not* built when a specification is prepared
    # directly -- which is what generation does -- and the fallback of 1 then
    # discarded `sdscale` without saying so: a model asking for a population
    # spread a fifth of the default generated data with the default spread.
    #
    # Reading it from `pars` costs a name lookup and makes the two routes agree.
    scale <- vapply(augmented_indices, function(row) {
      entry <- which(table$matrix == "T0MEANS" & table$row == row & table$col == 1L)
      if (!length(entry) || is.na(table$param[entry[1L]])) return(1)
      match_row <- which(!is.na(expanded$pars$param) &
        as.character(expanded$pars$param) == as.character(table$param[entry[1L]]))
      if (!length(match_row)) return(1)
      value <- suppressWarnings(as.numeric(expanded$pars$sdscale[match_row[1L]]))
      if (!is.finite(value)) 1 else value
    }, numeric(1L))
    random_sd_scale <- scale
  }
  if (length(random_sd_scale) != length(augmented_indices)) {
    stop("Prepared random-effect covariance metadata does not match the augmented state layout.", call. = FALSE)
  }
  # A cell repurposed as a population parameter must not inherit the
  # TI-predictor effect flags of whatever the user originally wrote there.
  # `T0VARredundancies()` (R/ctFit.R) fixes the free T0VAR cells that indvarying
  # T0MEANS makes redundant and clears their `<TIpred>_effect` columns along with
  # `param`, `transform` and `indvarying`; this function then hands those same
  # cells to the population SD and correlation parameters below, and clears the
  # flags again rather than trust that every cell arrived by that route. A stale
  # flag made `.ctJuliaTIEffects()` mint a TI-predictor coefficient for every
  # population parameter -- something the generated Stan model has no
  # counterpart for, since its population block is `rawpopsdbase`/`sqrtpcov`
  # and takes no TI effects. The consequences were one per backend: the engine silently
  # estimated those extra coefficients, and `.ctBackendPriorSpec()` refused the
  # fit outright because its Stan-derived layout was that many parameters short.
  # The cells `.ctJuliaPadMatrix()` created are already FALSE; these are the
  # ones that existed before.
  effect_columns <- grep("_effect$", names(table), value = TRUE)
  covariance_rows <- list()
  # The name of the parameter each varying state carries, read from its T0MEANS
  # cell before that cell is rewritten below. Everything downstream that has to
  # say which parameter a random effect belongs to -- the summary's popsd and
  # rawpopcorr rows, the carrier state names in ctKalman output, ctSubjectPars
  # -- gets it from here rather than re-deriving the augmentation.
  varying_names <- vapply(augmented_indices, function(row) {
    entry <- which(table$matrix == "T0MEANS" & table$row == row & table$col == 1L)
    if (!length(entry) || is.na(table$param[entry[1L]])) NA_character_ else
      as.character(table$param[entry[1L]])
  }, character(1L))
  # Match Stan's unconstrained parameter order exactly: all population scales,
  # then lower-triangular correlation coordinates column by column.
  for (position in seq_along(augmented_indices)) {
    row <- augmented_indices[position]
    col <- row
    index <- which(table$matrix == "T0VAR" & table$row == row & table$col == col)
    length(index) == 1L || stop("Internal Julia augmentation error: missing T0VAR entry.", call. = FALSE)
    # What the model says about this population sd, if anything. RAWPOPVAR is the
    # specification surface (see R/ctModelRawPopVar.R); a number there fixes the
    # cell and a label names the parameter, in place of the positional
    # `julia_popcov_i_j` this used to invent.
    spec <- .ctModelRawPopVarEntry(model, varying_names[position])
    fixedvalue <- .ctModelRawPopVarValue(spec)
    if (is.finite(fixedvalue)) {
      if (fixedvalue < 0) {
        stop("RAWPOPVAR['", varying_names[position], "', '",
          varying_names[position], "'] is ", fixedvalue,
          ". A population standard deviation cannot be negative.",
          call. = FALSE)
      }
      # Converted from the raw parameter scale a RAWPOPVAR entry is written on to
      # the state scale this cell is in, which is `k_i` and nothing else.
      #
      # The free branch produces `k_i * raw_sd`, so a requested raw spread `v`
      # needs `v * k_i` here. `k_i` is exact arithmetic -- the carrier state's
      # own `multiplier * meanscale` -- so the number a user writes is the raw
      # population sd on this route, and the same number is the raw population
      # sd on the Laplace route, which has no `k_i` and reads it directly.
      #
      # It used to be divided by the derivative of the parameter's own
      # transform as well, making the entry a spread on the parameter's
      # *natural* scale. That was exact only for a linear transform and a
      # first-order match otherwise, and it meant the specification surface and
      # the fit disagreed about what the number was: a fit reports the raw
      # population sd (see `.ctBackendRawPopCov`), so a model with a nonlinear
      # transform could not be generated from and fitted with one number.
      table$param[index] <- NA_character_
      table$parnumber[index] <- NA_integer_
      table$value[index] <- fixedvalue * t0means_state_scale[position]
      table$transform[index] <- NA_character_
    } else {
      next_parameter <- next_parameter + 1L
      table$param[index] <- if (is.na(spec) || !nzchar(spec))
        sprintf("julia_popcov_%d_%d", row, col) else spec
      table$parnumber[index] <- next_parameter
      table$value[index] <- NA_real_
      table$transform[index] <- sprintf("%.17g * (1e-10 + %.17g * log1p_exp(2 * param[%d] - 1))",
        t0means_state_scale[position], random_sd_scale[position], next_parameter)
    }
    # Whether the sd is fixed by RAWPOPVAR or free, this row is a population
    # covariance cell and carries no TI predictor effect of its own.
    if (length(effect_columns)) table[index, effect_columns] <- "FALSE"
    covariance_rows[[length(covariance_rows) + 1L]] <- data.frame(
      row = row, col = col,
      parameter = if (is.finite(fixedvalue)) NA_integer_ else next_parameter,
      type = "sd", param = varying_names[position],
      # The factor folded into the sd transform above, kept so the summary can
      # divide it back out: T0cov is in state units, and the random-effects
      # summary needs the raw-parameter sd to perturb the raw vector by (see
      # .ctBackendRandomEffectSummary).
      scale = t0means_state_scale[position]
    )
  }
  if (length(augmented_indices) > 1L) for (column_position in seq_len(length(augmented_indices) - 1L)) {
    for (row_position in (column_position + 1L):length(augmented_indices)) {
      row <- augmented_indices[row_position]
      col <- augmented_indices[column_position]
      index <- which(table$matrix == "T0VAR" & table$row == row & table$col == col)
      length(index) == 1L || stop("Internal Julia augmentation error: missing T0VAR entry.", call. = FALSE)
      spec <- .ctModelRawPopVarEntry(model, varying_names[row_position],
        varying_names[column_position])
      fixedvalue <- .ctModelRawPopVarValue(spec)
      if (is.finite(fixedvalue)) {
        table$param[index] <- NA_character_
        table$parnumber[index] <- NA_integer_
        table$value[index] <- fixedvalue
        table$transform[index] <- NA_character_
      } else {
        next_parameter <- next_parameter + 1L
        table$param[index] <- if (is.na(spec) || !nzchar(spec))
          sprintf("julia_popcov_%d_%d", row, col) else spec
        table$parnumber[index] <- next_parameter
        table$value[index] <- NA_real_
        # Identity, like every other covariance off-diagonal now: the
        # (-1, 1) map lives in `constraincorsqrt1`, which this cell reaches
        # through T0VAR. Leaving the squash here applied it twice and moved a
        # fitted likelihood by 0.88 on a two-latent model.
        table$transform[index] <- sprintf("param[%d]", next_parameter)
      }
      # As for the sd cells above: a population correlation carries no TI
      # predictor effect of its own, fixed by RAWPOPVAR or not.
      if (length(effect_columns)) table[index, effect_columns] <- "FALSE"
      covariance_rows[[length(covariance_rows) + 1L]] <- data.frame(
        row = row, col = col,
        parameter = if (is.finite(fixedvalue)) NA_integer_ else next_parameter,
        type = "correlation",
        param = paste0(varying_names[row_position], "__", varying_names[column_position]),
        scale = 1
      )
    }
  }
  rewritten <- .ctJuliaRewrittenCells(table)
  # `augmented_indices` (= Stan's `intoverpopindvaryingindex`) is every state
  # with population-varying T0VAR: both the newly-created carrier states for
  # non-T0MEANS random effects (DRIFT/CINT/etc., always appended contiguously
  # after the original states -- see `extralatents` in
  # ctModelWriter.R::.ctModelIntOverPop) AND any *original* state whose
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
  # dimension exceeded 4. `table` has already been padded to
  # `nlatent_augmented` for DRIFT/DIFFUSION/JAx/Jtd above, so
  # `.ctJuliaDerrind` now applies Stan's further reduction (`derrind`'s first
  # two steps: drop original states with structurally zero, drift-uncoupled
  # diffusion) on top of this exclusion rather than instead of it -- see that
  # function's comment for what "structurally zero" must mean and why the
  # drift-coupling closure cannot be skipped.
  list(parameter_table = table, nlatent = original_nlatent,
    nlatent_augmented = nlatent_augmented,
    dynamic_state_indices = .ctJuliaDerrind(table, original_nlatent, nlatent_augmented),
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
    # Constancy within subject is not re-tested here: `.ctJuliaPrepare` calls
    # `.ctJuliaValidateTIConstancy` on the same frame first, and unlike this it
    # runs on the prepared-data path too. What is left is the completeness
    # check above, which that one deliberately does not make -- it ignores NAs
    # so that a partly-recorded predictor is reported as missing rather than as
    # varying.
    values[i, ] <- subject_values[1L, ]
  }
  values
}

# A subject-level TI predictor matrix that leaves missingness as `NA` rather
# than erroring. `.ctJuliaTIData()` above keeps erroring for every other
# caller -- this is a separate function rather than a changed one, so that
# nothing about the existing (non-sampling) callers changes.
#
# Used only by the sampling path's imputation model (see
# `.ctJuliaTIMissingSpec` below), which needs to know *which* cells are
# missing, not merely that some are.
.ctJuliaTIDataAllowMissing <- function(dat, model) {
  subject_ids <- unique(dat[[model$subjectIDname]])
  if (!model$n.TIpred) return(matrix(numeric(), nrow = length(subject_ids), ncol = 0L))
  values <- matrix(NA_real_, nrow = length(subject_ids), ncol = model$n.TIpred,
    dimnames = list(NULL, model$TIpredNames))
  for (i in seq_along(subject_ids)) {
    rows <- dat[[model$subjectIDname]] == subject_ids[i]
    subject_values <- as.matrix(dat[rows, model$TIpredNames, drop = FALSE])
    for (j in seq_len(ncol(subject_values))) {
      observed <- subject_values[!is.na(subject_values[, j]), j]
      if (length(observed)) values[i, j] <- observed[1L]
    }
  }
  values
}

# Placeholder written into a missing TI predictor cell before the base
# predictor matrix is sent to the engine. Never read as data: every cell
# holding it is overwritten by the sampled parameter's current value on
# every evaluation (`_ctsem_tipred_vector`, parameter_transforms.jl). This is
# the same finite, absurd-on-sight value the engine already uses for exactly
# this "must be overwritten before read" role (`UNSET_PARAMETER`,
# inst/julia/.../parameters.jl) -- reused for the same reason there, not
# introduced as a new sentinel, and never tested for by the engine (which
# would be the mistake SPEC-tipred-sampling.md warns against; the engine is
# told the missing positions explicitly instead).
.ctJuliaTIMissingPlaceholder <- 99999

# Per-subject means and SDs of the manifest variables, from the observed
# rows only. A subject with zero (mean) or one (SD) observation of a given
# manifest gets `NA` for that cell, which is what makes it drop out of that
# subject's imputation model rather than contribute a fabricated value
# (Charles's 2026-09-03 decision: manifest means *and* SDs join the
# predictor terms, with the same per-subject dropout rule).
.ctJuliaTISubjectManifestSummary <- function(dat, model) {
  subject_ids <- unique(dat[[model$subjectIDname]])
  ybar <- matrix(NA_real_, length(subject_ids), model$n.manifest,
    dimnames = list(NULL, paste0("ybar_", model$manifestNames)))
  ysd <- matrix(NA_real_, length(subject_ids), model$n.manifest,
    dimnames = list(NULL, paste0("ysd_", model$manifestNames)))
  for (i in seq_along(subject_ids)) {
    rows <- dat[[model$subjectIDname]] == subject_ids[i]
    for (m in seq_len(model$n.manifest)) {
      y <- dat[rows, model$manifestNames[m]]
      y <- y[!is.na(y)]
      if (length(y)) ybar[i, m] <- mean(y)
      if (length(y) > 1L) ysd[i, m] <- stats::sd(y)
    }
  }
  list(ybar = ybar, ysd = ysd)
}

# Roughly five complete cases per coefficient (Charles's 2026-09-03 decision):
# an imputation regression with more terms than this ratio supports is not
# trustworthy, so it falls back to a simpler one instead.
.ctJuliaTIMinCasesPerCoef <- 5

# Fit one TI predictor's imputation regression on complete cases, falling
# back progressively when the data cannot support the fuller model: full
# (the other TI predictors, plus per-subject manifest means and SDs if
# `includeOutcome`) -> predictors only -> marginal (the predictor's own mean
# and SD). Coefficients and residual scale are fit once, here, and then held
# fixed for the whole sampling run -- they are not estimated by the sampler
# (Charles's 2026-09-03 decision): the missing value is still a sampled
# parameter with a proper conditional distribution, so its posterior still
# has spread, which is the part that matters.
#
# Warns, once per predictor, naming the fallback actually taken -- but only
# when a fallback actually happened, i.e. the tier used is not the best tier
# the request and the data shape could in principle support. A single
# predictor with `includeOutcome=FALSE` has no "predictors" tier to try at
# all, so going straight to marginal there is the plan, not a fallback, and
# draws no warning.
.ctJuliaTIImputationFit <- function(predictor_name, target, others, ybar, ysd,
  includeOutcome) {
  min_cases <- .ctJuliaTIMinCasesPerCoef
  fit_tier <- function(design) {
    ncoef <- ncol(design) + 1L
    d <- data.frame(.y = target, design, check.names = FALSE)
    cc <- stats::complete.cases(d)
    if (sum(cc) < min_cases * ncoef) return(NULL)
    fit <- try(stats::lm(.y ~ ., data = d[cc, , drop = FALSE]), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    fit
  }

  outcome_design <- cbind(ybar, ysd)
  if (ncol(outcome_design)) {
    outcome_design <- outcome_design[, apply(outcome_design, 2, function(x) any(!is.na(x))), drop = FALSE]
  }
  full_design <- if (includeOutcome) cbind(others, outcome_design) else others

  best_tier <- if (includeOutcome && ncol(outcome_design)) "full" else
    if (ncol(others)) "predictors" else "marginal"

  marginal_cc <- !is.na(target)
  marginal_fit <- stats::lm(target[marginal_cc] ~ 1)

  fit <- NULL; tier <- NULL
  if (identical(best_tier, "full") && ncol(full_design)) {
    fit <- fit_tier(full_design)
    if (!is.null(fit)) tier <- "full"
  }
  if (is.null(fit) && ncol(others)) {
    fit <- fit_tier(others)
    if (!is.null(fit)) tier <- "predictors"
  }
  if (is.null(fit)) {
    fit <- marginal_fit
    tier <- "marginal"
  }

  if (!identical(tier, best_tier)) {
    warning(sprintf(paste0("TI predictor '%s': not enough complete cases for the %s ",
      "imputation model (roughly %d complete cases needed per coefficient) -- ",
      "using the %s fallback instead."), predictor_name, best_tier, min_cases, tier),
      call. = FALSE)
  }

  list(fit = fit, marginal_fit = marginal_fit, tier = tier,
    design_names = if (identical(tier, "marginal")) character(0) else
      colnames(if (identical(tier, "full")) full_design else others))
}

# The conditional mean for one subject's missing cell, from an already-fit
# `.ctJuliaTIImputationFit()` result. A subject missing one of the fitted
# model's own inputs (e.g. a second TI predictor that is *also* missing for
# them) cannot use it -- rather than fabricate that input, this falls back
# to the marginal conditional for this subject alone, the same dropout
# principle applied throughout this file.
.ctJuliaTIPredictMu <- function(imputation, newdata_row) {
  if (identical(imputation$tier, "marginal")) {
    return(unname(stats::coef(imputation$fit)[["(Intercept)"]]))
  }
  vals <- newdata_row[imputation$design_names]
  if (anyNA(vals)) {
    return(unname(stats::coef(imputation$marginal_fit)[["(Intercept)"]]))
  }
  nd <- as.data.frame(as.list(vals))
  unname(as.numeric(stats::predict(imputation$fit, newdata = nd)))
}

# The sampling path's imputation spec: one row per missing TI predictor cell,
# with a pre-computed conditional mean and SD (`parameter`, the raw index
# that will sample it, is assigned later in `.ctJuliaPrepare`, once the rest
# of the raw vector's length is known). Empty when nothing is missing.
.ctJuliaTIMissingSpec <- function(dat, model, raw_tipred, includeOutcome = TRUE) {
  missing_mask <- is.na(raw_tipred)
  empty <- data.frame(subject = integer(), predictor = integer(),
    mu = numeric(), sigma = numeric())
  if (!any(missing_mask)) return(empty)

  manifest_summary <- .ctJuliaTISubjectManifestSummary(dat, model)
  rows <- list()
  for (j in seq_len(model$n.TIpred)) {
    if (!any(missing_mask[, j])) next
    others <- if (model$n.TIpred > 1L) raw_tipred[, -j, drop = FALSE] else
      matrix(numeric(), nrow(raw_tipred), 0L)
    imputation <- .ctJuliaTIImputationFit(model$TIpredNames[j], raw_tipred[, j],
      others, manifest_summary$ybar, manifest_summary$ysd, includeOutcome)
    sigma <- stats::sigma(imputation$fit)
    for (i in which(missing_mask[, j])) {
      newdata_row <- c(others[i, , drop = TRUE], manifest_summary$ybar[i, , drop = TRUE],
        manifest_summary$ysd[i, , drop = TRUE])
      mu <- .ctJuliaTIPredictMu(imputation, newdata_row)
      rows[[length(rows) + 1L]] <- data.frame(subject = as.integer(i),
        predictor = as.integer(j), mu = as.numeric(mu), sigma = as.numeric(sigma))
    }
  }
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
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

# `offset` is where the TI-predictor coefficients start in the raw vector. It
# defaults to the last parameter the table itself uses, which is right for the
# augmented route because that route puts its population-covariance parameters
# *in* the table. The Laplace route's population parameters are not in any
# matrix cell, so it passes the end of its own block instead; without that the
# coefficients would silently alias the population scales.
.ctJuliaTIEffects <- function(table, model, offset = NULL) {
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
  # The leading zero is the count of a model with nothing free, and it is not
  # hypothetical: `ctGenerate` resolves every free parameter to a value before
  # preparing, so on the generation path every `parnumber` is NA and `max` had
  # nothing to take a maximum over -- it warned and returned -Inf, on an
  # ordinary call. Same guard as `.ctBackendNpar()` and `.ctGenerateJulia()`
  # already use for the same reason.
  coefficient <- if (is.null(offset)) max(c(0L, table$parnumber), na.rm = TRUE) else
    as.integer(offset)
  # Coefficient index by `.ctTipredEffectKey()`: two parameters given the same
  # effect name share one. The stan path keeps the same map in
  # R/ctModelWriter.R. The indices themselves differ between the two -- stan
  # allocates as it walks the cells, this walks predictors then parameters --
  # but which effects are *the same coefficient* has to match, or the two
  # backends fit different models from one specification.
  keys <- list()
  for (predictor in seq_along(model$TIpredNames)) {
    column <- effect_columns[predictor]
    if (!column %in% available) next
    # `Free`, not `Active`: an effect fixed to a value carries the value and
    # needs no coefficient, so counting it would claim a parameter that is not
    # estimated. R/ctModelWriter.R has always used `Free` here; this used
    # `Active` and was inert only because nothing reached it with a fixed
    # effect -- measured, on the fitting path `.ctJuliaUnsupported()` refuses
    # a fixed effect before this runs, and on the generation path every
    # `parnumber` is NA so `direct` excludes every row and this returns
    # nothing. Both of those are other functions' behaviour, and the julia
    # refusal is worded "v1 does not support", so the day it does the wrong
    # test here would silently free an effect the user fixed.
    freeeffects <- direct & .ctTipredEffectFree(table[[column]])
    parameters <- sort(unique(table$parnumber[freeeffects]))
    for (parameter in parameters) {
      # Cells sharing a parameter are required to state the same effects
      # (R/ctModelSpecCheck.R), so the first of them carries the label.
      first <- which(freeeffects & table$parnumber %in% parameter)[1L]
      key <- .ctTipredEffectKey(table[[column]][first], predictor, parameter)
      if (is.null(keys[[key]])) {
        coefficient <- coefficient + 1L
        keys[[key]] <- coefficient
      }
      entries[[length(entries) + 1L]] <- data.frame(
        parameter = as.integer(parameter), predictor = as.integer(predictor),
        coefficient = as.integer(keys[[key]])
      )
    }
  }
  if (!length(entries)) return(data.frame(parameter = integer(), predictor = integer(), coefficient = integer()))
  do.call(rbind, entries)
}

# `intoverpop` deliberately has no default. It selects which *model* is
# prepared -- random effects as latent states, or integrated by Laplace -- and a
# default meant a caller could omit it and silently get the other one. That is
# what happened to prediction and to cross-validation, both of which rebuild a
# specification from a fit and neither of which said which kind of fit it was.
# Making it mandatory turns that from a silent wrong answer into a stop at the
# call site.
.ctJuliaPrepare <- function(datalong, model, prepared_data = NULL, project = NULL,
  priors = FALSE, intoverpop, optimize = TRUE, tipredMissingIncludeOutcome = TRUE) {
  # "none" prepares exactly as "laplace" does. The Laplace specification is what
  # *describes* the random effects -- which raw parameters vary, at which level,
  # with which population scale -- and that description is needed whether they
  # are integrated out or sampled. The two differ only in what the fit then does
  # with it, so they must not differ in how the model is built.
  intoverpop <- match.arg(as.character(intoverpop)[1L],
    c("augmented", "laplace", "none"))
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
  # A missing TI predictor cell, on the sampling path only. `optimize=TRUE`
  # is untouched here: `.ctPrepareData()` (`prepared_data`, when supplied)
  # already regression-imputed it into a real number before this function
  # ever saw it (`ctData.R`; that is what Charles asked for and it is
  # cheaper), so nothing below runs and `tipred_data` is used exactly as
  # prepared -- as it always was.
  #
  # `optimize=FALSE` is the new path, SPEC-tipred-sampling.md, and missingness
  # is detected from the raw data (`dat`) rather than from whatever
  # `prepared_data$tipredsdata` holds, because on this path Stan would have
  # written its own 99999 sentinel there and this engine does not read that
  # convention (see `.ctJuliaTIMissingPlaceholder` for the one it does use,
  # and why it is not the same thing). Only `intoverpop='augmented'` is
  # supported yet: sampling a missing predictor value needs
  # `gradient_method='forward'`, which the Laplace inner solve and the joint
  # sampler do not (yet) offer an equivalent for.
  ti_missing <- NULL
  if (!isTRUE(optimize) && model$n.TIpred > 0) {
    raw_tipred <- .ctJuliaTIDataAllowMissing(dat, model)
    missing_mask <- is.na(raw_tipred)
    n_missing <- sum(missing_mask)
    if (n_missing > 0L) {
      if (identical(intoverpop, "augmented")) {
        ti_missing <- .ctJuliaTIMissingSpec(dat, model, raw_tipred,
          includeOutcome = tipredMissingIncludeOutcome)
        tipred_data[missing_mask] <- .ctJuliaTIMissingPlaceholder
      } else {
        stop("Julia backend cannot sample missing TI predictor values (",
          n_missing, " missing) with intoverpop='", intoverpop, "'; it can ",
          "currently do so only with intoverpop='augmented'. Impute them ",
          "before fitting, drop those subjects, use intoverpop='augmented', ",
          "or use backend='stan'.", call. = FALSE)
      }
    }
  }
  max_timestep <- if (!is.null(prepared_data$maxtimestep)) {
    as.numeric(prepared_data$maxtimestep)[1L]
  } else if (!is.null(model$nlcontrol$maxtimestep)) {
    as.numeric(model$nlcontrol$maxtimestep)[1L]
  } else 999999
  if (!is.finite(max_timestep) || max_timestep <= 0) stop("Julia maxtimestep must be a positive finite number.", call. = FALSE)
  # The automatic substep policy, when asked for. The mesh itself is decided
  # later by `.ctJuliaAutoSubsteps` at the starting values and replaces
  # `max_timestep` on the spec, so every objective built from the spec -- the
  # fit, prediction, the state-explicit target -- integrates the same way.
  # `floor` is the maxtimestep rule the mesh never goes below.
  substeps <- if (identical(model$nlcontrol$nsubsteps, "auto")) list(
    tol = as.numeric(.ctJuliaOr(model$nlcontrol$substeptol, 0.01))[1L],
    max = as.integer(.ctJuliaOr(model$nlcontrol$maxsubsteps, 64L))[1L],
    floor = max_timestep) else NULL
  # The state-explicit path's step between substeps; see `_ctsem_transition`.
  transition <- as.character(.ctJuliaOr(model$nlcontrol$transition, "exponential"))[1L]
  laplace <- NULL
  if (intoverpop %in% c("laplace", "none")) {
    # No state augmentation at all: the model the engine filters is the plain
    # per-subject one, and the random effects are described alongside it.
    parameter_table <- .ctJuliaParameterTable(model)
    laplace <- .ctJuliaLaplaceSpec(model, parameter_table, prepared_data, dat)
    if (!laplace$nrandom) {
      stop("intoverpop='", intoverpop, "' was requested but no parameters are marked ",
        "indvarying, so there is nothing to integrate over. Mark parameters as ",
        "varying in the model, or leave intoverpop at its default.", call. = FALSE)
    }
    augmented <- list(nlatent = model$n.latent, nlatent_augmented = model$n.latent,
      dynamic_state_indices = .ctJuliaDerrind(parameter_table, model$n.latent, model$n.latent),
      random_effects = data.frame(),
      rewritten_cells = .ctJuliaRewrittenCells(parameter_table))
    ti_effects <- .ctJuliaTIEffects(parameter_table, model, offset = laplace$npar)
  } else {
    augmented <- .ctJuliaAugmentRandomEffects(model)
    parameter_table <- augmented$parameter_table
    ti_effects <- .ctJuliaTIEffects(parameter_table, model)
  }
  # `laplace$npar` already counts every level's scales and correlations. Taking
  # the maximum over the *level-one* index vectors instead sized the raw vector
  # to the subject level alone, so an outer level's scale sat past the end of
  # it: every evaluation threw a bounds error, was swallowed by the optimizer's
  # invalid-point guard, and the level stayed pinned at its starting value while
  # reporting a plausible-looking number.
  # `laplace$npar` counts the model parameters and every level's scales and
  # correlations; the TI-predictor coefficients sit after all of that, so the
  # raw vector is as long as whichever reaches furthest. Taking `laplace$npar`
  # alone left the coefficients past the end of it -- the same failure the
  # level scales had, one block further along.
  # The leading zero is the count of a model with nothing free. Every index
  # vector here is NA-filled for a fixed cell, so a fully fixed model -- what
  # `ctGenerate` prepares, having resolved every free parameter to a value --
  # leaves `max` nothing to take a maximum over: it warns and returns -Inf.
  # The pieces are still loose here -- this is where the spec is assembled --
  # so they are handed over as one rather than the count being written out a
  # sixth time.
  # Sampled TI-predictor values get the last block of the raw vector, one
  # new index per missing cell, past everything above -- the same "reserve N
  # indices past the current end" move `.ctJuliaTIEffects()`'s `offset` makes
  # for TI-effect coefficients, and for the same reason: nothing here is a
  # matrix cell, so nothing places it automatically.
  if (!is.null(ti_missing) && nrow(ti_missing)) {
    offset <- .ctBackendNpar(list(parameter_table = parameter_table,
      laplace = laplace, ti_effects = ti_effects))
    ti_missing$parameter <- offset + seq_len(nrow(ti_missing))
  }
  npar <- .ctBackendNpar(list(parameter_table = parameter_table,
    laplace = laplace, ti_effects = ti_effects, ti_missing = ti_missing))
  .ctJuliaCheckLayout(parameter_table, laplace, ti_effects, npar, ti_missing = ti_missing)
  prior_spec <- if (!isTRUE(priors)) NULL else if (!is.null(laplace)) {
    .ctBackendLaplacePriorSpec(prepared_data, laplace, npar)
  } else .ctBackendPriorSpec(prepared_data, npar)
  list(
    class = "ctJuliaModel",
    intoverpop = intoverpop,
    laplace = laplace,
    model = model,
    data = dat,
    parameter_table = parameter_table,
    subject_starts = as.integer(subject_starts),
    times = as.numeric(dat[[model$timeName]]),
    # Doubles with NaN for missing, never a union type.
    #
    # `as.matrix` on a data frame of integers with NAs marshals to
    # `Matrix{Union{Missing,Int64}}`, and on doubles with NAs to
    # `Matrix{Union{Missing,Float64}}` -- a *different element type per
    # dataset*, and none of them the `Matrix{Float64}` the precompile workload
    # built. The filter is specialised on that type, so every real fit
    # recompiled it however well the model shape matched: in a pure Julia
    # session the captured shape's first evaluate took 0.017 s, through R it
    # took twelve seconds.
    #
    # `_ctsem_observed` is `!ismissing(x) && isfinite(x)`, so NaN already means
    # missing to the engine and nothing is lost by saying it that way. A union
    # element type is also boxed on every access in the hot loop.
    #
    # This takes the first evaluation from 12.8 s to 3.5 s. It does not close the
    # gap: the same shape in a pure Julia session is 0.017 s, so about 3.5 s of
    # shape compilation survives. That remainder has a separate cause -- the
    # captured shapes hold unsimplified transform strings (`0 + 10 * (param[1] *
    # 1 + 0)`) where the fit path sends simplified ones (`10 * param[1]`), so the
    # cached closure types differ from the ones in the image and the pipeline
    # specialises again. The fix belongs in tools/generate-precompile-shapes.R.
    manifest_data = {
      d <- t(as.matrix(dat[, model$manifestNames, drop = FALSE]))
      storage.mode(d) <- "double"
      d[is.na(d)] <- NaN
      d
    },
    tdpred_data = t(tdpred_data),
    tipred_data = as.matrix(tipred_data),
    ti_effects = ti_effects,
    # One row per sampled (missing) TI-predictor cell: `subject`, `predictor`
    # (both 1-based, matching `tipred_data`'s rows/columns), the raw-parameter
    # `parameter` that samples it, and its pre-computed conditional `mu`/
    # `sigma`. `NULL` for every model with nothing missing.
    ti_missing = ti_missing,
    priors = prior_spec,
    max_timestep = max_timestep,
    substeps = substeps,
    transition = transition,
    # A discrete-time model is the same filter with a different discretization:
    # DRIFT, CINT and DIFFUSION are already the one-step quantities, so the
    # exponential, the Lyapunov solve and the intercept solve all collapse.
    continuoustime = isTRUE(model$continuoustime),
    # Which covariance construction to use, as the integer code both backends
    # share -- 0 for the unconstrained correlation square root, 2 for
    # covmattransform='z'. Carried on the spec so a fit rebuilt from a saved
    # object constructs its covariances the way it was estimated.
    covmatcode = .ctCovMatCode(model),
    TDpredNames = model$TDpredNames,
    TIpredNames = model$TIpredNames,
    # 0 Gaussian, 1 binary, 2 ordinal, with the category count alongside for
    # the ordinal ones. Carried on the spec so a fit rebuilt from a saved
    # object knows its own measurement model.
    manifesttype = if (is.null(model$manifesttype)) integer(0) else
      as.integer(model$manifesttype),
    ncategories = if (is.null(model$ncategories)) integer(0) else
      as.integer(model$ncategories),
    censormin = if (is.null(model$censormin)) numeric(0) else
      as.numeric(model$censormin),
    censormax = if (is.null(model$censormax)) numeric(0) else
      as.numeric(model$censormax),
    nlatent = augmented$nlatent,
    nlatent_augmented = augmented$nlatent_augmented,
    dynamic_state_indices = augmented$dynamic_state_indices,
    random_effects = augmented$random_effects,
    rewritten_cells = augmented$rewritten_cells,
    project = project,
    engine = .ctJuliaEngineVersion()
  )
}

# Say so, once, when a model shape is about to be compiled for.
#
# `ctsem_transforms_cached()` reports whether every transform expression in this
# model already has a closure. If they all do, the compiled filter is reused and
# the first evaluation is immediate; if any is new, the specialisation happens
# on the first evaluation and takes tens of seconds.
.ctJuliaAnnounceCompilation <- function(module, table) {
  free <- .ctJuliaNoNA(as.integer(table$parnumber), 0L) > 0L
  # A free cell with no transform of its own compiles `param[n]`, exactly as
  # `ekf_from_columns` writes it -- so the query has to say the same thing.
  regular <- .ctJuliaNoNA(as.character(table$transform), "")[free]
  numbers <- .ctJuliaNoNA(as.integer(table$parnumber), 0L)[free]
  regular <- ifelse(nzchar(regular), regular, paste0("param[", numbers, "]"))
  complex <- c(.ctJuliaNoNA(as.character(table$predicttransform), ""),
    .ctJuliaNoNA(as.character(table$updatetransform), ""),
    .ctJuliaNoNA(as.character(table$tdtransform), ""))
  regular <- unique(regular[nzchar(regular)])
  complex <- unique(complex[nzchar(complex)])
  if (!length(regular) && !length(complex)) return(invisible(FALSE))
  cached <- isTRUE(tryCatch(
    .ctBackendJuliaValue(module$ctsem_transforms_cached(
      .ctJuliaVector(if (length(regular)) regular else ""),
      .ctJuliaVector(if (length(complex)) complex else ""))),
    error = function(e) FALSE))
  if (cached) return(invisible(FALSE))
  # No duration promised. It ranges from a few seconds when the precompiled
  # image covers the shape to a couple of minutes when it does not, and a
  # message that says "20-60 seconds" is wrong at both ends -- alarming for the
  # quick case and misleading for the slow one.
  message("Compiling the julia engine for this model shape (once per shape ",
    "per session).")
  invisible(TRUE)
}

# How many raw parameters a model has.
#
# The parameter table alone undercounts, and that has already caused a bug: the
# random-effect population scales sit above it in the raw vector and the
# TI-predictor coefficients above those, so `max(parnumber)` gives a number that
# looks plausible and is short by however many random effects the model has. The
# engine then refuses the vector with "level 1 references population scale 17
# but the parameter vector holds 13".
#
# Written out in five places before this, in three variants: some guarded a
# fully fixed model with a leading `0L` (without which `max` warns and returns
# -Inf), some omitted it and tested `is.finite` afterwards instead, and only one
# unwrapped `model_spec` for callers holding a fit rather than a spec. Variants
# that differ are worse than plain duplication, because the difference is
# invisible at each site.
#
# Accepts a spec, a prepared model, or a fit.
.ctBackendNpar <- function(x) {
  spec <- if (!is.null(x$model_spec)) x$model_spec else x
  n <- suppressWarnings(max(c(0L, spec$parameter_table$parnumber,
    spec$laplace$npar, spec$ti_effects$coefficient, spec$ti_missing$parameter),
    na.rm = TRUE))
  if (!is.finite(n)) 0L else as.integer(n)
}

# The integer code for a model's covariance construction, shared with the stan
# path's `standata$choleskymats` so one model means one construction on either
# backend. 'rawcorr_indep' and 'cholesky' are refused above rather than mapped:
# they were accepted and then ignored here once already.
.ctCovMatCode <- function(model) {
  tf <- if (is.null(model$covmattransform)) "rawcorr" else
    as.character(model$covmattransform)[1L]
  if (identical(tf, "z")) 2L else 0L
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
  # Omitted when zero: the engine defaults to it, and this keeps a spec built
  # before the setting existed working unchanged.
  covmatcode <- if (is.null(spec$covmatcode)) 0L else as.integer(spec$covmatcode)
  if (covmatcode != 0L) arguments$covmatcode <- covmatcode
  # Only when something is actually binary: an empty vector lets the engine
  # skip the branch, and a zero-length vector deadlocks the bridge, so the two
  # reasons to omit it agree.
  if (any(spec$manifesttype > 0)) {
    arguments$manifesttype <- .ctJuliaVector(as.integer(spec$manifesttype))
  }
  if (any(spec$manifesttype %in% 2)) {
    arguments$ncategories <- .ctJuliaVector(as.integer(spec$ncategories))
  }
  if (any(spec$manifesttype %in% 4)) {
    arguments$censormin <- .ctJuliaVector(as.numeric(spec$censormin))
    arguments$censormax <- .ctJuliaVector(as.numeric(spec$censormax))
  }
  # A model shape Julia has not seen mints new closure types for its transform
  # expressions, and the whole filter specialises again for them -- tens of
  # seconds, once, and indistinguishable from a hang if nothing says so. The
  # engine is asked first so the message comes before the wait rather than
  # after it.
  .ctJuliaAnnounceCompilation(module, table)
  params <- do.call(module$ekf_from_columns, arguments)
  # .ctJuliaVector, not juliaPut, for the vectors: JuliaConnectoR marshals a
  # length-one R vector as a *scalar*, so a single-subject model would hand the
  # objective constructor an Int where it wants an AbstractVector. That never
  # showed up while every caller was a whole fitted dataset; ctPredict() on one
  # subject is a caller where it does.
  objective_args <- list(params, .ctJuliaVector(spec$subject_starts),
    .ctJuliaVector(spec$times), JuliaConnectoR::juliaPut(spec$manifest_data),
    JuliaConnectoR::juliaPut(spec$tdpred_data), JuliaConnectoR::juliaPut(spec$tipred_data),
    .ctJuliaSubstepArgument(spec$max_timestep))
  if (!is.null(spec$priors) && length(spec$priors$index)) {
    objective_args$prior_index <- .ctJuliaVector(spec$priors$index)
    objective_args$prior_scale <- .ctJuliaVector(spec$priors$scale)
    objective_args$prior_weight <- spec$priors$weight
  }
  # Sampled TI-predictor values, one entry per missing cell -- omitted
  # entirely (not sent as empty vectors) for every model with none, which is
  # what keeps that the zero-cost path on the Julia side too.
  if (!is.null(spec$ti_missing) && nrow(spec$ti_missing)) {
    objective_args$ti_missing_subject <- .ctJuliaVector(as.integer(spec$ti_missing$subject))
    objective_args$ti_missing_predictor <- .ctJuliaVector(as.integer(spec$ti_missing$predictor))
    objective_args$ti_missing_parameter <- .ctJuliaVector(as.integer(spec$ti_missing$parameter))
    objective_args$ti_missing_mu <- .ctJuliaVector(as.numeric(spec$ti_missing$mu))
    objective_args$ti_missing_sigma <- .ctJuliaVector(as.numeric(spec$ti_missing$sigma))
  }
  objective <- do.call(module$ctsem_objective, objective_args)
  # The Laplace route wraps the ordinary objective rather than replacing it:
  # the process likelihood, the parameter layer and the TI-predictor effects
  # are all still the same code, evaluated per subject at a shifted parameter
  # vector. Only zero-length index vectors are withheld, because JuliaConnectoR
  # deadlocks marshalling one -- `cor_index` is empty whenever a model has a
  # single random effect.
  if (!is.null(spec$laplace)) {
    levels <- spec$laplace$levels
    grab <- function(field) unlist(lapply(levels, function(x) x[[field]]), use.names = FALSE)
    laplace_args <- list(objective,
      re_index = .ctJuliaVector(as.integer(grab("re_index"))),
      sd_index = .ctJuliaVector(as.integer(grab("sd_index"))),
      sd_scale = .ctJuliaVector(as.numeric(grab("sd_scale"))))
    if (length(grab("cor_index"))) {
      laplace_args$cor_index <- .ctJuliaVector(as.integer(grab("cor_index")))
    }
    if (length(levels) > 1L) {
      # Concatenated innermost level first, split on the far side by the
      # per-level counts. Flat vectors because the bridge marshals those and
      # not nested structures.
      laplace_args$level_nre <- .ctJuliaVector(as.integer(vapply(levels,
        function(x) x$nrandom, integer(1))))
      laplace_args$group <- .ctJuliaVector(as.integer(grab("group")))
      laplace_args$level_ngroups <- .ctJuliaVector(as.integer(vapply(levels,
        function(x) x$ngroups, integer(1))))
    }
    objective <- do.call(module$ctsem_laplace_objective, laplace_args)
  }
  assign(key, objective, envir = .ct_julia_cache$objectives)
  objective
}

# The state-explicit objective ------------------------------------------------
#
# `intoverstates=FALSE`: the target is the joint density of the parameters and
# the innovations that build the latent states, over `x = [theta; z]`. The
# engine's optimiser and sampler both take it exactly as they take the marginal
# one -- they ask for a value and a gradient at a vector, and this answers --
# so nothing about the fit's control flow changes, only which objective it is
# handed.
#
# Cached beside the marginal objective and keyed on the same spec plus the
# parameter count, because building it walks the design to count innovations.

.ctJuliaJointObjective <- function(object, npar) {
  # A fit, a classed model, or the bare spec `.ctJuliaPrepare` returns -- the
  # fit path holds the last of those and classing it at every call site is how
  # a caller ends up passing the wrong one.
  spec <- if (inherits(object, "ctJuliaFit")) object$model_spec else object
  if (!inherits(spec, "ctJuliaModel")) {
    if (!is.list(spec) || is.null(spec$parameter_table)) {
      stop("object must be a ctJuliaModel, a ctJuliaFit, or a prepared spec.",
        call. = FALSE)
    }
    spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  }
  if (!is.null(spec$laplace)) {
    stop("intoverstates=FALSE cannot be combined with intoverpop='laplace' or ",
      "'none': the Laplace route wraps the likelihood in a per-subject inner ",
      "problem over the random effects, and the state path replaces the ",
      "likelihood itself. Use the default intoverpop, which carries random ",
      "effects as augmented states and so is sampled along with them.",
      call. = FALSE)
  }
  npar <- max(1L, as.integer(npar)[1L])
  transition <- .ctJuliaOr(spec$transition, "exponential")
  key <- paste0(.ctJuliaObjectiveKey(spec), "|joint|", npar, "|", transition)
  if (exists(key, envir = .ct_julia_cache$objectives, inherits = FALSE)) {
    return(get(key, envir = .ct_julia_cache$objectives, inherits = FALSE))
  }
  module <- .ctJuliaModule(spec$project)
  objective <- module$ctsem_joint_objective(.ctJuliaObjective(spec), npar,
    transition = transition)
  assign(key, objective, envir = .ct_julia_cache$objectives)
  objective
}

# How many innovations this design needs, which is what the parameter vector is
# extended by.
.ctJuliaStateDimension <- function(spec) {
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  .ctBackendStateDimension(handle)
}

# The latent trajectory an estimate implies, as rows by latents -- the thing a
# state-explicit fit has that a marginal one has to run a smoother for.
.ctJuliaJointStates <- function(spec, joint, x) {
  module <- .ctJuliaModule(spec$project)
  states <- .ctBackendJuliaValue(module$ctsem_joint_states(joint,
    .ctJuliaNumericVector(as.numeric(x))))
  nlatent <- length(spec$model$latentNames)
  out <- t(matrix(as.numeric(states), nrow = nlatent))
  colnames(out) <- spec$model$latentNames
  out
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
summary.ctJuliaFit <- function(object, timeinterval = 1, digits = 3, parmatrices = TRUE,
  priorcheck = TRUE, residualcov = TRUE, ...) {
  # `priorcheck` is accepted and ignored rather than rejected: it is part of
  # summary.ctStanFit's signature, and a script that summarises whichever fit it
  # was handed should not fail on the argument. What it reports -- posterior
  # means and sds against ctsem's normal(0,1) raw priors -- would need the Stan
  # model's own prior block, which this backend does not carry.
  .ctBackendSummary(object, timeinterval = timeinterval, digits = digits,
    parmatrices = parmatrices, residualcov = residualcov, ...)
}

#' @export
ctExtract.ctJuliaFit <- function(object, subjectMatrices = FALSE, cores = 2,
  nsamples = "all", subjects = "all", state = NULL, ...) {
  # `cores` was accepted and dropped. It is the engine's subject-chunk ceiling
  # here, not a number of R processes -- there is no cluster on this path -- and
  # it is restored afterwards so an extract does not leave the session
  # reconfigured. Note that `subjectMatrices=TRUE` is dominated by moving the
  # filter output back across the bridge rather than by computing it, so this
  # bounds the work rather than speeding it up much.
  #
  # `state` reached the engine raw before, so the shorthands 'mean' and
  # 'asymptotic' arrived as a string and died in the engine with a Julia
  # MethodError; only an explicit vector worked, and it came back unlabelled.
  # Resolve it here, as ctSummaryMatrices does, so both spell the same point the
  # same way.
  resolved <- .ctResolveState(object, state)
  out <- .ctBackendWithMaxChunks(cores,
    .ctBackendExtract(object, subjectMatrices = subjectMatrices, nsamples = nsamples,
      subjects = subjects, state = resolved$state, ...))
  # Recorded, not messaged: see ctExtract.ctStanFit.
  attr(out, "evaluatedAt") <- resolved$label
  out
}

#' @export
ctSummaryMatrices.ctJuliaFit <- function(fit, calcfunc = quantile,
  calcfuncargs = list(probs = 0.5), timeinterval = 1, state = NULL, ...) {
  # 'T0MEANS' (the default), 'mean', 'asymptotic', or a state vector. Resolved
  # here rather than passed on raw so that the message below can name the point
  # in the same words whichever form the caller used.
  resolved <- .ctResolveState(fit, state)
  out <- .ctBackendSummaryMatrices(fit, calcfunc = calcfunc, calcfuncargs = calcfuncargs,
    timeinterval = timeinterval, state = resolved$state, ...)
  .ctContextMessage(fit, resolved$label)
  .ctContextAttach(out, fit)
}

# Run `expr` with the engine's subject-chunk ceiling set to `chunks`, and put
# the previous ceiling back afterwards.
#
# `_CTSEM_MAX_CHUNKS` is session-global Julia state, and three call sites wrote
# it and none restored it: after a `cores=8` fit the session read 8, after a
# `cores=1` fit it read 1, and every later ctKalman(), ctJuliaEvaluate(),
# ctExtract() or ctLOO() in that session inherited whichever fit came last. That
# is performance-only, but it makes a timing depend on history, which is exactly
# what makes one impossible to reproduce.
.ctBackendWithMaxChunks <- function(chunks, expr) {
  previous <- .ctBackendSetMaxChunks(chunks)
  on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  expr
}

# Set the ceiling, and report the one that was there.
#
# Three places cap the ceiling for the duration of some work and put it back:
# the wrapper above, the uncertainty phase, and ctLaplaceCheck. Each had its own
# copy of read-set-restore and the copies had drifted -- one set only for
# `chunks >= 1`, one only for `chunks > 1`, one for any non-NULL `cores`; two
# wrapped the set in `try` and the third let a failure abort the caller. The
# restore half stays the caller's own `on.exit`, because two of the three cap
# for the rest of their function rather than for a single expression.
#
# NA back means the previous ceiling could not be read, and so that there is
# nothing to put back. It is read only when there is a ceiling to set, because
# reading it is a `juliaEval` and `juliaEval` *starts* a session when none is
# running -- capping nothing must not be the thing that launches the engine.
.ctBackendSetMaxChunks <- function(chunks) {
  chunks <- suppressWarnings(as.integer(chunks)[1L])
  if (is.na(chunks) || chunks < 1L) return(NA_integer_)
  previous <- tryCatch(as.integer(.ctBackendJuliaValue(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM.ctsem_max_chunks().max_chunks"))), error = function(e) NA_integer_)
  try(JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!", chunks),
    silent = TRUE)
  previous
}

.ctBackendRestoreMaxChunks <- function(previous) {
  previous <- suppressWarnings(as.integer(previous)[1L])
  if (!is.na(previous)) {
    try(JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!",
      previous), silent = TRUE)
  }
  invisible(NULL)
}

# The Julia session's thread count, or NA -- without starting a session.
#
# `juliaEval()` starts a server rather than reporting on one, so a running
# session is asked directly and a session that has not begun is predicted from
# the environment instead. Julia also accepts "auto" and "N,M" (worker threads,
# interactive threads) there; "auto" is never short of anything, and any other
# spelling this cannot parse is reported as unknown rather than guessed at.
.ctBackendSessionThreads <- function() {
  if (.ctJuliaSessionRunning()) {
    return(tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
      error = function(e) NA_integer_))
  }
  existing <- Sys.getenv("JULIA_NUM_THREADS", unset = "")
  if (!nzchar(existing)) return(1L)
  suppressWarnings(as.integer(sub(",.*$", "", existing))[1L])
}

# Say so when `cores` is more than the session can give, once per session per
# pair. Separate from `.ctBackendReportChunks()` below on purpose: that one is
# about the chunk tuner declining headroom it measured, and this is about
# headroom that never existed, which is a different thing to say and a different
# thing to do about it.
.ctBackendReportThreads <- function(cores, threads) {
  key <- paste0(cores, ":", threads)
  if (key %in% .ct_julia_cache$threads_reported) return(invisible(NULL))
  .ct_julia_cache$threads_reported <- c(.ct_julia_cache$threads_reported, key)
  message("cores = ", cores, " requested, ", threads, " used: Julia's thread ",
    "count is fixed at session start. ctJuliaSetup(threads = ", cores,
    ", force = TRUE) restarts it with ", cores, ".")
  invisible(NULL)
}

# Resolve `cores` against that thread count, and never do it in silence.
#
# Julia fixes `Threads.nthreads()` at process start, so `cores` is a request
# only an unstarted session can grant. Three situations, and the third is the
# one this exists for:
#
#   * No session yet and the thread count is ours to set -- set it, and the fit
#     gets what it asked for.
#   * No session yet but JULIA_NUM_THREADS was set deliberately, by
#     ctJuliaSetup(threads=), by the user, or by a cluster scheduler -- that
#     wins, and it may be narrower than `cores`.
#   * A session is already running -- its count is fixed and cannot be raised.
#
# The third cost a whole benchmark pass. A harness that calls `ctJuliaSetup()`
# once up front, to pay the engine load before timing anything, pins the session
# to one thread; every later `cores = n` fit then ran on one subject chunk,
# serially, and said nothing. It was caught only because the gradient counts for
# `cores = 1` and `cores = 4` came back identical to the digit.
#
# `options(ctsem.julia.restart = TRUE)` opts in to fixing it rather than saying
# it. Not the default, for two measured reasons: a restart discards the
# session's compiled model shapes -- ~5 s to restart, plus ~15 s of
# respecialisation on the next fit, on a one-latent model -- and `cores`
# defaults to `getOption("mc.cores", 2)`, so a default fit would otherwise kill
# and rebuild a deliberately narrow session nobody asked it to widen.
#
# `threads` is a parameter so the gate can be exercised at a count the test
# machine does not have to be restarted into. `report = FALSE` is `fit = FALSE`:
# a call that only builds a specification still wants the thread count requested
# while that is possible, since preparing starts the session, but has nothing to
# say about cores it was never going to use and no reason to restart anything.
.ctBackendResolveThreads <- function(cores, threads = NULL, report = TRUE) {
  cores <- suppressWarnings(as.integer(cores)[1L])
  if (is.na(cores) || cores < 2L) return(invisible(NA_integer_))
  restart <- isTRUE(report) && isTRUE(getOption("ctsem.julia.restart", FALSE))
  if (is.null(threads) && !.ctJuliaSessionRunning()) {
    existing <- Sys.getenv("JULIA_NUM_THREADS", unset = "")
    # An existing value is respected unless *this* function set it for an
    # earlier fit. Without that distinction a `cores=8` fit left the variable
    # behind, and the next `ctFit(cores=2)` in a restarted session started
    # eight threads while asking for two -- the subject loop still honoured
    # `cores`, but the process held cores the user had not asked for. A value
    # from ctJuliaSetup(threads=) or from the user's own environment is
    # deliberate and still wins, unless the restart option says otherwise --
    # there is no session to lose here, so honouring it costs nothing.
    ours <- nzchar(existing) &&
      identical(existing, .ct_julia_cache$threads_from_cores)
    if (!nzchar(existing) || ours || restart) {
      Sys.setenv(JULIA_NUM_THREADS = as.character(cores))
      .ct_julia_cache$threads_from_cores <- as.character(cores)
      return(invisible(cores))
    }
  }
  if (is.null(threads)) threads <- .ctBackendSessionThreads()
  threads <- suppressWarnings(as.integer(threads)[1L])
  if (is.na(threads) || threads >= cores || !isTRUE(report)) {
    return(invisible(threads))
  }
  if (restart) {
    # Said before it happens, because it takes seconds and ends a process the
    # caller may have set up. The engine environment is ctJuliaSetup()'s own
    # `project`, and .ctJuliaClearSession() forgets it, so it is carried over.
    message("Restarting the Julia session at ", cores, " threads for cores = ",
      cores, ".")
    project <- .ct_julia_cache$project
    got <- tryCatch({
      suppressWarnings(ctJuliaSetup(project = project, threads = cores,
        force = TRUE))
      .ct_julia_cache$threads_from_cores <- as.character(cores)
      .ctBackendSessionThreads()
    }, error = function(e) NA_integer_)
    if (!is.na(got)) {
      if (got >= cores) return(invisible(got))
      threads <- got
    }
  }
  .ctBackendReportThreads(cores, threads)
  invisible(threads)
}

# Say so when the chunk tuner used materially fewer chunks than `cores` allowed.
#
# `cores` is a ceiling, not an instruction: `ctsem_tune_chunks!` times a ladder
# of chunk counts within it and pins the fastest, because the subject loop is
# not monotone in the count -- one model measured 3.7x *slower* on 23 threads
# than on one. That is the right behaviour, and it was entirely silent, so
# `ctFit(cores = 12)` could run on two chunks with nothing said: no way to tell
# that ten cores went unused, and no way to tell whether asking for more would
# have helped or hurt.
#
# Two things keep it from becoming noise. It fires only where the tuner itself
# left the headroom -- the ceiling is `min(cores, threads)`, so a session with
# fewer threads than `cores` says nothing here rather than blaming the tuner
# for a limit it never saw -- and each limit-and-count pair is said once per
# session, so a simulation study looping a hundred fits gets one line.
#
# That first exclusion left the worst case silent for a while: a `cores = 4` fit
# in a one-thread session is not the tuner declining headroom, it is headroom
# that never existed, and nothing said so. The exclusion is still right -- the
# tuner is not at fault and naming it would misattribute -- so the case is
# reported by `.ctBackendResolveThreads()` above instead, before the fit and
# against `threads`, which is where the fix is.
#
# `threads` is the session's thread count and is asked for when not supplied.
# It is a parameter so the gate can be exercised at a thread count the test
# machine does not have to be restarted into: the alternative was a test that
# skipped itself whenever an earlier file had left the session narrow, which is
# a test that reports nothing and looks like a pass.
.ctBackendReportChunks <- function(cores, picked, threads = NULL) {
  cores <- suppressWarnings(as.integer(cores)[1L])
  picked <- suppressWarnings(as.integer(picked)[1L])
  if (is.na(cores) || is.na(picked) || picked < 1L) return(invisible(NULL))
  # Cheap test first. A fit that used what it asked for is the usual case and
  # must not pay a bridge call to establish that.
  if (picked * 2L > cores || cores - picked < 2L) return(invisible(NULL))
  # Guarded, because `juliaEval` *starts* a session when none is running rather
  # than failing, and a function whose whole job is to report must never be the
  # thing that launches the engine. In its own call site a fit has just
  # finished, so the session is up and this is one cheap bridge call.
  if (is.null(threads)) {
    threads <- if (!.ctJuliaSessionRunning()) NA_integer_ else
      tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
        error = function(e) NA_integer_)
  }
  threads <- suppressWarnings(as.integer(threads)[1L])
  limit <- if (is.na(threads)) cores else min(cores, threads)
  if (picked * 2L > limit || limit - picked < 2L) return(invisible(NULL))
  key <- paste0(limit, ":", picked)
  if (key %in% .ct_julia_cache$chunks_reported) return(invisible(NULL))
  .ct_julia_cache$chunks_reported <- c(.ct_julia_cache$chunks_reported, key)
  message("cores = ", cores, " requested, ", picked, " used.")
  invisible(NULL)
}

# The engine's per-iteration record, as a data frame.
#
# Flat vectors cross the bridge; the shape is rebuilt here. A missing or
# malformed trace is NULL rather than an error: nothing downstream requires it,
# and an old cached engine environment predating this will simply not send one.
#' @keywords internal
.ctBackendTrace <- function(trace) {
  if (is.null(trace)) return(NULL)
  columns <- lapply(trace, as.numeric)
  columns <- columns[vapply(columns, length, integer(1)) > 0L]
  if (!length(columns) || is.null(columns$iteration)) return(NULL)
  n <- length(columns$iteration)
  if (!all(vapply(columns, length, integer(1)) == n)) return(NULL)
  out <- as.data.frame(columns, stringsAsFactors = FALSE)
  out$iteration <- as.integer(out$iteration)
  out
}

# Run the engine's optimizer over a prepared specification.
#
# Factored out of .ctFitJuliaBackend() because cross-validation re-optimises the
# same model against held-out data (see .ctBackendLOO) and must do it exactly
# the way a fit does -- same tolerances, same gradient method, same thread cap.
# A second copy of this call would be a second set of defaults to keep in step.
# The engine takes the substep policy in one slot: `maxtimestep` as a number,
# or a mesh -- one substep count per row of the data -- as an integer vector.
# `.ctJuliaVector` keeps a one-row mesh a vector rather than the scalar
# JuliaConnectoR would make of it.
.ctJuliaSubstepArgument <- function(rule) {
  if (is.integer(rule) || length(rule) > 1L) .ctJuliaVector(as.integer(rule)) else as.numeric(rule)[1L]
}

# Choose the substep mesh for `spec` at the parameter values `values`, and
# return the spec carrying it together with a summary. The engine measures how
# nonlinear each interval was (substep_mesh.jl); a linear model comes back with
# the floor and no filter pass. A non-finite likelihood at `values` leaves the
# spec as it was, and the summary says so.
.ctJuliaAutoSubsteps <- function(spec, values) {
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(structure(spec, class = c("ctJuliaModel", "ctFitModel")))
  arguments <- list(objective, .ctJuliaNumericVector(values), tol = as.numeric(spec$substeps$tol),
    max_substeps = as.integer(spec$substeps$max), floor_rule = as.numeric(spec$substeps$floor))
  # Re-meshing at the optimum: the engine starts from one step per interval and
  # may find the filter not finite there on a stiff model; the mesh the fit was
  # found with is then its fallback.
  if (is.integer(spec$max_timestep)) arguments$fallback <- .ctJuliaVector(spec$max_timestep)
  chosen <- JuliaConnectoR::juliaGet(do.call(module$ctsem_auto_substeps, arguments))
  summary <- list(intervals = as.integer(chosen$intervals), refined = as.integer(chosen$refined),
    max_substeps = as.integer(chosen$max_substeps), total = as.integer(chosen$total),
    passes = as.integer(chosen$passes), finite = isTRUE(chosen$finite),
    tol = as.numeric(spec$substeps$tol), refit = FALSE)
  if (summary$finite) spec$max_timestep <- as.integer(chosen$mesh)
  list(spec = spec, summary = summary)
}

.ctJuliaSubstepMessage <- function(s) {
  if (!isTRUE(s$finite)) return("Substeps: likelihood not finite at the starting values; maxtimestep rule kept.")
  paste0("Substeps: ", s$refined, " of ", s$intervals, " intervals refined, max ", s$max_substeps,
    if (isTRUE(s$refit)) "; refit after the mesh moved at the optimum." else ".")
}

.ctJuliaOptimise <- function(model_spec, start, optimcontrol = list(),
  gradient = "adjoint", cores = 1L, verbose = 0L, maxiter = NULL,
  callback = NULL, objective = NULL, progress_label = NULL,
  progress_budget = FALSE) {
  spec <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  # A caller may hand in the objective to maximise. `intoverstates=FALSE` does,
  # passing the joint one over `[theta; z]`; everything below is unchanged by
  # that, because `ctsem_optimize` is typed on what the two have in common.
  if (is.null(objective)) objective <- .ctJuliaObjective(spec)
  module <- .ctJuliaModule(model_spec$project)
  # The shared optimcontrol vocabulary (see the table at the top of R/ctFit.R).
  # `tol` is the objective criterion on both backends -- mize's `abs_tol` there,
  # the engine's `f_tol` here -- and `g_tol` is the gradient norm on both. The
  # defaults are this backend's own and are the mirror image of stan's: stop on
  # the gradient, with the objective criterion off unless asked for.
  # `maxiter` as an argument is the caller overriding the cap for one stage (the
  # prior warm-up does), and beats the user's setting for that stage only.
  # Whether anything is reported at all, decided before the list below because
  # two of its elements need the same answer: `progress` runs the engine's
  # reporter and `progress_sink` says where its line goes. See the note on
  # `progress` at the foot of the list for why it is not simply `verbose`.
  reporting <- isTRUE(.ctJuliaOr(optimcontrol$progress,
    verbose > 0L || .ctProgressConsole()))
  common <- list(
    maxiter = as.integer(.ctJuliaOr(maxiter, .ctJuliaOr(optimcontrol$maxiter, 1000L))),
    g_tol = .ctJuliaOr(optimcontrol$g_tol, 1e-8),
    f_tol = .ctJuliaOr(optimcontrol$tol, 0),
    x_tol = .ctJuliaOr(optimcontrol$x_tol, 0),
    verbose = verbose > 0L,
    # Overwrite one line in place when someone is watching, and print
    # occasional separate lines when the output is going to a file or a knitr
    # chunk, where a carriage return is not a cursor movement. `verbose = 2`
    # keeps the history too, because at that point the point is the history.
    progress_overwrite = .ctProgressOverwrite(verbose),
    # Delivered as an R message rather than printed by the engine; see
    # `.ctProgressSink()`. Only when something is being reported at all -- a
    # silent fit should not pay a bridge callback per iteration.
    progress_sink = if (reporting) .ctProgressSink(
      .ctProgressOverwrite(verbose)) else NULL,
    # Names the stage on the progress line, defaulting to the engine's own
    # "optimise" when unset. `carefulfit` runs an optimisation before the fit's
    # and both carried that same label, so the counter ran up to the warm-up's
    # cap and then restarted from one -- which reads as a run that finished and
    # began again. The callback had the same fault, fixed separately; this is
    # its printed twin.
    progress_label = if (is.null(progress_label)) "optimise" else
      as.character(progress_label)[1L],
    # Whether `maxiter` is the plan or a safety limit, which decides whether
    # the counter carries a denominator at all. The warm-up's cap of 10 is
    # always reached, so "7/10" is honest there; the fit's cap of 1000 is not a
    # target and showing "297/1000" invited reading a completion fraction that
    # meant nothing. See `_progress_optimise` in progress.jl.
    progress_budget = isTRUE(progress_budget),
    # `verbose = 2` asks for the history, and the engine's non-overwriting
    # cadence is 5s against 0.4s for the overwriting one -- so raising the
    # level printed *fewer* progress lines than leaving it alone, and on any
    # fit shorter than five seconds it printed none at all. Asking for more
    # detail must not give less of it.
    # Keyed on the level, not on `!.ctProgressOverwrite()`: that is also FALSE
    # for a `verbose = 1` run whose output goes to a log file, where the slow
    # cadence is right and hundreds of lines are not.
    progress_every = if (is.numeric(verbose) && length(verbose) == 1L &&
      !is.na(verbose) && verbose >= 2) 0.4 else 0,
    # On when someone is watching, which is not the same question as how
    # verbose to be. A default fit used to print two lines and then nothing at
    # all however long it ran, because progress was tied to `verbose` and
    # `verbose` defaults to 0 -- so the reporting existed and almost nobody
    # saw it. Keyed on the same console detection the overwriting uses, so a
    # script or a knitr chunk still gets nothing, and overridable with
    # `optimcontrol$progress`.
    progress = reporting)
  # A live callback into R, for a front end that wants to draw the trace as it
  # happens rather than read it afterwards. The engine calls it on the same time
  # cadence as the printed line, not once per iteration: measured through
  # JuliaConnectoR a callback costs 0.5 ms, which is nothing occasionally and
  # 17% of a three-second fit at every iteration.
  failure <- NULL
  if (!is.null(callback)) {
    if (!is.function(callback)) {
      stop("optimcontrol$callback must be a function of (iteration, total, ",
        "objective, gradient_norm).", call. = FALSE)
    }
    # Caught here rather than in the engine, because an error thrown out of an
    # R callback does not reach the engine at all: it aborts before a reply is
    # sent, JuliaConnectoR reports "Message type not supported (yet)", and the
    # bridge is left desynchronised -- the fit is lost to a fault in a
    # *reporting* function. A callback is a convenience and must never be able
    # to do that, so it is wrapped to swallow everything and disable itself.
    #
    # The message is stored rather than warned immediately: `options(warn = 2)`
    # would turn the warning into exactly the error this exists to prevent.
    alive <- TRUE
    common$progress_callback <- function(iteration, total, objective,
      gradient_norm) {
      if (alive) {
        tryCatch(callback(iteration, total, objective, gradient_norm),
          error = function(e) {
            alive <<- FALSE
            failure <<- conditionMessage(e)
          })
      }
      NULL
    }
  }
  # Exposed because it is the one optimiser knob that measurably changed both
  # speed and whether the gradient criterion was met; the engine's default is
  # documented at `_CTSEM_LBFGS_MEMORY`.
  if (!is.null(optimcontrol$lbfgs_memory)) {
    common$lbfgs_memory <- as.integer(optimcontrol$lbfgs_memory)[1L]
  }
  # `cores` is the ceiling; `ctsem_tune_chunks!` measures the count to use
  # within it, and the fit records what it picked. Restored afterwards so the
  # session does not carry this fit's ceiling into the next thing that runs.
  result <- .ctBackendWithMaxChunks(cores, {
    if (!is.null(model_spec$laplace)) {
      # `gradient` selects how the *process* likelihood's gradient is taken and
      # does not apply here: the Laplace objective's gradient is a forward sweep
      # over that reverse pass regardless.
      JuliaConnectoR::juliaGet(do.call(module$ctsem_laplace_optimize,
        c(list(objective, .ctJuliaNumericVector(start)), common)))
    } else {
      JuliaConnectoR::juliaGet(do.call(module$ctsem_optimize,
        c(list(objective, .ctJuliaNumericVector(start)), common,
          list(gradient_method = gradient))))
    }
  })
  # What the tuner settled on, when that is well short of what was asked for.
  .ctBackendReportChunks(cores, result$chunks)
  if (!is.null(failure)) {
    warning("The progress callback failed and was disabled after the first ",
      "error; the fit itself is unaffected. The error was: ", failure,
      call. = FALSE)
  }
  result
}

.ctFitJuliaBackend <- function(datalong, model, prepared_data = NULL, inits = NULL, cores = 1L,
  optimcontrol = list(), verbose = 0L, fit = TRUE,
  priors = FALSE, intoverpop = "augmented", optimize = TRUE, chains = 4L,
  iter = 2000L, control = list(), intoverstates = TRUE) {
  .ctJuliaInterruptSafe(.ctFitJuliaBackendImpl(datalong = datalong,
    model = model, prepared_data = prepared_data, inits = inits, cores = cores,
    optimcontrol = optimcontrol,
    verbose = verbose, fit = fit, priors = priors, intoverpop = intoverpop,
    optimize = optimize, chains = chains, iter = iter, control = control,
    intoverstates = intoverstates))
}

#' @keywords internal
.ctFitJuliaBackendImpl <- function(datalong, model, prepared_data = NULL,
  inits = NULL, cores = 1L,
  optimcontrol = list(), verbose = 0L, fit = TRUE,
  priors = FALSE, intoverpop = "augmented", optimize = TRUE, chains = 4L,
  iter = 2000L, control = list(), intoverstates = TRUE) {
  # The engine environment and a session restart are ctJuliaSetup()'s own
  # `project` and `force = TRUE`; the per-fit copies of both were a second
  # spelling of an exported function, so a fit uses whatever session is set up.
  # `cores` splits the engine's subject loop. It is requested as a Julia thread
  # count before the session starts (which is the only time that can be set),
  # and capped per fit afterwards, so a session started with more threads is not
  # forced to use them all. A session started with fewer cannot be widened, and
  # `.ctBackendResolveThreads()` says so rather than running serially in silence.
  requested <- cores
  cores <- suppressWarnings(as.integer(cores)[1L])
  if (is.na(cores)) {
    # Silently falling back to one core is how `cores='maxneeded'` went
    # unnoticed: a fit that should have used the machine used a single thread
    # and said nothing. ctFit() resolves that keyword now, so anything still
    # arriving unparseable here is a caller mistake worth hearing about.
    warning("cores=", deparse(requested)[1L], " is not a number of cores; ",
      "using 1. Pass an integer.", call. = FALSE)
    cores <- 1L
  }
  cores <- max(1L, cores)
  # Requested as a Julia thread count while that is still possible, and said
  # out loud when it is not: a session already running at fewer threads than
  # `cores` caps the subject loop at its own count, and used to do so silently.
  # Before `.ctJuliaPrepare()` deliberately -- the opt-in restart ends the
  # process, and nothing prepared here may be alive across that.
  .ctBackendResolveThreads(cores, report = isTRUE(fit))
  gradient <- .ctJuliaOr(optimcontrol$gradient, "adjoint")
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
    priors = priors, intoverpop = intoverpop, optimize = optimize,
    tipredMissingIncludeOutcome = .ctJuliaOr(optimcontrol$tipredMissingIncludeOutcome, TRUE))
  if (!is.null(model_spec$ti_missing) && nrow(model_spec$ti_missing)) {
    # The state-explicit route (`intoverstates=FALSE`) samples the latent
    # trajectory through a different objective (`CTSEMJointObjective`,
    # state_sampling.jl) that this feature has not touched at all -- not
    # "not verified", genuinely absent. Refusing here is the same choice as
    # the `intoverpop='laplace'/'none'` refusal inside `.ctJuliaPrepare()`,
    # for the same reason: no silent gap in the gradient.
    if (isFALSE(intoverstates)) {
      stop("Julia backend cannot sample missing TI predictor values (",
        nrow(model_spec$ti_missing), " missing) with intoverstates=FALSE; ",
        "the state-explicit route does not yet support this. Impute them ",
        "before fitting, drop those subjects, use intoverstates=TRUE, or ",
        "use backend='stan'.", call. = FALSE)
    }
  }
  if (!fit) return(structure(model_spec, class = c("ctJuliaModel", "ctFitModel")))
  # A sampled TI-predictor value used to have no adjoint cotangent, so a fit
  # here was refused unless `gradient='forward'` was requested by name (see
  # the git history of this function, and SPEC-tipred-sampling.md). The
  # reverse pass now covers it -- `_ctsem_ti_pullback!`'s `TIMissingRecipe`
  # method adds the TI-effect product-rule term and `ctsem_adjoint_gradient`
  # adds `_ctsem_ti_missing_loglik_gradient!`'s term for the imputation
  # log-density itself, both cross-checked against ForwardDiff and FiniteDiff
  # in test_ti_missing_predictor.jl -- so `gradient` is left exactly as
  # resolved above: 'adjoint' by default, 'forward' still available on
  # request, neither one singled out for this kind of model any more.

  # `optimize=FALSE` fits by sampling. Which sampler is decided by
  # `intoverpop`, which says what has already been integrated out; see
  # `.ctJuliaSampleFit`.
  if (!isTRUE(optimize)) {
    return(.ctJuliaSampleFit(model_spec, datalong = datalong, model = model,
      inits = inits, cores = cores,
      optimcontrol = optimcontrol, chains = chains, iter = iter,
      control = control, priors = priors, intoverpop = intoverpop,
      gradient = gradient, verbose = verbose,
      intoverstates = intoverstates))
  }

  npar <- .ctBackendNpar(model_spec)
  # A fully fixed model has nothing to maximise over. Without the zero above,
  # `max` warned and returned -Inf, and `rnorm(-Inf, ...)` then failed with
  # "invalid arguments", which says nothing about the model. Refused rather
  # than run: a zero-length start would put the engine's L-BFGS on a
  # zero-dimensional problem, which nothing here tests. Evaluating a fixed
  # model's likelihood is a fair thing to want, and `fit = FALSE` still
  # returns the prepared model to evaluate.
  if (npar < 1L) {
    stop("This model has no free parameters, so there is nothing to optimise. ",
      "Free a parameter, or use fit = FALSE to prepare the model without ",
      "fitting it.", call. = FALSE)
  }
  start <- .ctJuliaInitialValues(npar, inits,
    initsd = .ctJuliaOr(optimcontrol$initsd, .01))
  # Starting values read off the data, for the diagonals whose defaults are
  # guesses about the data's scale. See R/ctDataStart.R for what is derived and
  # why; `optimcontrol$datastart = FALSE` restores the fixed start. Supplied
  # `inits` are never overridden -- a starting value the caller chose is the
  # one thing here that is not a guess.
  datastart <- .ctJuliaOr(optimcontrol$datastart, TRUE)
  if (isTRUE(datastart) && is.null(inits) && !is.null(datalong)) {
    derived <- try(.ctDataStart(datalong, model, model_spec, npar), silent = TRUE)
    if (!inherits(derived, "try-error") && !is.null(derived)) {
      use <- is.finite(derived)
      # Exactly, not jittered: where the data decided the value, two runs of
      # the same fit should start in the same place.
      start[use] <- derived[use]
      if (verbose > 0) message("Starting values derived from the data for ",
        sum(use), " of ", npar, " parameters.")
    }
  }
  # The state-explicit target, when asked for: the same optimiser over a longer
  # vector. The innovations start at zero, which is both their prior mode and
  # the trajectory the parameters alone imply -- there is nothing better to
  # start them at and nothing arbitrary about it.
  jointobjective <- NULL
  nstate <- 0L
  if (!isTRUE(intoverstates)) {
    jointobjective <- .ctJuliaJointObjective(model_spec, npar)
    nstate <- .ctJuliaStateDimension(model_spec)
    start <- c(start, numeric(nstate))
  }
  # The prior-warmed spec, or NULL when priors cannot be mapped onto this
  # model's raw layout. `.ctBackendLaplacePriorSpec` refuses shapes it cannot
  # map rather than mis-assigning priors across levels, so this is allowed to
  # fail and simply leave the fit unwarmed.
  warmspec <- function() {
    if (isTRUE(priors) || is.null(prepared_data)) return(NULL)
    spec <- model_spec
    spec$priors <- try(if (!is.null(model_spec$laplace))
      .ctBackendLaplacePriorSpec(prepared_data, model_spec$laplace, npar)
      else .ctBackendPriorSpec(prepared_data, npar), silent = TRUE)
    if (inherits(spec$priors, "try-error") || !length(spec$priors$index))
      return(NULL)
    spec
  }
  # `carefulfit`, the same argument `stanoptimis` takes and with the same
  # meaning: a rough first pass with the priors on, to get starting values,
  # when priors are otherwise off. Stan's version caps that pass at 50
  # iterations and loosens its tolerance; this one caps it at 10, which the
  # measurements below settled, and lets the cap do the work rather than also
  # loosening the tolerance.
  #
  # It does not make fits faster. Measured over 800 fits across both engine
  # routes and four measurement types, total iterations against a plain fit
  # came to 0.91-1.14 at ten prior iterations, 1.09-1.50 at twenty and
  # 1.47-1.64 at forty. The second stage does converge in fewer iterations --
  # 40 to 26 on ordinal Laplace at a cap of twenty -- but not by enough to pay
  # for the first stage.
  #
  # What it buys is where the fit lands. Over those 800 fits the warmed start
  # was never worse in any cell, and three cells were better: an ordinal
  # Laplace fit that failed at -2498.24 and warmed reaches -2207.06; two binary
  # fits by about two log units; and a mixed-indicator model that *converged*
  # at -1922.81 unwarmed and at -1833.42 warmed. That last one is the case this
  # is for -- it reports success while returning a random-effect SD of 7.23
  # against a truth of 0.5, and it was the single replication inflating that
  # condition's RMSE to 0.88 in the recorded study. Warmed, the same data give
  # 0.544.
  #
  # A longer pass is not a safer one, and the cap is 10 rather than 20 or 50
  # because of it. Over 720 fits of three measurement types on both routes,
  # judged on the random-effect SD against a truth of 0.5 rather than on the
  # log likelihood -- a higher likelihood in the wrong basin is still the wrong
  # answer:
  #
  #   cap 10   240/240 converged, RMSE 0.190, worst error 0.547
  #   cap 20   240/240 converged, RMSE 0.550, worst error 7.991
  #   off      236/240 converged, RMSE 0.478, worst error 6.732
  #
  # Both caps fix convergence and the two are identical to four decimals on
  # every condition but one. All of the difference is the mixed-indicator model
  # under state augmentation, where a cap of 20 is *worse than not warming up
  # at all* -- it leaves the local optimum with the SD of 8.0 in place, where a
  # cap of 10 removes it outright. Median wall time is the same either way: 20
  # buys back its extra prior iterations in the second stage and no more.
  #
  # The prior pass pulls the start toward the prior mode, and past about ten
  # iterations that is what it hands the likelihood.
  careful <- optimcontrol$carefulfit
  if (is.null(careful)) careful <- TRUE
  warmiter <- if (isTRUE(careful)) 10L else
    if (is.numeric(careful) && length(careful) == 1L && careful >= 1)
      as.integer(careful) else 0L
  # `stanoptimis` turns `carefulfit` off when starting values were supplied,
  # since the point of the pass is to produce some. Overriding a starting value
  # the caller chose would be worse than surprising.
  if (!is.null(inits) && !identical(inits, "random")) warmiter <- 0L
  # And not at all on the state-explicit target. The warm-up exists to place
  # the *population* parameters from the priors; the innovations already
  # start at their own prior mode, and running a second optimisation over
  # the whole extended vector to rediscover that would cost as much as the
  # fit it is warming.
  if (!isTRUE(intoverstates)) warmiter <- 0L
  if (warmiter >= 1) {
    spec <- warmspec()
    if (!is.null(spec)) {
      # No callback here, deliberately. This stage is a starting-value device,
      # not the fit: it optimises a *different* objective (the posterior rather
      # than the likelihood) and its result is used only as `start` below.
      #
      # Passing the user's callback through gave it two runs of iteration
      # numbers -- 1..10 from here, then 1..N from the fit -- so a front end
      # drawing a live trace saw the counter restart, and
      # `max(seen$iterations)` reported this stage's cap rather than the fit's
      # iteration count. That broke the contract the callback shares with
      # `fit$trace` and `fit$estimate$iterations`, both of which describe the
      # fit alone. Reported as test-julia-trace.R expecting 8 and seeing 10,
      # which is exactly `warmiter`.
      warmed <- try(.ctJuliaOptimise(spec, start, optimcontrol = optimcontrol,
        maxiter = as.integer(warmiter),
        gradient = gradient, cores = cores, verbose = verbose,
        callback = NULL, progress_label = "prior warm-up",
        # This stage runs its cap and stops; the cap is the plan, not a limit
        # it is trying to stay under, so its counter keeps a denominator.
        progress_budget = TRUE), silent = TRUE)
      # A warm start is only a starting value: if it produced numbers the fit
      # can use, use them, and otherwise start where we would have anyway.
      if (!inherits(warmed, "try-error")) {
        values <- as.numeric(warmed$minimizer)
        if (length(values) == npar && all(is.finite(values))) start <- values
      }
    }
  }
  # Automatic substeps: choose the mesh at the starting values, fit, choose
  # it again at the optimum, and refit once from there if it moved. Frozen
  # within each fit so the objective stays smooth; see substep_mesh.jl. The
  # state-explicit target has one innovation per substep, so its objective
  # and start vector are rebuilt whenever the mesh changes.
  substeps <- NULL
  remesh <- function(values) {
    chosen <- .ctJuliaAutoSubsteps(model_spec, values)
    model_spec <<- chosen$spec
    if (!is.null(jointobjective)) {
      jointobjective <<- .ctJuliaJointObjective(model_spec, npar)
      nstate <<- .ctJuliaStateDimension(model_spec)
    }
    chosen$summary
  }
  if (!is.null(model_spec$substeps)) {
    substeps <- remesh(start[seq_len(npar)])
    if (!is.null(jointobjective)) start <- c(start[seq_len(npar)], numeric(nstate))
  }
  result <- .ctJuliaOptimise(model_spec, start, optimcontrol = optimcontrol,
    gradient = gradient, cores = cores, verbose = verbose,
    callback = optimcontrol$callback, objective = jointobjective)
  if (!is.null(substeps) && isTRUE(substeps$finite)) {
    before <- model_spec$max_timestep
    optimum <- as.numeric(result$minimizer)[seq_len(npar)]
    substeps <- remesh(optimum)
    if (!identical(model_spec$max_timestep, before)) {
      substeps$refit <- TRUE
      restart <- if (is.null(jointobjective)) optimum else c(optimum, numeric(nstate))
      result <- .ctJuliaOptimise(model_spec, restart, optimcontrol = optimcontrol,
        gradient = gradient, cores = cores, verbose = verbose,
        callback = optimcontrol$callback, objective = jointobjective)
    }
  }
  if (!is.null(substeps)) message(.ctJuliaSubstepMessage(substeps))
  # There is deliberately no second, after-the-fact prior restart here.
  #
  # An earlier version retried a non-converged fit from a full prior
  # optimisation. `carefulfit` above makes that redundant: it warms every fit
  # from the priors already, and over 720 fits at a cap of ten it converged
  # 240 out of 240 in each condition, so the retry had nothing left to rescue.
  # What it did instead was move answers. A fit deliberately capped at one
  # iteration from supplied starting values came back from raw 24 at raw 12.8 --
  # not the fit that was asked for, and reported without comment. Two
  # mechanisms for one job, where the second can only act in cases the first
  # did not fix, is a way to be surprised rather than a safety net.
  # The engine maximises the log posterior, so its `maximum_loglik` is the log
  # posterior and the per-subject objectives (which carry no prior term) sum to
  # the log likelihood. Without priors the two are the same number; with them
  # they are not, and the summary reports both, as it does for Stan.
  subject_loglik <- as.numeric(result$subject_loglik)
  loglik <- if (length(subject_loglik)) sum(subject_loglik) else
    as.numeric(result$maximum_loglik)
  # `[theta; z]` comes back as one vector and is split here, so that
  # `estimate$raw` means the same thing on every fit: the population
  # parameters, and nothing else. Everything downstream -- the summary, the
  # transforms, prediction -- reads that and needs no notion of a state
  # block. The trajectory is kept beside it rather than folded in.
  minimizer <- as.numeric(result$minimizer)
  gradientvec <- as.numeric(result$gradient)
  states <- NULL
  innovations <- NULL
  if (!isTRUE(intoverstates)) {
    innovations <- minimizer[npar + seq_len(nstate)]
    states <- .ctJuliaJointStates(model_spec, jointobjective, minimizer)
    minimizer <- minimizer[seq_len(npar)]
    gradientvec <- gradientvec[seq_len(npar)]
  }
  # `$data` is not set here: it is attached by ctFit.R after this call
  # returns, as the sentinel-cleaned `standata` structure (matching a stan
  # fit's `$data`/`$standata` split -- see R/ctFit.R). Setting the verbatim
  # `datalong` here would only be overwritten there in the normal ctFit()
  # entry point, and `model_spec$data` (the long data.frame `.ctJuliaPrepare()`
  # kept) already covers what code inside this package needs from it, via
  # `.ctBackendSpec()`.
  out <- list(backend = "julia", model = model, model_spec = model_spec,
    estimate = list(raw = minimizer,
      loglik = loglik,
      logposterior = as.numeric(result$maximum_loglik),
      gradient = gradientvec,
      subject_loglik = result$subject_loglik, converged = isTRUE(result$converged),
      iterations = as.integer(result$iterations),
      # What `ctsem_tune_chunks!` measured as the best subject-chunk count
      # within the `cores` ceiling. Carried onto the fit because the uncertainty
      # phase reads it rather than re-deriving it from `cores`: the subject loop
      # is not monotone in the chunk count, which is the whole reason the tuner
      # exists.
      chunks = if (is.null(result$chunks)) NA_integer_ else as.integer(result$chunks),
      # Evaluation counts, because "how many times did it call the likelihood"
      # is the first question about a fit that took longer than expected, and
      # it was previously only obtainable by timing one evaluation and dividing.
      f_calls = if (is.null(result$f_calls)) NA_integer_ else as.integer(result$f_calls),
      g_calls = if (is.null(result$g_calls)) NA_integer_ else as.integer(result$g_calls),
      # The gradient at the estimate, and the tolerance it was judged against.
      # `converged` is one bit and a fit that stops just short of a tolerance
      # looks the same as one that never moved; these are what tell them apart.
      gradient_norm = if (is.null(result$gradient_norm)) NA_real_ else
        as.numeric(result$gradient_norm),
      gradient_tolerance = if (is.null(result$scaled_tolerance)) NA_real_ else
        as.numeric(result$scaled_tolerance),
      # A parameter that reached the flat region of its transform, and which.
      # Invisible in the gradient -- a saturated transform reports a gradient
      # of zero, which passes every tolerance -- so without this the warning
      # below described such a fit by its gradient alone and read as though it
      # had passed.
      #
      # `saturated` is a statement about identification, not about
      # convergence: the usual cause is a population standard deviation with
      # no individual differences behind it, which is a finding. `overshot` is
      # the half of it that *is* a convergence failure -- the optimizer
      # overstepped and the point is not a maximum -- and it is what
      # `converged` is keyed on. See `_ctsem_overshot` in the engine.
      saturated = isTRUE(result$saturated),
      # 0 from the engine means "none"; never an empty vector, which deadlocks
      # the bridge. Named where names are available, because "raw parameter 10"
      # is not something a user can act on and "rawcor_mm2__mm1" is.
      saturated_parameters = .ctJuliaSaturatedNames(result, model_spec, npar),
      overshot = isTRUE(result$overshot),
      overshoot_gain = if (is.null(result$overshoot_gain)) NA_real_ else
        as.numeric(result$overshoot_gain),
      # Whether the first pass with priors ran, and how long it was allowed.
      carefulfit = warmiter >= 1, carefulfit_iterations = as.integer(warmiter),
      # Which line search produced the answer. "hagerzhang+backtracking" means
      # Hager-Zhang stopped short and the fit was finished by the fallback.
      linesearch = if (is.null(result$linesearch)) NA_character_ else
        as.character(result$linesearch),
      stalled = isTRUE(result$stalled),
      # What nsubsteps = 'auto' decided: intervals, how many were refined, the
      # largest count, the total, and whether the fit was redone after the
      # mesh moved at the optimum. NULL unless it was asked for.
      substeps = substeps,
      # State-explicit fits only. `states` is the trajectory at the
      # estimate, rows by latents, and `innovations` the standard normal
      # vector it was built from -- the same thing `ctGenerate` is handed,
      # so a fitted trajectory can be replayed.
      #
      # `loglik_type` says what `loglik` is. On the marginal route it is a
      # marginal log likelihood; here it is the *joint* density of the data
      # and the states, which is a different quantity and is not comparable
      # with one -- an information criterion computed across the two would
      # be meaningless, and this is what says so.
      states = states, innovations = innovations,
      loglik_type = if (isTRUE(intoverstates)) "marginal" else "joint"),
    engine = model_spec$engine,
    # Every iteration, recorded whatever `verbose` said. It costs a push onto a
    # vector in Julia and one transfer at the end, and the fit whose trace turns
    # out to be worth looking at is exactly the one nobody thought to turn
    # reporting on for.
    trace = .ctBackendTrace(result$trace),
    args = list(backend = "julia",
      optimcontrol = optimcontrol, cores = cores, priors = priors,
      intoverpop = intoverpop, intoverstates = isTRUE(intoverstates)))
  # An optimizer that ends where it started has not fitted anything, whatever
  # its convergence flags say -- and Optim's own verdict is a disjunction that a
  # failed first line search satisfies trivially. Saying so here is what stops a
  # simulation study from averaging over starting values it never left.
  if (isTRUE(result$stalled)) {
    warning("The optimizer made no progress from its starting values, and the ",
      "gradient there is not zero. Treat this fit as failed: check the starting ",
      "values, and see fit$estimate$stalled.", call. = FALSE)
  } else if (isTRUE(result$overshot)) {
    # Reported separately because the gradient says nothing useful here: the
    # coordinate's transform is flat, so the gradient underflows to zero and
    # the point passes every tolerance while not being a maximum at all. The
    # engine establishes that by pulling the coordinate back and finding the
    # objective improves -- `overshoot_gain` is by how much.
    warning("The optimizer overstepped into a region where the transform of ",
      paste(out$estimate$saturated_parameters, collapse = ", "),
      " is flat to machine precision, and the objective improves by ",
      signif(as.numeric(result$overshoot_gain), 3), " when that parameter is ",
      "pulled back, so this is not a maximum. Treat this fit as failed and ",
      "check the starting values. See fit$estimate$overshot.", call. = FALSE)
  } else if (!isTRUE(result$converged)) {
    warning("The optimizer stopped without meeting its convergence criterion: ",
      "largest gradient ", signif(as.numeric(result$gradient_norm), 3),
      " against a tolerance of ",
      signif(as.numeric(result$scaled_tolerance), 3),
      ". The estimate may still be usable -- compare the two, and see ",
      "fit$estimate$gradient_norm.", call. = FALSE)
  }
  if (!is.null(model_spec$laplace)) {
    # The inner solve is part of the objective, so its status is part of
    # whether the fit means anything. Kept on the fit rather than printed and
    # discarded, since a subject whose mode did not converge contributes a term
    # that is not the integral it is supposed to approximate.
    out$laplace <- list(
      nrandom = model_spec$laplace$nrandom,
      param = model_spec$laplace$param,
      nlevels = model_spec$laplace$nlevels,
      levels = lapply(model_spec$laplace$levels, function(x)
        list(name = x$name, param = x$param, nrandom = x$nrandom,
          ngroups = x$ngroups)),
      linesearch = if (is.null(result$linesearch)) NA_character_ else
        as.character(result$linesearch),
      inner_converged = isTRUE(result$inner_converged),
      inner_iterations = as.integer(result$inner_iterations),
      # Two different things. `hessian_repaired` is true if *any* Newton iterate
      # for that unit needed its curvature shifted, which is ordinary behaviour
      # for a nonlinear model on the way to a mode. `mode_repaired` is true if
      # the curvature at the reported mode needed it, which is the one that
      # says the approximation there is questionable.
      hessian_repaired = as.logical(result$hessian_repaired),
      mode_repaired = if (is.null(result$mode_repaired)) NA
        else as.logical(result$mode_repaired))
    if (!isTRUE(result$inner_converged)) {
      warning("The random-effect mode did not converge for every subject; ",
        "see fit$laplace$inner_converged.", call. = FALSE)
    }
    # A correlation sitting on its cap means the data could not locate it --
    # usually a level with too few groups for the number of effects asked of
    # it. The fit is still a fit; the estimate for that pair is not.
    boundary <- .ctJuliaBoundaryReport(out, model_spec, result)
    out$laplace$boundary <- boundary
    if (nrow(boundary)) {
      warning("Random-effect correlation", if (nrow(boundary) > 1) "s" else "",
        " reached the boundary and ", if (nrow(boundary) > 1) "were" else "was",
        " capped at |r| = ", signif(boundary$correlation[1L], 3), ": ",
        paste(sprintf("%s (%s level)", boundary$param, boundary$level),
          collapse = ", "),
        ". This usually means the level has too few groups to identify that ",
        "correlation. See fit$laplace$boundary.", call. = FALSE)
    }
  }
  class(out) <- c("ctJuliaFit", "ctFit")

  # Uncertainty is part of fitting, not a separate step the user has to know to
  # take -- `stanoptimis()` finishes every optimized Stan fit the same way, and
  # a `summary()` that silently reported point estimates only, purely because of
  # the backend chosen, is the difference this exists to remove. The control
  # names are `stanoptimis()`'s, so `optimcontrol` means the same thing to both
  # backends; `optimcontrol$estonly` skips it, as it does there.
  # Not for an optimised state-explicit fit. The only curvature available
  # there is the profile's, and with about as many innovations as
  # observations the profile is nearly flat: the states re-optimise to
  # absorb almost any change in the parameters, which is what profiling
  # does and is why it is not an observed information. The term that
  # identifies the parameters -- the log determinant the Laplace marginal
  # carries -- is exactly the one profiling drops. Measured on fifteen
  # subjects of six Gaussian rows with four free parameters, the largest
  # eigenvalue of the profiled curvature was 0.05.
  #
  # Reporting intervals from that would be numbers that look like standard
  # errors and are not, which is worse than saying so. The matrix is kept
  # on the fit for anyone who wants to look, and an explicit
  # ctOptimUncertainty(fit) still computes from it.
  if (!isTRUE(intoverstates)) {
    profiled <- try(.ctBackendJointHessian(out, out$estimate$raw),
      silent = TRUE)
    if (!inherits(profiled, "try-error")) out$estimate$hessian_profile <- profiled
    message("No standard errors: an optimised intoverstates=FALSE fit has ",
      "only the profile curvature, which the states flatten. Use ",
      "optimize=FALSE to sample them, or intoverstates=TRUE.")
  } else if (!isTRUE(optimcontrol$estonly)) {
    uncertainty <- .ctJuliaOr(optimcontrol$uncertainty, "hessian")
    out <- ctOptimUncertainty(fit = out, uncertainty = uncertainty,
      draws = .ctJuliaOr(optimcontrol$uncertaintyDraws, "auto"),
      finishsamples = .ctJuliaOr(optimcontrol$finishsamples, 1000L),
      cores = cores, control = .ctJuliaOr(optimcontrol$uncertaintyControl, list()),
      verbose = verbose)
  }

  # The draws pushed through the model's transforms, once, exactly as the Stan
  # path stores `stanfit$transformedpars` at fit time. Every summary, extract
  # and system-matrix collapse reads this rather than asking the engine again.
  #
  # `.ctBackendConstrained`, not `.ctBackendConstrain`: `ctOptimUncertainty()`
  # has already done exactly this work for exactly these draws a few lines
  # above, and calling the uncached form here repeated it in full. That was
  # 4.5 s of a 16.7 s fit -- 27% of it -- computed twice and thrown away once,
  # and it is the single largest reason a julia fit was slower end to end than
  # the Stan path whose numerics it beats. The cached form recomputes only when
  # the draws differ, which is what happens when uncertainty was skipped and
  # `out$transformedpars` is still NULL.
  out$transformedpars <- .ctBackendConstrained(out)

  # The filter output at the estimate, cached as the Stan path caches
  # `stanfit$kalman`: summary()'s standardised residual covariance reads it, and
  # recomputing it per summary call would repeat a whole filter pass.
  #
  # The Laplace route filters each subject at its own estimated random effects,
  # which is the smoothed equivalent of what the augmented route's carrier
  # states give -- see `.ctBackendKalmanRaw`, which says so at the point of use.
  # Only the prior prediction errors, not the whole filter pass. Every summary
  # that reads a cached filter reads exactly `errprior`; the other seventeen
  # arrays were transferred and stored so that one could be. See
  # `.ctBackendPriorErrors`.
  out$priorerrors <- .ctBackendPriorErrors(out)

  # What the fit can honestly say about its own credibility. Not a verdict on
  # whether the answer is sensible -- that is a judgement about the model and
  # the data -- but a statement of which directions the data does not determine,
  # and therefore which reported intervals do not mean what they appear to.
  rawnames <- .ctBackendRawParameterNames(out, length(out$estimate$raw))
  # `fit`/`at` let it tell a random-effect block trading its scale off against
  # its correlations -- where the covariances are determined, and fixing a
  # value throws them away -- from a direction the data says nothing about.
  out$identifiability <- .ctBackendIdentifiability(out$uncertainty$hessian,
    rawnames, fit = out, at = out$estimate$raw)
  # `$uncertainty$intervalcheck` is attached by `.ctBackendUncertainty()`, so
  # it describes whichever method ran; it is only warned about here. A separate
  # question from identifiability: a direction can be flat enough to ruin every
  # interval that touches it without being flat enough to be called
  # unidentified, and that case reaches the user as a plausible fit with
  # suspiciously wide intervals.
  out$collapsedScales <- .ctBackendCollapsedScales(out)
  .ctBackendIdentifyWarn(out$identifiability, out$collapsedScales,
    out$uncertainty$intervalcheck)
  out
}

#' @export
print.ctJuliaFit <- function(x, ...) {
  cat("ctsem Julia fit\n")
  cat("  log likelihood:", format(x$estimate$loglik), "\n")
  cat("  converged:", x$estimate$converged, " iterations:", x$estimate$iterations, "\n")
  # A sampled fit's headline is not the optimiser's. That optimisation ran only
  # to place the sampler and build its metric, so a fit whose chains never
  # agreed prints `converged: TRUE` on the line above and is still worthless.
  if (!is.null(x$sample)) {
    cat("  sampled: ", x$sample$chains, " chains x ", x$sample$draws,
      " draws; chains converged: ",
      if (is.null(x$sample$converged) || is.na(x$sample$converged)) "unknown" else
        as.character(isTRUE(x$sample$converged)), "\n", sep = "")
    if (length(x$sample$diagnosis)) {
      cat("  ", paste(utils::head(x$sample$diagnosis, 3), collapse = "; "),
        ". See fit$sample.\n", sep = "")
    }
  }
  # One line, only when there is something to say. A reported interval much
  # wider than the curvature at the estimate supports is not visible anywhere
  # in the numbers themselves -- it was found once only by fitting the same
  # data twice -- so it is said here, where a single fit is looked at.
  ivc <- x$uncertainty$intervalcheck
  if (!is.null(ivc) && isTRUE(ivc$nflagged > 0L)) {
    cat("  ", ivc$nflagged, " interval(s) far wider than the curvature supports: ",
      paste(utils::head(ivc$parameters, 4), collapse = ", "),
      ". See fit$uncertainty$intervalcheck.\n", sep = "")
  }
  # And the opposite, which needs saying here more than the above does: a
  # coordinate lying in a direction with no curvature keeps none of the
  # variance the projected inverse has left to give, so it prints as tightly
  # estimated. Nothing in the numbers themselves shows it.
  if (!is.null(ivc) && isTRUE(ivc$nunidentified > 0L)) {
    cat("  ", ivc$nunidentified,
      " parameter(s) with no curvature, so no interval: ",
      paste(utils::head(ivc$unidentified, 4), collapse = ", "),
      ". summary() reports these as NA.\n", sep = "")
  }
  invisible(x)
}

#' Plots for ctJuliaFit objects
#'
#' The \code{backend='julia'} counterpart of \code{\link{plot.ctStanFit}}, and a
#' subset of it: the two plot types that describe a fitted model rather than a
#' sampler. \code{'regression'} plots model implied regression coefficients
#' against the time interval via \code{\link{ctDiscretePars}}, and
#' \code{'kalman'} plots expectations via \code{\link{ctPredict}}. The remaining
#' types \code{plot.ctStanFit} offers -- prior/posterior densities, traces,
#' intervals -- all read Stan's own sample object, which an optimized julia fit
#' has no analogue of.
#'
#' @param x Fit object from \code{ctFit(..., backend='julia')}.
#' @param types Vector of plot types: 'all', 'regression', 'kalman'.
#' @param wait Logical. Pause between plots?
#' @param ... Passed through to \code{\link{ctDiscretePars}} and
#'   \code{\link{ctPredict}}. Beware clashes when \code{types='all'}.
#' @return Nothing. Generates plots.
#' @method plot ctJuliaFit
#' @examples
#' \donttest{
#' # Needs a working julia backend, so the example is inert where Julia is
#' # absent -- including on CRAN, whose check machines have none.
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   data <- ctstantestdat
#'   model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
#'     LAMBDA = matrix(1))
#'   fit <- ctFit(data, model, backend = 'julia')
#'   plot(fit, wait = FALSE)
#' }
#' }
#' @export
plot.ctJuliaFit <- function(x, types = "all", wait = TRUE, ...) {
  available <- c("regression", "kalman")
  if (identical(types[1L], "all")) {
    types <- available
    if (!isTRUE(.ctFitModelObject(x)$continuoustime)) types <- "kalman"
  }
  unknown <- setdiff(types, available)
  if (length(unknown)) {
    stop("plot types ", paste0("'", unknown, "'", collapse = ", "),
      " need Stan's sample object; backend='julia' offers ",
      paste0("'", available, "'", collapse = ", "), ".", call. = FALSE)
  }
  waitf <- function() {
    if (!isTRUE(wait) || !length(types)) return(TRUE)
    answer <- readline("Input [s] to stop, or leave blank and press [return] for next plot.")
    !answer %in% c("s", "S")
  }
  if ("regression" %in% types) {
    message("Plotting model implied regression coeffcients conditional on time interval using ctDiscretePars")
    print(ctDiscretePars(x, plot = TRUE, ...))
    types <- types[types != "regression"]
    if (!waitf()) return(invisible(NULL))
  }
  if ("kalman" %in% types) {
    message("Plotting expectations from ctPredict")
    print(ctPredict(x, plot = TRUE, ...))
  }
  invisible(NULL)
}

#' @export
coef.ctJuliaFit <- function(object, ...) object$estimate$raw

#' @export
logLik.ctJuliaFit <- function(object, ...) {
  # `object$data` is now the sentinel-cleaned `standata` structure (a list,
  # for parity with a stan fit -- see R/ctFit.R), not a data.frame, so
  # `nrow()` on it no longer gives the row count. `$standata$ndatapoints` is
  # the same count read directly off that structure; `model_spec$data` (the
  # long data.frame `.ctJuliaPrepare()` kept) is the fallback for a fit built
  # before that field existed.
  nobs <- if (!is.null(object$standata$ndatapoints)) as.integer(object$standata$ndatapoints) else
    nrow(.ctBackendSpec(object)$data)
  structure(object$estimate$loglik, df = length(object$estimate$raw), nobs = nobs, class = "logLik")
}
