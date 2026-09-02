# One-call setup for backend='julia' -----------------------------------------
#
# The Julia backend needs three things before it can run: the JuliaConnectoR
# bridge package, a Julia installation, and the vendored engine's Julia
# dependencies. Only the third was ever automatic. The other two each failed
# with an error telling the user to go away and do something -- install a
# package, download Julia, set JULIA_BINDIR, restart R -- so trying the backend
# for the first time cost four manual steps and two R sessions.
#
# All three are one call now, ctJuliaInstall(), and the two that can be missing
# offer themselves at the point of failure, so ctFit(backend='julia') on a bare
# machine is a question rather than an error. Nothing needs a restart: a package
# installed into the running session's library is visible immediately, and a
# Julia that ctsem installed lands in a directory .ctJuliaBin() searches, so
# later sessions find it with no environment variable to remember.
#
# Nothing here touches the network without consent. An interactive session is
# asked, with the download size and the destination named; a non-interactive one
# declines unless agree=TRUE or CTSEM_JULIA_AGREE=yes said otherwise in advance.
# That is what keeps R CMD check, and CRAN's own runs, offline.

# The Julia ctsem installs when it has to install one. Pinned rather than
# resolved from julialang.org's version index at run time, for two reasons: it
# is what allows the archive checksums below to ship with the package, and it is
# the version the vendored Manifest.toml was resolved under, so the engine's
# dependencies instantiate exactly as they were tested rather than being
# re-resolved on first use.
.ct_julia_version <- "1.12.7"

# sha256 and approximate size of the official archive for .ct_julia_version, per
# platform. A user-supplied version is downloaded unverified -- with a message
# saying so -- since these are the only hashes we know.
.ct_julia_archives <- list(
  `windows-x86_64` = list(dir = "winnt/x64", suffix = "win64", ext = "zip", mb = 267,
    sha256 = "ff5c7eb354c2fcb48401114a5fbcfe8e60181f95d9af42b266f975265a5bad47"),
  `linux-x86_64` = list(dir = "linux/x64", suffix = "linux-x86_64", ext = "tar.gz", mb = 277,
    sha256 = "4e7e9e776634d24835250de67cde39b0d4af15bc432eb20697e6be6c28ea69e8"),
  `linux-aarch64` = list(dir = "linux/aarch64", suffix = "linux-aarch64", ext = "tar.gz", mb = 293,
    sha256 = "9243c0b524c7f300883240a1ee5ea3916a30e070bff718acf8ccaee31a731ef2"),
  `macos-x86_64` = list(dir = "mac/x64", suffix = "mac64", ext = "tar.gz", mb = 259,
    sha256 = "a21a15c7b7d294a03482a3598b18cde48d37be86ac7a408e495e51fb3afe0157"),
  `macos-aarch64` = list(dir = "mac/aarch64", suffix = "macaarch64", ext = "tar.gz", mb = 220,
    sha256 = "af8fcedfb25b6b9c8c13d99695c38faa2d59bf0162b3f458c8b8b56f14c96919"))

# The oldest Julia the engine is known to work on. A Julia found on the user's
# machine that is older than this is treated as absent rather than used and
# failed on later, deep inside Pkg.
.ct_julia_minimum <- "1.10"

.ctJuliaExeName <- function() if (.Platform$OS.type == "windows") "julia.exe" else "julia"

.ctJuliaIsBinDir <- function(dir) {
  length(dir) == 1L && !is.na(dir) && nzchar(dir) &&
    file.exists(file.path(dir, .ctJuliaExeName()))
}

# Newest of a set of julia-<version>[...] directories. Sorting matters: the
# names are alphabetical only by accident, and "julia-1.9" sorts after
# "julia-1.12" as a string.
.ctJuliaNewestDir <- function(dirs) {
  if (!length(dirs)) return(NULL)
  named <- grepl("julia-\\d+(\\.\\d+)*", basename(dirs))
  if (!any(named)) return(dirs[[length(dirs)]])
  dirs <- dirs[named]
  version <- sub("^.*julia-(\\d+(\\.\\d+)*).*$", "\\1", basename(dirs))
  dirs[[order(numeric_version(version), decreasing = TRUE)[[1]]]]
}

# Where a ctsem-installed Julia lives: under R_user_dir(), never in the R
# library, the home directory, or anywhere else CRAN policy puts off limits.
.ctJuliaInstallRoot <- function() {
  gsub("\\", "/", file.path(tools::R_user_dir("ctsem", which = "data"), "julia"), fixed = TRUE)
}

.ctJuliaManagedBin <- function() {
  root <- .ctJuliaInstallRoot()
  if (!dir.exists(root)) return(NULL)
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[vapply(file.path(dirs, "bin"), .ctJuliaIsBinDir, logical(1))]
  newest <- .ctJuliaNewestDir(dirs)
  if (is.null(newest)) return(NULL)
  normalizePath(file.path(newest, "bin"), winslash = "/", mustWork = TRUE)
}

# juliaup's own installs, for a user who has juliaup but whose R session did not
# inherit its PATH -- which is the normal case for R started from a launcher
# rather than from a shell.
.ctJuliaJuliaupBin <- function() {
  roots <- unique(c(Sys.getenv("USERPROFILE", unset = ""), Sys.getenv("HOME", unset = ""),
    path.expand("~")))
  for (root in roots[nzchar(roots)]) {
    for (shim in c(file.path(root, ".juliaup", "bin"),
      file.path(root, ".julia", "juliaup", "bin"))) {
      if (.ctJuliaIsBinDir(shim)) return(normalizePath(shim, winslash = "/", mustWork = TRUE))
    }
    installs <- list.dirs(file.path(root, ".julia", "juliaup"), recursive = FALSE,
      full.names = TRUE)
    installs <- installs[vapply(file.path(installs, "bin"), .ctJuliaIsBinDir, logical(1))]
    newest <- .ctJuliaNewestDir(installs)
    if (!is.null(newest)) {
      return(normalizePath(file.path(newest, "bin"), winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

# The version a Julia binary directory actually reports, or NA. Asking the
# binary rather than reading its directory name is the only reliable form: a
# juliaup shim's directory name says nothing about what it launches.
.ctJuliaBinVersion <- function(bin) {
  if (!.ctJuliaIsBinDir(bin)) return(NA_character_)
  out <- tryCatch(suppressWarnings(
    system2(file.path(bin, .ctJuliaExeName()), "--version", stdout = TRUE, stderr = FALSE,
      timeout = 120)), error = function(e) character())
  hit <- regmatches(out, regexpr("\\d+\\.\\d+\\.\\d+", out))
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

.ctJuliaVersionOk <- function(version) {
  length(version) == 1L && !is.na(version) &&
    numeric_version(version) >= numeric_version(.ct_julia_minimum)
}

.ctJuliaPlatform <- function() {
  os <- if (.Platform$OS.type == "windows") "windows" else
    switch(Sys.info()[["sysname"]], Darwin = "macos", Linux = "linux", NA_character_)
  arch <- R.version$arch
  if (arch %in% c("arm64", "aarch64")) arch <- "aarch64"
  if (arch %in% c("x86_64", "x86-64", "amd64")) arch <- "x86_64"
  if (is.na(os) || !arch %in% c("x86_64", "aarch64")) return(NA_character_)
  paste0(os, "-", arch)
}

# URL, checksum and size of the official archive for one version and platform.
.ctJuliaArchive <- function(version = .ct_julia_version, platform = .ctJuliaPlatform()) {
  if (is.na(platform) || is.null(.ct_julia_archives[[platform]])) return(NULL)
  spec <- .ct_julia_archives[[platform]]
  minor <- sub("^(\\d+\\.\\d+).*$", "\\1", version)
  list(
    url = sprintf("https://julialang-s3.julialang.org/bin/%s/%s/julia-%s-%s.%s",
      spec$dir, minor, version, spec$suffix, spec$ext),
    ext = spec$ext, mb = spec$mb,
    sha256 = if (identical(version, .ct_julia_version)) spec$sha256 else NA_character_)
}

# Consent. isTRUE/isFALSE(agree) settle it outright; otherwise the environment
# variable, and otherwise the user -- but only if there is a user to ask.
.ctJuliaAgreed <- function(agree, prompt) {
  if (isTRUE(agree)) return(TRUE)
  if (isFALSE(agree)) return(FALSE)
  configured <- tolower(trimws(Sys.getenv("CTSEM_JULIA_AGREE", unset = "")))
  if (configured %in% c("yes", "true", "1")) return(TRUE)
  if (configured %in% c("no", "false", "0")) return(FALSE)
  if (!interactive()) return(FALSE)
  message(prompt)
  isTRUE(utils::askYesNo("Proceed?", default = TRUE))
}

# A library we may install into. install.packages() would otherwise fail, or
# raise a second prompt of its own, when .libPaths()[1] is the system library.
.ctJuliaWritableLib <- function() {
  first <- .libPaths()[[1]]
  if (dir.exists(first) && file.access(first, mode = 2) == 0L) return(first)
  personal <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(personal)) return(first)
  personal <- strsplit(personal, .Platform$path.sep, fixed = TRUE)[[1]][[1]]
  dir.create(personal, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(personal)) return(first)
  if (!personal %in% .libPaths()) .libPaths(c(personal, .libPaths()))
  personal
}

.ctJuliaInstallConnectoR <- function(agree = NULL, quiet = FALSE) {
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) return(TRUE)
  lib <- .ctJuliaWritableLib()
  agreed <- .ctJuliaAgreed(agree, paste0(
    "The julia backend needs the JuliaConnectoR package, which is not installed.\n",
    "  install: JuliaConnectoR, from CRAN (pure R, nothing to compile)\n",
    "  into:    ", lib))
  if (!agreed) return(FALSE)
  utils::install.packages("JuliaConnectoR", lib = lib, quiet = quiet)
  requireNamespace("JuliaConnectoR", quietly = TRUE)
}

.ctJuliaChecksumOk <- function(path, sha256) {
  if (is.na(sha256)) {
    message("No checksum is recorded for this Julia version, so it was not verified.")
    return(TRUE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    message("Package 'digest' is not installed, so the download could not be verified.")
    return(TRUE)
  }
  identical(digest::digest(path, algo = "sha256", file = TRUE), sha256)
}

# Download and unpack an official Julia into .ctJuliaInstallRoot(). Returns its
# bin directory, or NULL if consent was refused.
.ctJuliaInstallJulia <- function(version = NULL, agree = NULL, quiet = FALSE) {
  version <- .ctJuliaOr(version, Sys.getenv("CTSEM_JULIA_VERSION", unset = .ct_julia_version))
  archive <- .ctJuliaArchive(version)
  if (is.null(archive)) {
    stop("ctsem has no Julia download for this platform (", R.version$platform, ").\n",
      "  Install Julia (", .ct_julia_minimum, " or newer) from https://julialang.org/downloads/,\n",
      "  then either put its bin directory on PATH or set\n",
      "    Sys.setenv(JULIA_BINDIR = \"/path/to/julia/bin\")", call. = FALSE)
  }
  target <- file.path(.ctJuliaInstallRoot(), paste0("julia-", version))
  agreed <- .ctJuliaAgreed(agree, paste0(
    "The julia backend needs a Julia installation, and none was found.\n",
    "  download: ", archive$url, "\n",
    "            (the official build from julialang.org, about ", archive$mb, " MB)\n",
    "  into:     ", target, "\n",
    "  Nothing outside that directory is modified: no PATH, no shell profile."))
  if (!agreed) return(NULL)

  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  tarball <- file.path(tempdir(), basename(archive$url))
  on.exit(unlink(tarball), add = TRUE)
  # R's default 60 s applies to the whole transfer rather than to inactivity,
  # and this is a quarter of a gigabyte.
  timeout <- options(timeout = max(3600L, getOption("timeout")))
  on.exit(options(timeout), add = TRUE)
  if (!quiet) message("Downloading Julia ", version, " (~", archive$mb, " MB)...")
  status <- tryCatch(utils::download.file(archive$url, tarball, mode = "wb", quiet = quiet),
    error = function(e) e)
  # download.file() reports failure by return code as well as by condition,
  # depending on the method in use, and a failed transfer can still leave a
  # partial file behind.
  failed <- if (inherits(status, "error")) conditionMessage(status) else
    if (!identical(as.integer(status), 0L)) paste("download.file() returned", status) else
      if (!file.exists(tarball)) "the download produced no file" else NULL
  if (!is.null(failed)) {
    stop("Could not download Julia from ", archive$url, "\n  ", failed, call. = FALSE)
  }
  if (!.ctJuliaChecksumOk(tarball, archive$sha256)) {
    stop("The downloaded Julia archive does not match the sha256 checksum recorded in ctsem, ",
      "so it was not installed.\n  This usually means a truncated download; try again.",
      call. = FALSE)
  }

  if (!quiet) message("Unpacking...")
  .ctJuliaUnpack(tarball, target, archive$ext)
}

# Unpack a Julia archive into <target>, and return <target>/bin.
#
# Separate from the download so the fiddly half can be tested without a quarter
# of a gigabyte of network: everything that can go wrong here -- an archive
# whose contents are wrapped in a directory, an interrupted extraction, a lost
# executable bit -- goes wrong the same way for a two-file archive.
.ctJuliaUnpack <- function(archive_file, target, ext) {
  # Extract into a staging directory and rename, so an interrupted extraction
  # cannot leave a half-populated install that a later session would find and
  # believe in.
  staging <- paste0(target, "-partial")
  unlink(staging, recursive = TRUE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  if (identical(ext, "zip")) {
    utils::unzip(archive_file, exdir = staging)
  } else {
    utils::untar(archive_file, exdir = staging)
  }
  # The archives hold a single julia-<version>/ directory; hoist it so the
  # layout is <target>/bin/julia either way.
  inner <- list.dirs(staging, recursive = FALSE, full.names = TRUE)
  source <- if (length(inner) == 1L && !dir.exists(file.path(staging, "bin"))) {
    inner[[1]]
  } else staging
  unlink(target, recursive = TRUE)
  if (!file.rename(source, target)) {
    stop("Could not move the unpacked Julia into ", target, call. = FALSE)
  }
  # Not every extraction path preserves the executable bit, and a Julia that
  # cannot be executed fails later with a message about nothing in particular.
  if (.Platform$OS.type != "windows") {
    for (dir in c("bin", "libexec")) {
      files <- list.files(file.path(target, dir), full.names = TRUE, recursive = TRUE)
      if (length(files)) Sys.chmod(files, "0755", use_umask = FALSE)
    }
  }

  bin <- file.path(target, "bin")
  if (!.ctJuliaIsBinDir(bin)) {
    stop("The unpacked Julia does not contain ", file.path(bin, .ctJuliaExeName()), call. = FALSE)
  }
  normalizePath(bin, winslash = "/", mustWork = TRUE)
}

# Offered at the point of failure by .ctJuliaCheckAvailable(), so that the first
# thing a user does with the backend is the thing they meant to do rather than a
# setup errand.
.ctJuliaOfferJulia <- function() {
  bin <- tryCatch(.ctJuliaInstallJulia(), error = function(e) {
    message(conditionMessage(e))
    NULL
  })
  if (is.null(bin)) return(FALSE)
  Sys.setenv(JULIA_BINDIR = bin)
  TRUE
}

# Pre-flight for ctFit(backend='julia'): make sure the two things that can be
# missing are present, offering to supply them, and do it before the fit spends
# any time preparing data -- a user who declines should not have waited first.
#
# Deliberately does not start a Julia session, only checks that one could be:
# .ctFitJuliaBackend() still has to translate `cores` into JULIA_NUM_THREADS
# before the process exists, which is the only moment Julia will read it.
.ctJuliaEnsureInstalled <- function() {
  .ctJuliaRequire()
  if (is.null(.ctJuliaBin()) && !isTRUE(.ctJuliaOfferJulia())) {
    stop("Julia was not found, and backend='julia' needs it (", .ct_julia_minimum,
      " or newer).\n", .ctJuliaDeclined(), call. = FALSE)
  }
  invisible(TRUE)
}

# What to tell a user who was not asked, or who said no.
.ctJuliaDeclined <- function() {
  if (interactive()) {
    "  Run ctJuliaInstall() to set it up, or use backend='stan', which needs no external toolchain."
  } else {
    paste0("  Run ctJuliaInstall() in an interactive session, or ctJuliaInstall(agree = TRUE) ",
      "to consent from a script.\n",
      "  backend='stan' needs no external toolchain.")
  }
}

#' Set up the Julia backend
#'
#' Installs everything \code{ctFit(backend='julia')} needs, in one call and
#' without restarting R: the \pkg{JuliaConnectoR} bridge package, a Julia
#' installation, and the Julia dependencies of the engine that ships inside
#' ctsem. Any step that is already satisfied is skipped, so on a machine that
#' already has Julia this only prepares the engine.
#'
#' Each step that installs something asks first, naming the download and its
#' destination. A non-interactive session has nobody to ask, so nothing is
#' installed there unless \code{agree = TRUE} or the environment variable
#' \code{CTSEM_JULIA_AGREE=yes} gave permission in advance.
#'
#' Julia is downloaded only if none can be found. ctsem looks in
#' \code{JULIA_BINDIR}, then at a Julia it installed itself, then on the
#' \code{PATH}, then among juliaup's installations. A Julia that ctsem installs
#' goes under \code{tools::R_user_dir("ctsem", "data")} and later sessions find
#' it automatically -- there is no environment variable to set, and no PATH or
#' shell profile is modified. To use a different Julia afterwards, set
#' \code{JULIA_BINDIR}.
#'
#' @param threads Number of Julia threads for the engine's subject loop. Julia
#'   fixes this when its process starts, so it applies to the session this call
#'   creates; see \code{\link{ctJuliaSetup}}.
#' @param version Julia version to install. Defaults to the version ctsem's
#'   vendored engine was resolved against. The download is checksum-verified
#'   only for that default.
#' @param agree \code{TRUE} to consent to the installs without being asked,
#'   \code{FALSE} to refuse them. \code{NULL}, the default, asks in an
#'   interactive session and refuses otherwise.
#' @param force Reinstall Julia even if one is already available, and restart
#'   any running Julia session.
#' @param quiet Suppress download and installation progress output.
#' @return A Julia-engine status list, invisibly; see \code{\link{ctJuliaStatus}}.
#' @examples
#' \dontrun{
#' ctJuliaInstall()            # ask, then install whatever is missing
#' ctJuliaInstall(threads = 4) # and start the Julia session with four threads
#' }
#' @seealso \code{\link{ctJuliaSetup}}, \code{\link{ctJuliaStatus}}
#' @export
ctJuliaInstall <- function(threads = NULL, version = NULL, agree = NULL, force = FALSE,
  quiet = FALSE) {
  done <- character()
  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    if (!.ctJuliaInstallConnectoR(agree = agree, quiet = quiet)) {
      stop("The julia backend needs the JuliaConnectoR package.\n", .ctJuliaDeclined(),
        call. = FALSE)
    }
    done <- c(done, "JuliaConnectoR")
  }

  bin <- if (isTRUE(force)) NULL else .ctJuliaBin()
  if (!is.null(bin) && !.ctJuliaVersionOk(.ctJuliaBinVersion(bin))) bin <- NULL
  if (is.null(bin)) {
    bin <- .ctJuliaInstallJulia(version = version, agree = agree, quiet = quiet)
    if (is.null(bin)) {
      stop("The julia backend needs a Julia installation.\n", .ctJuliaDeclined(), call. = FALSE)
    }
    done <- c(done, paste0("Julia ", .ctJuliaOr(.ctJuliaBinVersion(bin), "")))
  }
  Sys.setenv(JULIA_BINDIR = bin)

  if (!quiet) {
    message("Preparing the ctsem Julia engine",
      if (!length(done)) " (Julia is already installed)" else "", "...")
  }
  status <- ctJuliaSetup(threads = threads, force = force)
  if (!quiet) {
    if (length(done)) message("Installed: ", paste(done, collapse = ", "), ".")
    message("The julia backend is ready. Use ctFit(..., backend = 'julia').")
  }
  invisible(status)
}
