# Every call into the Julia session goes through here ---------------------------
#
# JuliaConnectoR's `juliaCall()` wraps each call in
#
#     tryCatch(..., interrupt = function(e) killJulia())
#
# so Escape during a call kills the Julia process and then *returns NULL* -- the
# interrupt is caught there and never re-raised. The fit therefore carried on
# past the Escape with a NULL result, started a new Julia, noticed that its
# compiled objective belonged to the dead one, and recompiled the model shape:
# "Stopping Julia ... Starting Julia ... rebuilding ... Compiling", and no way
# back to the prompt short of waiting it out. ctsem's own interrupt handler
# never ran, because the condition never got past JuliaConnectoR's.
#
# So ctsem sends its calls itself, over JuliaConnectoR's own connection and in
# its own wire format, using the same reading and writing functions its
# `juliaCall()` is built from. The only difference is what an interrupt does:
#
#   * Escape while R waits for a reply -- which is where R spends essentially
#     all of any call worth interrupting -- returns to the prompt at once and
#     leaves Julia running. The request is complete and nothing of the reply
#     has been read, so the stream is still in step; the reply is simply owed.
#     R creates the session's interrupt file, and the engine's next
#     per-iteration checkpoint sees it and ends the call (interrupt.jl). Work
#     with no checkpoint in it -- compiling a model shape, precompiling the
#     engine -- runs to completion in the background, and what it compiled is
#     kept.
#   * The next call collects the owed reply first and discards it. If Julia is
#     still busy it says so, and Escape during that wait ends the Julia process
#     instead: that is the one case where waiting is the only alternative.
#   * Sending a request and reading a reply's body are not interruptible, so
#     Escape can never leave half a message on the socket. Both take
#     microseconds to milliseconds; the interrupt is held until they finish.
#
# Julia is ended outright only when it has to be, and by process id rather than
# JuliaConnectoR's netstat/lsof search for the port. A busy Julia is also ended
# when R exits, and ends itself if R disappears without exiting cleanly.
#
# JuliaConnectoR exports none of the pieces this needs. They are resolved
# defensively, as `.ctJuliaCommunicator()` already does for the socket tuning,
# and if any has moved the calls fall back to `JuliaConnectoR::juliaCall()`:
# the old behaviour on Escape, except that the fit is stopped rather than
# silently continued.

# JuliaConnectoR's wire functions, or NULL if any of them cannot be found.
.ctJuliaWire <- function() {
  wire <- .ct_julia_cache$wire
  if (!is.null(wire)) return(if (isFALSE(wire)) NULL else wire)
  wire <- tryCatch({
    ns <- asNamespace("JuliaConnectoR")
    pieces <- c(pkgLocal = "pkgLocal", ensure = "ensureJuliaConnection",
      CALL = "CALL_INDICATOR", RESULT = "RESULT_INDICATOR", FAIL = "FAIL_INDICATOR",
      STDOUT = "STDOUT_INDICATOR", STDERR = "STDERR_INDICATOR",
      writeString = "writeString", writeList = "writeList",
      readElement = "readElement", readString = "readString",
      readOutput = "readOutput", readCall = "readCall",
      answerCallback = "answerCallback", writeFailMessage = "writeFailMessage")
    out <- lapply(pieces, get, envir = ns, inherits = FALSE)
    if (!is.environment(out$pkgLocal)) stop("pkgLocal moved")
    out
  }, error = function(e) FALSE)
  .ct_julia_cache$wire <- wire
  if (isFALSE(wire)) NULL else wire
}

#' @keywords internal
.ctJuliaCall <- function(name, ..., .defer = FALSE) {
  wire <- .ctJuliaWire()
  if (is.null(wire)) return(.ctJuliaCallFallback(name, ...))
  .ctJuliaCheckSwallowed()
  if (!is.null(.ct_julia_cache$inflight)) {
    # A restore from an `on.exit` running while an interrupt unwinds. Its reply
    # is not needed, and waiting for the interrupted call here would hold up
    # the very return to the prompt that Escape asked for, so it goes after
    # the owed reply instead -- the order it would have run in anyway.
    if (isTRUE(.defer)) {
      .ct_julia_cache$deferred <- c(.ct_julia_cache$deferred,
        list(list(name = name, args = list(...))))
      return(invisible(NULL))
    }
    .ctJuliaSettle(wire)
  }
  .ctJuliaEnsureStarted(wire)
  result <- .ctJuliaExchange(wire, name, list(...))
  .ctJuliaReleaseRefs(wire)
  result
}

#' @keywords internal
.ctJuliaEval <- function(expr) .ctJuliaCall("RConnector.mainevalcmd", expr)

# `JuliaConnectoR::juliaPut()`, sent through `.ctJuliaCall()`.
#' @keywords internal
.ctJuliaPut <- function(x) {
  if (inherits(x, "JuliaProxy")) stop("Argument is already a Julia object", call. = FALSE)
  if (is.list(x) && !is.null(attr(x, "JLTYPE"))) return(.ctJuliaCall("identity", x))
  .ctJuliaCall("RConnector.EnforcedProxy", x)
}

# `JuliaConnectoR::juliaGet()`, sent through `.ctJuliaCall()`. The flag that
# asks for a full translation is reset in an `on.exit`, deferred if an
# interrupt left a reply owed, since the next call must not inherit it.
#' @keywords internal
.ctJuliaGet <- function(x) {
  if (!inherits(x, "JuliaProxy")) return(JuliaConnectoR::juliaGet(x))
  wire <- .ctJuliaWire()
  if (is.null(wire)) return(.ctJuliaCallFallback(NULL, x))
  .ctJuliaCall("RConnector.full_translation!", wire$pkgLocal$communicator, TRUE)
  on.exit(.ctJuliaCall("RConnector.full_translation!", wire$pkgLocal$communicator,
    FALSE, .defer = TRUE), add = TRUE)
  .ctJuliaCall("identity", x)
}

# The engine module, with every function sent through `.ctJuliaCall()`.
#
# `juliaImport()` builds an environment of R functions that each call
# `juliaCall()` with the Julia name recorded in their `JLFUN` attribute. The
# import itself is a few short calls, made with interrupts held; what it returns
# is rewrapped so the fitting calls made through it -- the long ones -- are ours.
#' @keywords internal
.ctJuliaImport <- function(module) {
  imported <- suspendInterrupts(.ctJuliaHeld(JuliaConnectoR::juliaImport(module)))
  if (is.null(.ctJuliaWire())) return(imported)
  out <- new.env(parent = emptyenv())
  for (name in ls(imported, all.names = TRUE)) {
    f <- get(name, envir = imported, inherits = FALSE)
    jl <- attr(f, "JLFUN")
    assign(name, if (is.character(jl) && is.null(attr(f, "JLTYPE"))) {
      local({
        target <- jl
        function(...) .ctJuliaCall(target, ...)
      })
    } else f, envir = out)
  }
  out
}

# Writing a request and reading a reply ---------------------------------------

.ctJuliaExchange <- function(wire, name, args) {
  con <- wire$pkgLocal$con
  failed <- FALSE
  withCallingHandlers(suspendInterrupts({
    writeBin(wire$CALL, con)
    wire$writeString(name)
    wire$writeList(args)
  }), warning = function(w) {
    # JuliaConnectoR only warns when a write fails, then waits for a reply
    # that cannot come. Matched on the call rather than the translated text.
    call <- conditionCall(w)
    if (is.call(call) && identical(call[[1L]], as.name("writeBin"))) {
      failed <<- TRUE
      invokeRestart("muffleWarning")
    }
  })
  if (failed) .ctJuliaLost(name)
  .ctJuliaAwait(wire, name)
}

# Read messages until the reply to `name` arrives. `discard` is the settling of
# an interrupted call: its output is dropped, a callback it makes is refused,
# and its result or failure is read and thrown away.
.ctJuliaAwait <- function(wire, name, discard = FALSE, announce = NULL) {
  repeat {
    type <- .ctJuliaNextMessage(wire, name, discard)
    done <- suspendInterrupts(.ctJuliaTakeMessage(wire, type, discard))
    if (is.list(done)) {
      if (identical(done$kind, "result")) return(done$value)
      if (!discard) stop(done$message, call. = FALSE, domain = NA)
      return(invisible(NULL))
    }
    if (!is.null(announce)) announce()
  }
}

# The first byte of the next message, and the one place a call can be
# interrupted. A read that returns nothing is either the socket's timeout or
# the end of the stream; a socket that is readable yet yields nothing is the
# end, which is how a Julia that has died is told apart from one that is busy.
# JuliaConnectoR's own reader treats both as "nothing yet" and spins at full CPU
# for the rest of the R session.
.ctJuliaNextMessage <- function(wire, name, discard) {
  con <- wire$pkgLocal$con
  byte <- withCallingHandlers({
    repeat {
      byte <- readBin(con, "raw", 1L)
      if (length(byte)) break
      if (isTRUE(tryCatch(socketSelect(list(con), timeout = 0), error = function(e) TRUE))) {
        byte <- readBin(con, "raw", 1L)
        if (length(byte)) break
        byte <- NULL
        break
      }
    }
    byte
  }, interrupt = function(cnd) {
    if (discard) {
      .ctJuliaKill()
      message("Julia stopped; the next call starts a new session.")
    } else {
      .ctJuliaAbandon(name)
    }
  })
  if (is.null(byte)) .ctJuliaLost(name)
  byte
}

.ctJuliaTakeMessage <- function(wire, type, discard) {
  if (identical(type, wire$RESULT)) {
    value <- wire$readElement()
    return(list(kind = "result", value = if (discard) NULL else value))
  }
  if (identical(type, wire$FAIL)) return(list(kind = "fail", message = wire$readString()))
  if (identical(type, wire$STDOUT) || identical(type, wire$STDERR)) {
    # An output message is framed as a string is: a length, then the bytes.
    if (discard) wire$readString()
    else wire$readOutput(writeTo = if (identical(type, wire$STDOUT)) stdout() else stderr())
    return(NULL)
  }
  if (identical(type, wire$CALL)) {
    call <- wire$readCall()
    if (discard) {
      # Refused rather than answered: the progress display or callback it is
      # for belongs to a fit that is gone. The engine's callers disable a
      # callback that fails and carry on to their next checkpoint.
      wire$writeFailMessage("the call was interrupted from R")
    } else {
      wire$answerCallback(get(call$name, envir = wire$pkgLocal$callbacks), call$args)
    }
    return(NULL)
  }
  .ctJuliaKill()
  stop("The Julia session sent a message ctsem could not read, so it was stopped. ",
    "The next call starts a new one.", call. = FALSE)
}

# Releasing proxies R has collected, as `juliaCall()` does after every call.
# The queue is one raw vector of 8-byte references, and a collection during the
# call can append to it, so only what was sent is removed.
.ctJuliaReleaseRefs <- function(wire) {
  local <- wire$pkgLocal
  refs <- local$finalizedRefs
  if (is.null(refs)) return(invisible(NULL))
  sent <- length(refs)
  try({
    ids <- .ctJuliaExchange(wire, "RConnector.decrefcounts",
      list(local$communicator, refs))
    if (length(ids)) rm(envir = local$callbacks, list = ids)
  }, silent = TRUE)
  rest <- local$finalizedRefs
  local$finalizedRefs <- if (length(rest) > sent) rest[-seq_len(sent)] else NULL
  invisible(NULL)
}

# Interrupts ------------------------------------------------------------------

# Escape while a reply is owed: ask the engine to stop, remember the reply, and
# let the interrupt carry on to the prompt.
.ctJuliaAbandon <- function(name) {
  .ct_julia_cache$inflight <- list(name = name, since = Sys.time())
  file <- .ct_julia_cache$interrupt_file
  if (!is.null(file)) try(file.create(file, showWarnings = FALSE), silent = TRUE)
  invisible(NULL)
}

# Collect an interrupted call's reply before anything else is sent.
.ctJuliaSettle <- function(wire) {
  inflight <- .ct_julia_cache$inflight
  con <- wire$pkgLocal$con
  if (is.null(con)) {
    .ct_julia_cache$inflight <- NULL
    .ct_julia_cache$deferred <- NULL
    return(invisible(FALSE))
  }
  said <- FALSE
  announce <- function() {
    if (said || difftime(Sys.time(), started, units = "secs") < 1) return()
    said <<- TRUE
    message("Waiting for Julia to finish the interrupted call; Esc stops Julia instead.")
  }
  started <- Sys.time()
  if (!isTRUE(tryCatch(socketSelect(list(con), timeout = 1), error = function(e) TRUE))) {
    announce()
  }
  .ctJuliaAwait(wire, inflight$name, discard = TRUE, announce = announce)
  .ct_julia_cache$inflight <- NULL
  file <- .ct_julia_cache$interrupt_file
  if (!is.null(file)) unlink(file)
  deferred <- .ct_julia_cache$deferred
  .ct_julia_cache$deferred <- NULL
  for (d in deferred) try(.ctJuliaExchange(wire, d$name, d$args), silent = TRUE)
  invisible(TRUE)
}

# JuliaConnectoR ends Julia and swallows the interrupt when Escape lands inside
# one of its own calls -- the fallback route, and proxy methods such as `$` on
# a returned struct. Its message is the only trace, so it is noted here and the
# operation stopped at the next call rather than continued on a new session.
.ctJuliaNoteSwallowed <- function(m) {
  if (startsWith(conditionMessage(m), "Stopping Julia")) {
    .ct_julia_cache$swallowed <- TRUE
    invokeRestart("muffleMessage")
  }
}

.ctJuliaCheckSwallowed <- function() {
  if (!isTRUE(.ct_julia_cache$swallowed)) return(invisible(FALSE))
  .ct_julia_cache$swallowed <- NULL
  .ctJuliaForgetSession()
  .ctJuliaInterruptNow()
}

# What R itself does on Escape: offer the interrupt to any handler, then return
# to the prompt.
.ctJuliaInterruptNow <- function() {
  signalCondition(structure(class = c("interrupt", "condition"),
    list(message = "", call = NULL)))
  invokeRestart("abort")
}

.ctJuliaHeld <- function(expr) withCallingHandlers(expr, message = .ctJuliaNoteSwallowed)

.ctJuliaCallFallback <- function(name, ...) {
  .ctJuliaCheckSwallowed()
  result <- .ctJuliaHeld(if (is.null(name)) JuliaConnectoR::juliaGet(...) else
    JuliaConnectoR::juliaCall(name, ...))
  .ctJuliaCheckSwallowed()
  result
}

# Ending a Julia session ------------------------------------------------------

# By process id, which the session reports when it starts. The connection is
# closed afterwards; JuliaConnectoR's goodbye into a dead socket is harmless.
.ctJuliaKill <- function() {
  pid <- .ct_julia_cache$pid
  .ctJuliaForgetSession()
  if (!is.null(pid) && !is.na(pid)) {
    try(tools::pskill(pid, tools::SIGKILL), silent = TRUE)
  }
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    suppressWarnings(try(JuliaConnectoR::stopJulia(), silent = TRUE))
  }
  invisible(NULL)
}

# The socket ended mid-call: Julia has gone.
.ctJuliaLost <- function(name) {
  .ctJuliaForgetSession()
  if (requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    suppressWarnings(try(JuliaConnectoR::stopJulia(), silent = TRUE))
  }
  stop("The Julia session ended while running ", sub("^.*[.:]", "", name),
    "; its process is no longer running. Calling again starts a new one.",
    call. = FALSE)
}

# Julia ends when R does -- JuliaConnectoR's own exit finalizer closes the
# socket -- but only once it next reads that socket, which a busy Julia does
# when its call finishes. Registered once per R session, on the cache itself.
.ctJuliaRegisterExit <- function() {
  if (isTRUE(.ct_julia_cache$exit_registered)) return(invisible(NULL))
  reg.finalizer(.ct_julia_cache, function(env) {
    if (!is.null(env$inflight) && !is.null(env$pid)) {
      try(tools::pskill(env$pid, tools::SIGKILL), silent = TRUE)
    }
  }, onexit = TRUE)
  .ct_julia_cache$exit_registered <- TRUE
  invisible(NULL)
}

# Starting a Julia session ----------------------------------------------------

# Start one if none is running, saying which Julia it is.
#
# Held against interrupts: JuliaConnectoR launches the process and then polls
# for the port it listens on, and Escape in that loop left the process waiting
# for a connection that would never come, for good. It takes a few seconds.
.ctJuliaEnsureStarted <- function(wire) {
  started <- !is.null(wire$pkgLocal$con)
  # A session something else started -- JuliaConnectoR::juliaSetupOk() starts
  # one -- is adopted once: its process id and interrupt file are still needed.
  if (started && !is.null(.ct_julia_cache$interrupt_file)) return(invisible(FALSE))
  external <- nzchar(Sys.getenv("JULIACONNECTOR_SERVER", unset = ""))
  if (!started) {
    if (!external) .ctJuliaAnnounce()
    suspendInterrupts(withCallingHandlers(wire$ensure(), message = function(m) {
      if (!external && startsWith(conditionMessage(m), "Starting Julia")) {
        invokeRestart("muffleMessage")
      }
    }))
  }
  .ct_julia_cache$pid <- if (external) NULL else tryCatch(
    as.integer(.ctJuliaExchange(wire, "RConnector.mainevalcmd", list("getpid()"))),
    error = function(e) NULL)
  .ct_julia_cache$interrupt_file <- gsub("\\", "/",
    tempfile("ctsem-julia-interrupt-"), fixed = TRUE)
  .ctJuliaRegisterExit()
  invisible(TRUE)
}

# "Starting Julia 1.12.5 with 2 threads ...", in place of JuliaConnectoR's
# "Starting Julia ...", and a line when a newer Julia worth having is known to
# exist. Offline throughout: the version is the binary's own, and what is newer
# comes from juliaup's cached release list and from the version ctsem pins.
.ctJuliaAnnounce <- function() {
  bin <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (!.ctJuliaIsBinDir(bin)) {
    bin <- tryCatch(.ctJuliaBin(), error = function(e) NULL)
    if (!is.null(bin)) Sys.setenv(JULIA_BINDIR = bin)
  }
  version <- if (is.null(bin)) NA_character_ else .ctJuliaCachedVersion(bin)
  threads <- .ctJuliaStartThreads()
  message("Starting Julia",
    if (!is.na(version)) paste0(" ", version) else "",
    if (!is.na(threads)) paste0(" with ", threads,
      if (identical(threads, "1")) " thread" else " threads") else "",
    " ...")
  newer <- .ctJuliaNewerRelease(version, bin)
  if (!is.null(newer) && !isTRUE(.ct_julia_cache$update_said)) {
    .ct_julia_cache$update_said <- TRUE
    message("  Julia ", newer$version, " is available", newer$how, ".")
  }
  invisible(NULL)
}

.ctJuliaCachedVersion <- function(bin) {
  known <- .ct_julia_cache$bin_versions
  if (!is.null(known[[bin]])) return(known[[bin]])
  version <- .ctJuliaBinVersion(bin)
  .ct_julia_cache$bin_versions <- c(known, stats::setNames(list(version), bin))
  version
}

# The thread count the process will start with, as Julia will read it.
.ctJuliaStartThreads <- function() {
  value <- trimws(Sys.getenv("JULIA_NUM_THREADS", unset = ""))
  if (!nzchar(value)) return("1")
  first <- sub(",.*$", "", value)
  if (identical(first, "auto")) return("auto")
  if (grepl("^[0-9]+$", first)) return(first)
  NA_character_
}

# The newest Julia worth offering over `version`, known without a network.
#
# A newer patch of the running series is always worth it: it cannot move the
# engine's resolved manifest. A newer minor series only when it is the series
# ctsem is tested on, `.ct_julia_version`'s -- a later one may well work, but
# nothing here has checked. The candidates are that pinned version itself,
# which is known to exist, and whatever juliaup's release list says is newest
# in either series; a Julia ctsem installed is offered only the pin, since
# that is what ctJuliaInstall() installs.
.ctJuliaNewerRelease <- function(version, bin) {
  if (length(version) != 1L || is.na(version) || is.null(bin)) return(NULL)
  current <- tryCatch(numeric_version(version), error = function(e) NULL)
  if (is.null(current)) return(NULL)
  series <- function(v) sub("^(\\d+\\.\\d+).*$", "\\1", v)
  running <- series(version)
  tested <- series(.ct_julia_version)
  managed <- startsWith(normalizePath(bin, winslash = "/", mustWork = FALSE),
    .ctJuliaInstallRoot())
  candidates <- .ct_julia_version
  if (!managed) {
    candidates <- c(candidates, .ctJuliaupLatest(running),
      if (numeric_version(tested) > numeric_version(running)) .ctJuliaupLatest(tested))
  }
  candidates <- candidates[!is.na(candidates)]
  keep <- vapply(candidates, function(v) series(v) %in% c(running, tested) &&
    numeric_version(v) > current, logical(1))
  candidates <- candidates[keep]
  if (!length(candidates)) return(NULL)
  best <- candidates[order(numeric_version(candidates), decreasing = TRUE)][[1L]]
  # No instruction for juliaup: `juliaup update` moves a channel to its newest
  # release, which for the default channel can be a series ctsem is not tested
  # on.
  how <- if (managed) "; ctJuliaInstall(force = TRUE) installs it" else
    if (series(best) != running) paste0(", the series ctsem is tested on") else ""
  list(version = best, how = how)
}

# The newest release in a minor series according to juliaup's version database,
# which juliaup itself keeps current in the background. NA when there is none.
.ctJuliaupLatest <- function(series) {
  roots <- unique(c(Sys.getenv("JULIAUP_DEPOT_PATH", unset = ""),
    file.path(c(Sys.getenv("USERPROFILE", unset = ""), Sys.getenv("HOME", unset = ""),
      path.expand("~")), ".julia")))
  roots <- roots[nzchar(roots)]
  files <- unlist(lapply(file.path(roots, "juliaup"), function(dir)
    list.files(dir, pattern = "^versiondb-.*\\.json$", full.names = TRUE)))
  if (!length(files)) return(NA_character_)
  text <- tryCatch(paste(readLines(files[[1L]], warn = FALSE), collapse = ""),
    error = function(e) "")
  # The channel named after the series, e.g. "1.12": {"Version": "1.12.7+0..."}.
  pattern <- paste0('"', gsub(".", "\\.", series, fixed = TRUE),
    '"\\s*:\\s*\\{\\s*"Version"\\s*:\\s*"(\\d+\\.\\d+\\.\\d+)')
  hit <- regmatches(text, regexec(pattern, text))[[1L]]
  if (length(hit) < 2L) NA_character_ else hit[[2L]]
}
