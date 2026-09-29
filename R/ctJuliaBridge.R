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
      writeString = "writeString", writeList = "writeList", writeElement = "writeElement",
      readElement = "readElement", readCall = "readCall")
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
  # Evaluated before anything is written. An argument is often a Julia call of
  # its own -- `.ctJuliaVector()` marshals through `.ctJuliaPut()` -- and left
  # as a promise it ran when the argument list was being written, putting its
  # request inside this one's. Both ends then waited for each other for good.
  args <- list(...)
  if (!is.null(.ct_julia_cache$inflight)) {
    # A restore from an `on.exit` running while an interrupt unwinds. Its reply
    # is not needed, and waiting for the interrupted call here would hold up
    # the very return to the prompt that Escape asked for, so it goes after
    # the owed reply instead -- the order it would have run in anyway.
    if (isTRUE(.defer)) {
      .ct_julia_cache$deferred <- c(.ct_julia_cache$deferred,
        list(list(name = name, args = args)))
      return(invisible(NULL))
    }
    .ctJuliaSettle(wire)
  }
  .ctJuliaEnsureStarted(wire)
  result <- .ctJuliaExchange(wire, name, args)
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
  force(args)
  request <- .ctJuliaBuffered(wire, function(buffer) {
    writeBin(wire$CALL, buffer)
    wire$writeString(name)
    wire$writeList(args)
  })
  .ctJuliaSend(wire, request, name)
  .ctJuliaAwait(wire, name)
}

# A whole message as bytes, written nowhere yet.
#
# JuliaConnectoR's writers write to `pkgLocal$con` a piece at a time, so a
# message used to reach the socket in pieces: an argument it could not
# translate stopped the writing halfway, and Julia then read the next request
# as the rest of this one. `write` runs with them pointed at a buffer instead,
# put back however it ends. Nothing here touches the socket, so it needs no
# protection from Escape either. `write` must not call Julia.
.ctJuliaBuffered <- function(wire, write) {
  local <- wire$pkgLocal
  socket <- local$con
  buffer <- rawConnection(raw(0), open = "wb")
  local$con <- buffer
  on.exit({
    local$con <- socket
    close(buffer)
  }, add = TRUE)
  write(buffer)
  rawConnectionValue(buffer)
}

# Write a message to Julia in one piece, or say by name that Julia has gone.
#
# A write to a peer that has exited fails differently by platform. Windows
# warns ("problem writing to connection"). Linux raises SIGPIPE, which R turns
# into an *error* ("ignoring SIGPIPE signal") on the second write after the
# peer went -- the kernel accepts the first -- and only a third warns. Both are
# the same fact, and the error, escaping as an ordinary one, left the dead
# connection in place for the next call to write into again. A write the
# kernel does accept surfaces as the end of the stream at the next read, which
# `.ctJuliaNextMessage()` treats the same way.
.ctJuliaSend <- function(wire, bytes, name) {
  failed <- FALSE
  tryCatch(withCallingHandlers(suspendInterrupts(writeBin(bytes, wire$pkgLocal$con)),
    warning = function(w) {
      failed <<- TRUE
      invokeRestart("muffleWarning")
    }), error = function(e) failed <<- TRUE)
  if (failed) .ctJuliaLost(name)
  invisible(NULL)
}

# Answer a callback Julia made during `name`: run it, then send its result, or
# its failure, as one message. Refused outright with `refuse`, the reason. The
# callback runs before anything is buffered because it may call Julia itself;
# a result that cannot be translated is sent as a failure, since Julia is
# waiting for an answer of some kind.
.ctJuliaAnswer <- function(wire, call, name, refuse = NULL) {
  failure <- refuse
  result <- NULL
  if (is.null(failure)) {
    fun <- get(call$name, envir = wire$pkgLocal$callbacks)
    result <- tryCatch(do.call(fun, call$args),
      error = function(e) { failure <<- as.character(e); NULL })
  }
  fail <- function(message) .ctJuliaBuffered(wire, function(buffer) {
    writeBin(wire$FAIL, buffer)
    wire$writeString(message)
  })
  bytes <- if (!is.null(failure)) fail(failure) else tryCatch(
    .ctJuliaBuffered(wire, function(buffer) {
      writeBin(wire$RESULT, buffer)
      wire$writeElement(result)
    }), error = function(e) fail(conditionMessage(e)))
  .ctJuliaSend(wire, bytes, name)
}

# Read messages until the reply to `name` arrives. `discard` is the settling of
# an interrupted call: its output is dropped, a callback it makes is refused,
# and its result or failure is read and thrown away.
.ctJuliaAwait <- function(wire, name, discard = FALSE, announce = NULL) {
  repeat {
    type <- .ctJuliaNextMessage(wire, name, discard)
    done <- suspendInterrupts(.ctJuliaTakeMessage(wire, type, discard, name))
    if (is.list(done)) {
      if (identical(done$kind, "result")) return(done$value)
      if (discard) return(invisible(NULL))
      # Julia interrupted by something other than R -- a terminal's Ctrl-C
      # reaching it before `_ctsem_detach_console()` could, or a session
      # started without it. The reply is whole, so the session is fine; what
      # was asked for is an interrupt, not an engine error.
      if (grepl("InterruptException", done$message, fixed = TRUE)) .ctJuliaInterruptNow()
      stop(done$message, call. = FALSE, domain = NA)
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

# The body of a message whose first byte is `type`.
#
# Output and failure messages are read here rather than by JuliaConnectoR's
# readers, which loop until they have the bytes they asked for and so spin
# forever, at full CPU and with interrupts held, on a stream that ends partway
# through a message. A Julia that dies does that: the relay of its stderr is
# the likeliest thing to be mid-message when it goes, and on Linux the next
# call after its engine was killed read the first byte of one and hung. A
# result or a callback is still read by JuliaConnectoR, whose element format is
# its own; one cut off mid-message remains out of reach, and is far rarer.
#
# A pause partway through a message is not an ending, however long: the relay
# stops after a message's marker whenever Julia's main thread is busy, and a
# five-second limit on that once declared a session loading the engine dead
# (`.ctJuliaReadBytes()`). Only the end of the stream is. The wait holds
# interrupts, as JuliaConnectoR's readers always did, since Escape partway
# through a message would leave the stream out of step.
.ctJuliaTakeMessage <- function(wire, type, discard, name) {
  if (identical(type, wire$RESULT)) {
    value <- wire$readElement()
    return(list(kind = "result", value = if (discard) NULL else value))
  }
  if (identical(type, wire$FAIL)) {
    message <- .ctJuliaReadText(wire$pkgLocal$con, wait = Inf)
    if (is.null(message)) .ctJuliaLost(name)
    return(list(kind = "fail", message = message))
  }
  if (identical(type, wire$STDOUT) || identical(type, wire$STDERR)) {
    output <- .ctJuliaReadOutput(wire$pkgLocal$con, wait = Inf)
    if (is.null(output)) .ctJuliaLost(name)
    if (!discard) cat(output, file = if (identical(type, wire$STDOUT)) stdout() else stderr())
    return(NULL)
  }
  if (identical(type, wire$CALL)) {
    call <- wire$readCall()
    # Refused rather than answered when discarding: the progress display or
    # callback it is for belongs to a fit that is gone. The engine's callers
    # disable a callback that fails and carry on to their next checkpoint.
    .ctJuliaAnswer(wire, call, name,
      refuse = if (discard) "the call was interrupted from R")
    return(NULL)
  }
  .ctJuliaKill()
  stop("The Julia session sent a message ctsem could not read, so it was stopped. ",
    "The next call starts a new one.", call. = FALSE)
}

# A length-prefixed string as JuliaConnectoR frames one -- four bytes of length,
# then the bytes -- or NULL if the stream ends first, or goes silent for `wait`
# seconds (`.ctJuliaReadBytes()`).
.ctJuliaReadText <- function(connection, wait = 5) {
  len <- .ctJuliaReadBytes(connection, 4L, wait = wait)
  if (is.null(len)) return(NULL)
  n <- readBin(len, "integer", size = 4L)
  body <- if (n > 0L) .ctJuliaReadBytes(connection, n, wait = wait) else raw(0)
  if (is.null(body)) return(NULL)
  text <- tryCatch(rawToChar(body), error = function(e) "")
  Encoding(text) <- "UTF-8"
  text
}

# The text of one relayed stdout or stderr message, or NULL as for
# `.ctJuliaReadText()`. Escape sequences are stripped, as JuliaConnectoR's
# readOutput does.
.ctJuliaReadOutput <- function(connection, wait = 5) {
  text <- .ctJuliaReadText(connection, wait = wait)
  if (is.null(text)) return(NULL)
  gsub("\033(?:[@-Z\\\\-_]|\\[[0-?]*[ -/]*[@-~])", "", text)
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
  # A checkpointed loop answers within a fraction of a second, and then there
  # is nothing to say. Anything slower is work going on in the background --
  # compiling, typically, which can take minutes -- and that is worth one line.
  con <- tryCatch(get("pkgLocal", envir = asNamespace("JuliaConnectoR"))$con,
    error = function(e) NULL)
  answered <- !is.null(con) &&
    isTRUE(tryCatch(socketSelect(list(con), timeout = 0.3), error = function(e) TRUE))
  if (!answered) {
    message("Julia is finishing the interrupted step in the background; the next ",
      "julia call waits for it. ctJuliaKill(\"session\") stops it now.")
  }
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
  # Only an id read from the connection still open; see .ctJuliaEnsureStarted().
  pid <- if (identical(.ctJuliaSessionStamp(), .ct_julia_cache$pid_session)) {
    .ct_julia_cache$pid
  }
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
    # Idle, it ends when JuliaConnectoR closes the socket; either way the
    # record of it goes (ctJuliaProcesses()).
    if (!is.null(env$pid)) try(.ctJuliaUnregister(env$pid), silent = TRUE)
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
  # Keyed on the connection itself, as the objective cache is: a process id
  # read from a session that has since been replaced names a process that has
  # gone, and Windows hands ids out again, so killing by it could end something
  # else entirely.
  current <- .ctJuliaSessionStamp()
  if (started && !is.null(.ct_julia_cache$interrupt_file) &&
      identical(current, .ct_julia_cache$pid_session)) {
    return(invisible(FALSE))
  }
  external <- nzchar(Sys.getenv("JULIACONNECTOR_SERVER", unset = ""))
  if (!started) {
    if (!external) {
      .ctJuliaProvision()
      .ctJuliaAnnounce()
    }
    suspendInterrupts(withCallingHandlers(wire$ensure(), message = function(m) {
      if (!external && startsWith(conditionMessage(m), "Starting Julia")) {
        invokeRestart("muffleMessage")
      }
    }))
  }
  identity <- if (external) NULL else tryCatch(as.integer(strsplit(as.character(
    .ctJuliaExchange(wire, "RConnector.mainevalcmd",
      list("string(getpid(), \" \", Threads.nthreads())"))), " ")[[1L]]),
    error = function(e) NULL)
  .ct_julia_cache$pid <- if (length(identity) == 2L) identity[[1L]] else NULL
  .ct_julia_cache$pid_session <- .ctJuliaSessionStamp()
  if (!is.null(.ct_julia_cache$interrupt_file)) unlink(.ct_julia_cache$interrupt_file)
  .ct_julia_cache$interrupt_file <- gsub("\\", "/",
    tempfile("ctsem-julia-interrupt-"), fixed = TRUE)
  .ctJuliaRegisterExit()
  # Written down, so ctJuliaProcesses() can tell this engine from any other
  # Julia, and an R session that ends without cleaning up leaves a trace.
  if (length(identity) == 2L) .ctJuliaRegister(identity[[1L]], identity[[2L]])
  if (!external) .ctJuliaWarnOrphans()
  invisible(TRUE)
}

# How many threads a session starts with ---------------------------------------
#
# Julia fixes its thread count when the process starts; a running process cannot
# gain threads. A session started at two could therefore honour a later
# `cores = 8` only by restarting, which discards every model shape it has
# compiled. So a session is started with as many threads as this R process may
# use, and each call is held to its own `cores` by the engine's chunk ceiling,
# as it always was: the threads are there when asked for, idle when not.
#
# `parallelly::availableCores()` is "may use": every core on a desktop, the
# allocation under Slurm, PBS or a cgroup limit, `mc.cores` when that is set,
# and 2 under `R CMD check`. Measured on a 24-thread Windows desktop, a bare
# Julia at 24 threads against 2: the same start time, 30 MB more memory, and no
# CPU at all while idle.
.ctJuliaWidth <- function() {
  n <- tryCatch(suppressWarnings(as.integer(parallelly::availableCores())[1L]),
    error = function(e) NA_integer_)
  if (is.na(n) || n < 1L) 2L else n
}

# Set the thread count for a session about to start, unless someone chose it:
# ctJuliaSetup(threads=), the user's environment and a scheduler all win. At
# least `cores` wide when a call asks for more than the default. The value set
# is remembered, so a later start can tell its own setting from a deliberate one.
.ctJuliaProvision <- function(cores = NA_integer_) {
  existing <- Sys.getenv("JULIA_NUM_THREADS", unset = "")
  if (nzchar(existing) && !identical(existing, .ct_julia_cache$threads_from_cores)) {
    return(invisible(suppressWarnings(as.integer(sub(",.*$", "", existing)))))
  }
  width <- max(.ctJuliaWidth(), suppressWarnings(as.integer(cores)[1L]), na.rm = TRUE)
  Sys.setenv(JULIA_NUM_THREADS = as.character(width))
  .ct_julia_cache$threads_from_cores <- as.character(width)
  invisible(width)
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
# that is what ctJuliaInstall() installs. `pinned` and `latest` are arguments
# so the rule can be tested without either.
.ctJuliaNewerRelease <- function(version, bin, managed = NULL,
  pinned = .ct_julia_version, latest = .ctJuliaupLatest) {
  if (length(version) != 1L || is.na(version) || is.null(bin)) return(NULL)
  current <- tryCatch(numeric_version(version), error = function(e) NULL)
  if (is.null(current)) return(NULL)
  series <- function(v) sub("^(\\d+\\.\\d+).*$", "\\1", v)
  running <- series(version)
  tested <- series(pinned)
  if (is.null(managed)) {
    managed <- startsWith(normalizePath(bin, winslash = "/", mustWork = FALSE),
      .ctJuliaInstallRoot())
  }
  candidates <- pinned
  if (!managed) {
    candidates <- c(candidates, latest(running),
      if (numeric_version(tested) > numeric_version(running)) latest(tested))
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
    if (series(best) != running) paste0(" (ctsem is tested on ", tested, ")") else ""
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
