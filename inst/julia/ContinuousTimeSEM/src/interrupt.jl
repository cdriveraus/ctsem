# Stopping a call from R.
#
# Escape in R cannot reach a call that is running here. The bridge is one
# request and one reply on a socket, and this process reads the socket only
# between calls, so a request to stop sent that way would sit unread until the
# call it was meant to stop had finished. Signals are no better: Windows has no
# way to deliver one to a process that shares no console with R.
#
# So R leaves a file where this process looks. `ctsem_set_interrupt!` names it
# once per session; on Escape, R creates it and returns to its prompt without
# waiting for the reply, and the next checkpoint here sees it and throws an
# `InterruptException`. That unwinds the call to the bridge, which sends the
# failure as the reply, and R reads and discards that reply before its next
# request. Nothing is killed, so every method this session has compiled is kept.
#
# The checkpoints are the per-iteration hooks every long loop already calls --
# the optimiser's trace record, the progress and callback clocks -- so they cost
# one clock read, and a file lookup at most every 0.2 s. Work with no iteration
# in it, compilation above all, is not interruptible here and finishes first.
#
# The same check notices when R itself has gone. A process whose R has exited
# or crashed would otherwise run its call to the end, possibly for hours,
# because only an idle bridge notices the socket closing; now it exits at the
# next checkpoint instead.
#
# An `InterruptException` because the engine already treats one as the user
# stopping the run: the catches around trial evaluations rethrow it
# (`_ctsem_must_propagate`) rather than scoring the point as invalid.

"""File whose existence means R has asked the running call to stop; empty for none."""
const _CTSEM_INTERRUPT_FILE = Ref("")

"""The R process this session serves; 0 when unknown."""
const _CTSEM_PARENT_PID = Ref(0)

"""When the file may next be looked for. Atomic because chains run on threads."""
const _CTSEM_INTERRUPT_NEXT = Threads.Atomic{Float64}(0.0)

"""Whether the last look found a request to stop."""
const _CTSEM_INTERRUPT_SEEN = Threads.Atomic{Bool}(false)

"""
    ctsem_set_interrupt!(path, parent)

Name the file R creates to stop a running call, and the R process whose exit
should end this one. Returns this process's id, which R needs to stop the
process outright when a call cannot be interrupted.
"""
function ctsem_set_interrupt!(path::AbstractString, parent::Integer)
    _CTSEM_INTERRUPT_FILE[] = String(path)
    _CTSEM_PARENT_PID[] = Int(parent)
    _CTSEM_INTERRUPT_NEXT[] = 0.0
    _CTSEM_INTERRUPT_SEEN[] = false
    return Int(getpid())
end

"""
    _ctsem_parent_alive()

Whether the R process is still running. Answers `true` whenever that cannot be
established for certain: exiting on a wrong `false` would take a live R
session's engine away mid-call, while a wrong `true` only leaves things as they
were before this check existed.
"""
function _ctsem_parent_alive()
    pid = _CTSEM_PARENT_PID[]
    pid <= 0 && return true
    if Sys.iswindows()
        # PROCESS_QUERY_LIMITED_INFORMATION. A null handle with
        # ERROR_INVALID_PARAMETER is the one answer that means no such process.
        handle = ccall((:OpenProcess, "kernel32"), stdcall, Ptr{Cvoid},
            (UInt32, Int32, UInt32), 0x1000, 0, UInt32(pid))
        if handle == C_NULL
            return Libc.GetLastError() != 87
        end
        code = Ref{UInt32}(0)
        ok = ccall((:GetExitCodeProcess, "kernel32"), stdcall, Int32,
            (Ptr{Cvoid}, Ref{UInt32}), handle, code)
        ccall((:CloseHandle, "kernel32"), stdcall, Int32, (Ptr{Cvoid},), handle)
        # STILL_ACTIVE
        return ok == 0 || code[] == 259
    end
    ccall(:kill, Cint, (Cint, Cint), pid, 0) == 0 && return true
    return Libc.errno() != Libc.ESRCH
end

"""
    _ctsem_interrupted()

Whether R has asked the running call to stop. A negative answer is reused for
0.2 s; a positive one is re-checked every time, so the call R sends after
removing the file is not stopped by what the previous one saw.
"""
function _ctsem_interrupted()
    path = _CTSEM_INTERRUPT_FILE[]
    isempty(path) && return false
    now = time()
    if !_CTSEM_INTERRUPT_SEEN[] && now < _CTSEM_INTERRUPT_NEXT[]
        return false
    end
    _CTSEM_INTERRUPT_NEXT[] = now + 0.2
    if !_ctsem_parent_alive()
        exit(1)
    end
    seen = isfile(path)
    _CTSEM_INTERRUPT_SEEN[] = seen
    return seen
end

"""Throw an `InterruptException` if R has asked the running call to stop."""
function _ctsem_interrupt_check()
    _ctsem_interrupted() && throw(InterruptException())
    return nothing
end
