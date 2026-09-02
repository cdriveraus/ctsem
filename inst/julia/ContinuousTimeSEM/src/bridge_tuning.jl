# Socket options for the R bridge.
#
# JuliaConnectoR talks to R over a loopback TCP socket, and both ends assemble a
# message from many small `write`s: an indicator byte, a length, a name, a count,
# then each argument. That is the textbook write-write-read pattern, and on Linux
# it meets Nagle's algorithm on the sender and a delayed ACK on the receiver: the
# first small segment goes out, the second waits for an ACK that the kernel holds
# back for up to 40 ms. Measured on a bare socket with no Julia and no R package
# involved (23-core Linux box): one write per request 0.033 ms, two writes 40.7
# ms, eight writes 40.8 ms -- one stall per message, whatever the payload.
#
# That is where the bridge's cost was. A message cost 41 ms if one direction
# split, 82 ms if both did; `ctJuliaEvaluate` sends six and paid ~370 ms, which
# had been read as nine round trips of 42 ms each. It is not nine of anything.
#
# Two options remove it, and they are not interchangeable:
#
#   nagle(io, false)  disables Nagle on *this* socket, so Julia's replies leave
#                     immediately. Permanent, and works on every OS.
#
#   quickack(io, true) tells Linux to acknowledge R's segments at once instead of
#                     delaying, which is what releases R's Nagle. R sets no
#                     socket options and offers no way to, so this is the only
#                     side the R->Julia direction can be fixed from. TCP_QUICKACK
#                     is explicitly not sticky -- the kernel drops back to
#                     delayed ACKs after a handful of packets -- so it has to be
#                     re-armed, and it has to be re-armed *while R is mid
#                     message*, which no call from R can do. A Julia task can:
#                     it runs precisely when the main task is parked in the
#                     socket read, which is that window.
#
# Same box, per bridge message: 82 ms default, 41 ms with Nagle off, 2.5 ms with
# both. Non-Linux is unaffected either way -- Windows measures ~0.3 ms a message
# untuned, and `quickack` compiles to nothing off Linux.

using Sockets

# The re-arming task, so a second call does not start a second one.
const _CTSEM_QUICKACK_TASK = Ref{Any}(nothing)

"""
    ctsem_tune_bridge!(communicator; quickack_interval = 0.001)

Set socket options on the JuliaConnectoR connection held by `communicator`, and
on Linux start a task that keeps TCP_QUICKACK armed.

Returns `1` if the quickack task is running, `0` if only Nagle was disabled.
Idempotent: calling it again does not start a second task.

`communicator` is JuliaConnectoR's own `RConnector.CommunicatoR`, reached
through its `io` field by name rather than by type, so this file needs no
dependency on JuliaConnectoR. A `quickack_interval` of zero, or any non-Linux
platform, skips the task and leaves the socket merely Nagle-free.
"""
function ctsem_tune_bridge!(communicator; quickack_interval::Real = 0.001)
    io = getfield(communicator, :io)
    Sockets.nagle(io, false)
    if !(Sys.islinux() && quickack_interval > 0)
        return Int32(0)
    end
    Sockets.quickack(io, true)
    task = _CTSEM_QUICKACK_TASK[]
    if task !== nothing && !istaskdone(task)
        return Int32(1)
    end
    interval = float(quickack_interval)
    _CTSEM_QUICKACK_TASK[] = @async begin
        try
            while isopen(io)
                Sockets.quickack(io, true)
                sleep(interval)
            end
        catch
            # A closed or replaced socket ends the task; it must never take the
            # session down with it.
        finally
            _CTSEM_QUICKACK_TASK[] = nothing
        end
    end
    return Int32(1)
end

"""
    ctsem_bridge_tuned()

Whether the TCP_QUICKACK re-arming task started by [`ctsem_tune_bridge!`](@ref)
is currently running. For tests and diagnostics.
"""
function ctsem_bridge_tuned()
    task = _CTSEM_QUICKACK_TASK[]
    return task !== nothing && !istaskdone(task)
end
