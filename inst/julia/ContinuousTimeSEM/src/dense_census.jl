# TEMPORARY (branch sparseunits): records where dense unit-sized linear algebra
# runs, by call chain. Removed before merge.
const _DENSE_CENSUS = Dict{String,Vector{Int}}()
const _DENSE_CENSUS_LOCK = ReentrantLock()
const _DENSE_CENSUS_MIN = Ref(typemax(Int))
@noinline function _dense_census(kind::String, n::Integer)
    n >= _DENSE_CENSUS_MIN[] || return nothing
    frames = String[]
    for f in stacktrace()
        s = string(f.func)
        (s == "_dense_census" || s == kind || startswith(s, "_ctsem_jacobi") ||
            startswith(s, "_ctsem_symeig") || startswith(s, "_ctsem_cholesky") ||
            startswith(s, "_laplace_block_dense") || startswith(s, "#")) && continue
        push!(frames, string(s, ":", f.line))
        length(frames) >= 5 && break
    end
    key = string(kind, " <- ", join(frames, " <- "))
    lock(_DENSE_CENSUS_LOCK) do
        v = get!(_DENSE_CENSUS, key, [0, 0])
        v[1] += 1
        v[2] = max(v[2], Int(n))
    end
    return nothing
end
function ctsem_dense_census(min::Integer)
    _DENSE_CENSUS_MIN[] = min
    lock(_DENSE_CENSUS_LOCK) do
        out = [string(v[1], "x maxn=", v[2], "  ", k) for (k, v) in _DENSE_CENSUS]
        empty!(_DENSE_CENSUS)
        sort(out)
    end
end
export ctsem_dense_census
@noinline function _census_log_error(e, bt)
    path = get(ENV, "CTSEM_ERRLOG", "")
    isempty(path) && return nothing
    try
        r = e
        b = bt
        while r isa TaskFailedException
            st = Base.current_exceptions(r.task)
            r, b = st[end].exception, st[end].backtrace
        end
        s = sprint(showerror, r, b; context=:stacktrace_types_limited => Ref(true))
        open(path, "a") do io
            println(io, "==== thread ", Threads.threadid(), " ", typeof(r))
            println(io, first(s, 20000))
        end
    catch
    end
    return nothing
end
