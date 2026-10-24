# Evaluating an objective at many points at once.
#
# Tempered SMC moves its particles and DAIS evaluates its draws, each point
# independent of the others, and both spend nearly all their time there. One
# layout serves both (`_each_point!`):
#
#   * On the Laplace route the points go to the pool (`_laplace_parallel`),
#     each worker task holding one copy of the objective (`_laplace_draw_copy`:
#     its own modes and flags, everything else shared), and the pool's nesting
#     bounds the workers alive by the `cores` ceiling whether a point's
#     evaluation asks for helpers over its units or not.
#   * On a filter objective whose subjects are fewer than the pool -- a time
#     series is one -- the subject loop cannot fill it, so the points are
#     spread over deep copies of the objective, one a worker, each task
#     running its share with the chunk ceiling held at one so that the copies
#     together stay within the ceiling.
#   * Otherwise, one after another on the shared objective, each evaluation's
#     subject loop dividing the pool as it always does.
#
# Each task's copy and scratch are built inside a function of its own
# (`_points_share`): built inline in a spawned block, names the enclosing
# function also used were rebound, and every task wrote one shared NUTS
# workspace (CLAUDE.md, "Julia threads").

"""
    _points_mode(objective, parallel)

How a set of points is evaluated: `:serial`, one after another on the shared
objective; `:pool`, over the Laplace pool on copies; `:copies`, over deep
copies of a filter objective whose subjects are fewer than the pool. The
pool's width is the `cores` ceiling capped by the threads there are.
"""
function _points_mode(objective, parallel::Bool)
    (parallel && Threads.nthreads() > 1 && _laplace_pool_width() > 1) || return :serial
    objective isa CTSEMLaplaceObjective && return :pool
    if objective isa CTSEMObjective &&
        length(objective.subject_objectives) < _laplace_pool_width()
        return :copies
    end
    return :serial
end

"""A copy of `laplace` for one worker's points: fresh modes for the units it
will solve (the coreset's, or all), its own flags, the rest shared
(`_laplace_draw_copy`)."""
function _points_laplace_copy(laplace::CTSEMLaplaceObjective, coreset)
    cs = coreset === nothing ? CTSEMCoreset(collect(1:length(laplace.units.members)),
        ones(length(laplace.units.members))) : coreset
    return _laplace_draw_copy(laplace, cs)
end

"""No scratch: a point that needs none beside its objective."""
_points_no_scratch() = nothing

"""
    _each_point!(f, objective, n, mode; coreset=nothing, scratch=_points_no_scratch)

`f(i, objective_i, scratch_i)` for every point `i` in `1:n`, with the
objective it may evaluate and scratch of its own (`scratch()` builds one per
task). `mode` is `_points_mode`'s. On the pool a worker's copy solves the
coreset's units, or all where `coreset` is `nothing`. `f` writes only into its
point's own columns and slots.
"""
function _each_point!(f, objective, n::Int, mode::Symbol; coreset=nothing,
    scratch=_points_no_scratch)
    if mode == :pool && objective isa CTSEMLaplaceObjective
        try
            _laplace_parallel(i -> _points_pool_item(f, objective, coreset, scratch, i),
                objective, 1:n)
        finally
            # The calling task runs points too, and its storage outlives the
            # region: drop its copy rather than hold it for the session.
            delete!(task_local_storage(), _POINTS_KEY)
        end
        return nothing
    end
    if mode == :copies
        ncopies = min(_laplace_pool_width(), n)
        ceiling = _CTSEM_MAX_CHUNKS[]
        _CTSEM_MAX_CHUNKS[] = 1
        try
            tasks = [Threads.@spawn(_points_share(f, objective, c, ncopies, n, scratch))
                for c in 1:ncopies]
            foreach(wait, tasks)
        finally
            _CTSEM_MAX_CHUNKS[] = ceiling
        end
        return nothing
    end
    _points_share(f, objective, 1, 1, n, scratch; copy=false)
    return nothing
end

const _POINTS_KEY = :ctsem_point_scratch

"""
    _points_pool_item(f, objective, coreset, scratch, i)

Point `i` on the pool: on this worker task's copy of the objective and its
scratch, built the first time the task runs a point and kept in its storage
for the rest of the region. A function of its own, not the region's closure:
written inline, its locals were the enclosing function's, and the tasks
shared one (CLAUDE.md, "Julia threads").
"""
function _points_pool_item(f, objective, coreset, scratch, i::Int)
    _ctsem_interrupt_check()
    st = task_local_storage()
    mine = get(st, _POINTS_KEY, nothing)
    if mine === nothing
        mine = (obj=_points_laplace_copy(objective, coreset), scratch=scratch())
        st[_POINTS_KEY] = mine
    end
    f(i, mine.obj, mine.scratch)
    return true
end

"""
    _points_share(f, objective, first, stride, n, scratch; copy)

`f(i, objective_i, scratch_i)` for points `first:stride:n` in turn, on a deep
copy of `objective` (or the objective itself with `copy = false`) and scratch
that belongs to this call alone.
"""
function _points_share(f, objective, first::Int, stride::Int, n::Int, scratch;
    copy::Bool=true)
    obj = copy ? deepcopy(objective) : objective
    own = scratch()
    for i in first:stride:n
        _ctsem_interrupt_check()
        f(i, obj, own)
    end
    return nothing
end
