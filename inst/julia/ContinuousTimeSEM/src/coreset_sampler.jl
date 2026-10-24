# The importance sampler's coreset -------------------------------------------
#
# A Laplace fit's log density is a sum of per-unit terms, each an inner mode
# solve, and an importance sampler pays for every unit at every draw. Over the
# region the sampler visits the units' terms vary in about as many dimensions
# as it moves, so M of the N units with non-negative weights and a free
# intercept, fitted by least squares to the full sum, follow it, at M/N of a
# full evaluation a draw. `CTSEMISCoreset` is that density: R's
# `.ctCoresetLpg()` wraps it for `imis_is`, and every batch the sampler asks
# for is one call here, with the training rows (every unit's term at the draws
# evaluated in full, and at a few of each full gradient batch's draws every
# unit's score) kept where they are computed.
#
# Every decision is by the expected cost of the sampler's work in full
# evaluations, unit terms counted, never in seconds, so the same fit makes the
# same choices anywhere. It evaluates in full, as the sampler would without it,
# and learns from that until some coreset is expected to make a round cheaper;
# then gradients run on that coreset's units only, and each batch of draws does
# too where that is cheaper for its size: a few of its draws are checked in
# full, and the rest is evaluated in full when the check fails. A batch too
# small to pay for its check is evaluated in full, and the coreset's prediction
# at its draws, read off the unit terms, checks it for nothing. The size is
# the one whose expected round cost is least (`_is_coreset_choose`): its share
# of the round's values and gradients, its checks, and the chance a check
# fails times the full batch that costs, read from held-out errors and then
# from the checks the coreset actually fails. At the end it checks the draws
# that carry the weight and corrects them where needed, so what it hands back
# is the full density's. The selection itself is `_coreset_select`'s
# (coreset_select.jl), which the optimiser's stages and SMC share.
#
# Draws are columns, as everywhere in the engine.

export ctsem_is_coreset, ctsem_is_coreset_batch!, ctsem_is_coreset_gradbatch!,
    ctsem_is_coreset_fullgrad!, ctsem_is_coreset_finalize!, ctsem_is_coreset_record,
    ctsem_is_coreset_cost

"""Full evaluations the coreset density has spent, a gradient as two."""
ctsem_is_coreset_cost(cs) = cs.cost

"""
The tolerances, in log density units. A weight is off by exp(error), so an RMS
of 0.1 costs about 1% of an effective sample; `max` keeps any one draw within a
factor of 1.6. A round's checks may drift further before the coreset is
rebuilt (`drift`), so that it is not rebuilt on noise. Both are read where a
draw can carry weight: a sampler's first rounds are mostly tails, 50-150 log
units below its best draw, where a coreset within 0.1 near the top errs by a
unit or two, and a draw 50 down carries exp(-50) of the best one's weight. So a
draw counts fully within `band` of the best and e-fold less every `scale`
below that (`_is_relevance`), in the fit as in the test, and `max` reads only
draws of relevance at least `relevant`, about 55 below the best. `final` is the
end check's. `force` is how far a coreset's gradient may be from the full
sum's at held-out scored draws, in log density per sd each parameter moves
across the draws (`_is_force_error`), for its size to be priced at all: a
gradient off by e moves the values' error by about e for each sd the draws
move, so at a log density unit per sd a coreset cannot follow them. Between
that and a check's tolerance the value checks decide: held at 0.3, twice a
check's drift, the gate priced coresets of 84 units on the 250-unit strong
shell where 47-49 had served, at 0.82 full evaluations a gradient against
0.58-0.62. The sizes that failed there had gradients off by 0.55-11, those
that served 0.14-0.29 (review/SCORES-coreset-training-2026-10-24.md).
"""
const _IS_CORESET_TOL = (rms=0.1, max=0.5, drift=1.5, final=0.5, band=10.0,
    scale=20.0, relevant=0.1, force=1.0)

"""Scores: from the `per` draws of each full gradient batch of highest total,
where the weight sits, and from the newest `kept` such draws in all. A draw's
scores are a row per parameter over every unit, so they are what the training
matrix grows by; 30 draws were enough on the strong shell at 250 and 1,000
units (review/SCORES-coreset-training-2026-10-24.md). Drawn spread through the
batch instead, they were rows the fit weighs by relevance: DAIS's batches on
the strong shell at 1,000 units put their weight on one or two draws, and
spread scores came to almost nothing there."""
const _IS_SCORE_DRAWS = (per=4, kept=64)

"""Each draw's weight in a fit or a test: 1 within `band` of `top`, e-fold less
every `scale` below."""
_is_relevance(y, top) = [exp(-max(0.0, top - v - _IS_CORESET_TOL.band) /
    _IS_CORESET_TOL.scale) for v in y]

"""Whether errors `err` at draws of relevance `rel` pass the tolerances,
`mult` times over."""
function _is_passes(err::AbstractVector, rel::AbstractVector, mult::Real=1.0)
    (isempty(err) || sum(rel) <= 0) && return true
    rms = sqrt(sum(rel .* err .^ 2) / sum(rel))
    worst = 0.0
    for i in eachindex(err)
        rel[i] >= _IS_CORESET_TOL.relevant && (worst = max(worst, abs(err[i])))
    end
    return rms <= mult * _IS_CORESET_TOL.rms && worst <= mult * _IS_CORESET_TOL.max
end

"""What a sampler round is: its batch of draws, the gradients asked for before
the next (each as two values), and the checks a batch on a coreset takes."""
struct _ISRound
    values::Float64
    grads::Float64
    checks::Int
end

"""What a batch of `rd.values` draws is expected to cost on a coreset of `M`
of `N` units, and in full, in full evaluations: on the coreset its share
`M / N`, the checks evaluated in full, and with probability `fail` the rest in
full after a failed check. Measured on dev3, a coreset's draw took M/N of a
full one's time to within 5% at 100, 250 and 1,000 units."""
function _is_batch_cost(M::Real, N::Real, rd::_ISRound, fail::Real)
    n = rd.values
    k = min(rd.checks, n)
    return (coreset=M / N * n + k + fail * (n - k), full=n)
end

"""A sampler round's expected cost with a coreset of `M` of `N` units: the
batch whichever way is cheaper, and the gradients at its share where it passed
its last check and in full where it failed it, since a coreset serves
gradients only where it has just passed one."""
function _is_round_cost(M::Real, N::Real, rd::_ISRound, fail::Real)
    share = M / N
    bc = _is_batch_cost(M, N, rd, fail)
    return 2 * rd.grads * (share + fail * (1 - share)) + min(bc.coreset, bc.full)
end

"""A round in full."""
_is_full_round(rd::_ISRound) = rd.values + 2 * rd.grads

"""
How often a round's check would fail, from a coreset's errors `err` and
predictions `pred` at held-out draws of relevance `rel`. The check is simulated
as the sampler makes it: a round of `rd.values` draws resampled from these,
its checks the half the coreset ranks highest and an even spread of the rest.
Ranking by the coreset's own value is what finds a coreset too small: it
overrates the draws it puts on top, and a check drawn by relevance alone
missed that, chose sizes that failed in the sampler, and cost a full round each
time. On a private random stream, so that the sampler's own is the same with a
coreset as without.
"""
function _is_fail_rate(err::AbstractVector, rel::AbstractVector, pred::AbstractVector,
        rd::_ISRound)
    n = length(err)
    (n == 0 || sum(rel) <= 0) && return 0.0
    rng = Random.Xoshiro(20261021)
    k = max(1, trunc(Int, min(rd.checks, rd.values)))
    ntop = cld(k, 2)
    draws = max(k, trunc(Int, rd.values))
    failures = 0
    for _ in 1:200
        i = rand(rng, 1:n, draws)
        o = i[sortperm(-pred[i])]
        pick = collect(o[1:ntop])
        others = o[(ntop + 1):end]
        if k - ntop > 0 && !isempty(others)
            append!(pick, others[unique(round.(Int, range(1, length(others);
                length=k - ntop)))])
        end
        _is_passes(err[pick], rel[pick], _IS_CORESET_TOL.drift) || (failures += 1)
    end
    return failures / 200
end

"""The draws a check reads: the half of at most eight the coreset ranks
highest by `v`, where the weight will sit, and an even spread of the rest."""
function _is_checkpick(v::AbstractVector)
    k = min(8, length(v))
    top = sortperm(-v)[1:cld(k, 2)]
    rest = setdiff(1:length(v), top)
    m = fld(k, 2)
    pick = collect(top)
    if m > 0 && !isempty(rest)
        append!(pick, rest[unique(round.(Int, range(1, length(rest); length=m)))])
    end
    return pick
end

"""How many draws relevances `rel` amount to as the check weighs them, normalised
among themselves: Kish's effective size, `(sum rel)^2 / sum rel^2`."""
_is_effective(rel) = (s = sum(rel); s > 0 ? s^2 / sum(abs2, rel) : 0.0)
_is_key(x::AbstractVector) = [round(v; sigdigits=12) for v in x]

"""
The rows a coreset is fitted to: every unit's term at the draws evaluated in
full (`L`, units by draws, the first `nrows` in use), their totals (`y`), the
rows read as the newest checks (`focus`), and the score rows -- per row the
units' scores (`Z`, rows by units), their sum, its parameter, its draw's total
and a draw counter -- with the draws that gave them. Relevance is measured
from the best total of every row, or, where `batchtop`, of the newest batch's,
rows `batch_from` on (`_is_top`).
"""
mutable struct _ISRows
    N::Int
    L::Matrix{Float64}
    nrows::Int
    y::Vector{Float64}
    focus::Vector{Int}
    Z::Matrix{Float64}
    zsum::Vector{Float64}
    zpar::Vector{Int}
    ztotal::Vector{Float64}
    zdraw::Vector{Int}
    zdraws::Vector{Vector{Float64}}
    scoredraws::Int
    batchtop::Bool
    batch_from::Int
end

_ISRows(N::Integer; batchtop::Bool=false) = _ISRows(N, zeros(N, 0), 0, Float64[], Int[],
    zeros(0, N), Float64[], Int[], Float64[], Int[], Vector{Float64}[], 0, batchtop, 1)

"""
The total a draw's relevance is measured from (`_is_relevance`): the best of
every row, for a sampler whose weights compare every round's draws, as the
importance sampler's do; or, where the rows say `batchtop`, the best of the
newest batch's, for one whose weights are normalised within each batch, as
DAIS's are. Measured from the best of the run, DAIS's draws on the strong
shell at 1,000 units sat so far below the run's best draw that a held-out
set of hundreds amounted to two draws: its coresets were fitted and priced on
those two, chose a few units, and failed.
"""
_is_top(tr::_ISRows) = tr.batchtop && tr.batch_from <= tr.nrows ?
    maximum(view(tr.y, tr.batch_from:tr.nrows)) : maximum(tr.y)

"""A batch begins: the rows learned from here on are the newest batch's."""
_is_batch_start!(tr::_ISRows) = (tr.batch_from = tr.nrows + 1; nothing)

"""
    CTSEMISCoreset

The importance sampler's density with a coreset behind it: the training rows
(`_ISRows`), the coreset in force and the record of what it did
(`ctsem_is_coreset`).
"""
mutable struct CTSEMISCoreset
    laplace::CTSEMLaplaceObjective
    N::Int
    rows::_ISRows
    # The coreset in force: units, weights, intercept, and the expected round
    # cost and held-out failure rate it was chosen at.
    has_cf::Bool
    S::Vector{Int}
    w::Vector{Float64}
    c::Float64
    cf_cost::Float64
    cf_fail::Float64
    phase::Symbol
    cost::Float64
    sizes::Vector{Int}
    checks::Vector{Float64}
    final_checked::Int
    final_worst::Float64
    final_corrected::Int
    has_final::Bool
    # A round: a batch's size and the gradients asked for before the next,
    # each smoothed over the batches so far.
    values::Float64
    grads::Float64
    gradsnow::Float64
    batches::Int
    runs::Int
    fails::Int
    switches::Int
    runs_total::Int
    fails_total::Int
    free_total::Int
    trusted::Bool
    last_passed::Bool
    nchoose::Int
    proved::Int
    rows_at_choose::Int
    secs::Dict{Symbol,Float64}
    draws::Dict{Symbol,Float64}
    why::Dict{Symbol,Float64}
    known::Dict{Vector{Float64},Float64}
    trace::Vector{NTuple{8,Float64}}
    # The last full gradient batch's draws, values, totals and unit terms, which
    # a value batch at the same draws reads; the draws so served.
    lastgrad::Dict{Vector{Float64},Int}
    lastgrad_value::Vector{Float64}
    lastgrad_total::Vector{Float64}
    lastgrad_L::Matrix{Float64}
    reused::Int
end

"""
    ctsem_is_coreset(laplace; top=:run)

A coreset density over the Laplace objective `laplace`, for an importance
sampler: `ctsem_is_coreset_batch!` and `ctsem_is_coreset_gradbatch!` are its
value and gradient batches, `ctsem_is_coreset_fullgrad!` the full density's
gradient, `ctsem_is_coreset_finalize!` the end check, and
`ctsem_is_coreset_record` what it did. `top` says how the caller's weights are
normalised, which is what a draw's relevance is measured against
(`_is_top`): over the whole run (`:run`, the importance sampler) or within
each batch (`:batch`, DAIS).
"""
function ctsem_is_coreset(laplace::CTSEMLaplaceObjective; top::Symbol=:run)
    top in (:run, :batch) || throw(ArgumentError("top must be :run or :batch, not :$top"))
    N = length(laplace.units.members)
    secs = Dict(k => 0.0 for k in (:full, :coreset, :gradient_full, :gradient_coreset, :choose))
    draws = Dict(k => 0.0 for k in (:full, :coreset, :gradient_full, :gradient_coreset))
    why = Dict(k => 0.0 for k in (:training, :small, :check, :gradient_check, :failed_rest, :end))
    return CTSEMISCoreset(laplace, N, _ISRows(N; batchtop=top == :batch),
        false, Int[], Float64[], 0.0, NaN, NaN, :training, 0.0, Int[], Float64[],
        0, NaN, 0, false, NaN, NaN, 0.0, 0, 0, 0, 0, 0, 0, 0, false, true, 0, 0, 0,
        secs, draws, why, Dict{Vector{Float64},Float64}(), NTuple{8,Float64}[],
        Dict{Vector{Float64},Int}(), Float64[], Float64[], zeros(N, 0), 0)
end

"""The coresets chosen, of the choices tried (`cs.nchoose`, a row each of the
trace; code 1 is a choice made)."""
_is_chosen(cs::CTSEMISCoreset) = count(t -> t[3] == 1.0, cs.trace)

_is_smooth(old, new) = isnan(old) ? Float64(new) : 0.7 * old + 0.3 * new
_is_round(cs::CTSEMISCoreset) = _ISRound(cs.values,
    isnan(cs.grads) ? 0.3 * cs.values : cs.grads, 8)
_is_active_n(cs::CTSEMISCoreset) = count(>(0), cs.w)
_is_failnow(cs::CTSEMISCoreset) = (cs.fails + 2 * cs.cf_fail) / (cs.runs + 2)
_is_predict(cs::CTSEMISCoreset, L::AbstractMatrix) =
    [cs.c + sum(cs.w[i] * L[cs.S[i], j] for i in eachindex(cs.S); init=0.0)
     for j in 1:size(L, 2)]

function _is_tally!(cs::CTSEMISCoreset, kind::Symbol, secs::Real, n::Real=0)
    cs.secs[kind] += secs
    haskey(cs.draws, kind) && (cs.draws[kind] += n)
    return nothing
end

"""Every unit's term at the columns of `X`, through `ctsem_coreset_terms`."""
function _is_terms(cs::CTSEMISCoreset, X::AbstractMatrix)
    r = ctsem_coreset_terms(cs.laplace, X)
    return (value=r.value, total=r.total, L=r.unit_loglik)
end

"""The training rows of a batch's draws `idx`: those with every term finite."""
function _is_learn!(tr::_ISRows, t0, idx=1:length(t0.total); focus::Bool=false)
    ok = [j for j in idx if isfinite(t0.total[j]) && all(isfinite, view(t0.L, :, j))]
    isempty(ok) && return nothing
    if size(tr.L, 2) < tr.nrows + length(ok)
        grown = zeros(tr.N, max(2 * size(tr.L, 2), tr.nrows + length(ok), 64))
        copyto!(view(grown, :, 1:tr.nrows), view(tr.L, :, 1:tr.nrows))
        tr.L = grown
    end
    focus && append!(tr.focus, tr.nrows .+ (1:length(ok)))
    for j in ok
        tr.nrows += 1
        tr.L[:, tr.nrows] .= view(t0.L, :, j)
        push!(tr.y, t0.total[j])
    end
    return nothing
end

"""The full density's terms at `X`: those of the last gradient batch read
back, the rest evaluated, all counted against `why`."""
function _is_full!(cs::CTSEMISCoreset, X::AbstractMatrix, why::Symbol)
    n = size(X, 2)
    keys = [_is_key(view(X, :, j)) for j in 1:n]
    hit = [get(cs.lastgrad, k, 0) for k in keys]
    todo = findall(==(0), hit)
    if length(todo) == n
        t = time()
        t0 = _is_terms(cs, X)
        _is_tally!(cs, :full, time() - t, n)
    else
        value = Vector{Float64}(undef, n)
        total = Vector{Float64}(undef, n)
        L = Matrix{Float64}(undef, cs.N, n)
        for j in 1:n
            hit[j] == 0 && continue
            value[j] = cs.lastgrad_value[hit[j]]
            total[j] = cs.lastgrad_total[hit[j]]
            L[:, j] .= view(cs.lastgrad_L, :, hit[j])
        end
        if !isempty(todo)
            t = time()
            tt = _is_terms(cs, X[:, todo])
            _is_tally!(cs, :full, time() - t, length(todo))
            value[todo] .= tt.value
            total[todo] .= tt.total
            L[:, todo] .= tt.L
        end
        cs.reused += n - length(todo)
        t0 = (value=value, total=total, L=L)
    end
    cs.cost += length(todo)
    cs.why[why] += length(todo)
    for j in 1:n
        isfinite(t0.value[j]) && (cs.known[keys[j]] = t0.value[j])
    end
    return t0
end

"""The coreset's value at the columns of `X`, solving its units only."""
function _is_coreset_value!(cs::CTSEMISCoreset, X::AbstractMatrix)
    units = cs.S[cs.w .> 0]
    weights = cs.w[cs.w .> 0]
    t = time()
    v = ctsem_coreset_batch(cs.laplace, X, units, weights)
    _is_tally!(cs, :coreset, time() - t, size(X, 2))
    cs.cost += size(X, 2) * length(units) / cs.N
    return v .+ cs.c
end

"""
The score rows for a fit: each stored row (one parameter of one draw: the
units' scores, and their sum) times the parameter's spread over the draws that
gave scores -- a force in posterior sds, so no parameter's units weigh by its
scale alone -- and by the root of its draw's relevance against the value rows'
best (`top`). `nothing` with none.
"""
function _is_score_rows(tr::_ISRows, top::Real)
    isempty(tr.zsum) && return nothing
    p = length(tr.zdraws[1])
    spread = ones(p)
    if length(tr.zdraws) > 1
        for i in 1:p
            # One the draws do not move weighs nothing: its gradient is not
            # what a check reads (the importance sampler holds all but the
            # check's named directions at its centre).
            s = std([d[i] for d in tr.zdraws])
            spread[i] = isfinite(s) ? s : 1.0
        end
    end
    f = spread[tr.zpar] .* sqrt.(_is_relevance(tr.ztotal, top))
    return (Z=tr.Z .* f, z=tr.zsum .* f)
end

"""
The least-squares rows for weights on the training rows `rows` with
relevances `rel`: each centred by its relevance-weighted mean and scaled by
the root of its relevance (the free intercept), and the score rows `score`
beside them, which the intercept does not reach. Returns `(A, b, mx, my)`, the
rows by units and the means the intercept is read from.
"""
function _is_rows(tr::_ISRows, rows::AbstractVector{Int}, rel::AbstractVector,
        score)
    share = rel ./ sum(rel)
    root = sqrt.(rel)
    nr = length(rows)
    ns = score === nothing ? 0 : length(score.z)
    A = Matrix{Float64}(undef, nr + ns, tr.N)
    mx = zeros(tr.N)
    @inbounds for U in 1:tr.N
        m = 0.0
        for (i, row) in enumerate(rows)
            m += share[i] * tr.L[U, row]
        end
        mx[U] = m
        for (i, row) in enumerate(rows)
            A[i, U] = (tr.L[U, row] - m) * root[i]
        end
        for i in 1:ns
            A[nr + i, U] = score.Z[i, U]
        end
    end
    my = sum(share[i] * tr.y[row] for (i, row) in enumerate(rows))
    b = vcat([(tr.y[row] - my) * root[i] for (i, row) in enumerate(rows)],
        ns == 0 ? Float64[] : score.z)
    return A, b, mx, my
end

"""Weights on units `S` over every training row (`passed`'s refit)."""
function _is_weights(tr::_ISRows, S::Vector{Int})
    rows = collect(1:tr.nrows)
    A, b, mx, my = _is_rows(tr, rows, _is_relevance(tr.y, _is_top(tr)),
        _is_score_rows(tr, _is_top(tr)))
    As = A[:, S]
    w = _ctsem_nnls(transpose(As) * As, transpose(As) * b)
    return w, my - sum(mx[S] .* w)
end

"""The scored draws held out of a fit: every fourth, in the order scored."""
function _is_held_draws(zdraw::AbstractVector{Int})
    ids = unique(zdraw)
    return Set(ids[[i % 4 == 0 for i in eachindex(ids)]])
end

"""
How far a coreset (`units`, weights `w`) is from the full sum in gradient at
the held-out scored draws `zheld`: per draw the norm of the gradient's error,
each parameter's in log density per sd it moves across the draws (the score
rows `score`, which also carry the root of the draw's relevance below `top`),
as a root mean square over the draws by their relevance. 0 where none is held
out, which lets the values decide.
"""
function _is_force_error(score, tr::_ISRows, zheld, units::AbstractVector{Int},
        w::AbstractVector, top::Real)
    (score === nothing || isempty(zheld)) && return 0.0
    err2 = Dict{Int,Float64}()
    for r in eachindex(tr.zdraw)
        tr.zdraw[r] in zheld || continue
        e = -score.z[r]
        for k in eachindex(units)
            e += w[k] * score.Z[r, units[k]]
        end
        err2[tr.zdraw[r]] = get(err2, tr.zdraw[r], 0.0) + e^2
    end
    rel = Dict(tr.zdraw[r] => _is_relevance([tr.ztotal[r]], top)[1]
        for r in eachindex(tr.zdraw) if tr.zdraw[r] in zheld)
    total = sum(values(rel); init=0.0)
    total > 0 || return 0.0
    # Each draw's rows carry the root of its relevance, so its squared error
    # already counts by it.
    return sqrt(sum(values(err2); init=0.0) / total)
end

"""
    _is_coreset_choose(cs, rd, start)

The coreset whose expected round cost (`_is_round_cost`) is least, refitted on
every row, or why there is none: too few rows to fit (`:held`), none expected
to cost less than a full round (`:cost`), or weights all zero (`:empty`).
Sizes are tried along one forward-selection path (`_coreset_select`, refitted
by NNLS each step), which grows until what a size costs with checks that never
fail is more than the best found: nothing larger can be cheaper. Two sets are
held out: the newest of the rows the sampler's rounds checked in full
(`focus`), which are where the weight is moving to, and every fourth of the
rest. Older checked rows are fitted like any other: held out for good, a
coreset fitted mostly on the first rounds' tails could never reach them. Held
out deterministically, so that the sampler's random stream is the same with a
coreset as without. Score rows are held out from none, and the held-out rates
read values, which is what a check reads.

The failure rate a held-out set shows carries its uncertainty: the set counts
as as many observations as its draws amount to under the check's own weights
(`_is_effective`), under a uniform prior. A check weighs its draws' errors by
their relevance normalised among themselves, so draws that all sit far below
the top inform it as fully as draws at the top. Counted instead by their
relevance mass, how many draws near the top they amount to, they counted as
none, and every size was priced at the prior's 0.5: the cheapest won and
failed, 27 of 28 choices at 1,000 units on the strong shell while DAIS's
draws, which vary in every direction, went through the coreset, and 11 or 12
of 13 for DAIS's own coreset at 250 units. A sampler's first rounds are tails
with a draw or two near the top between them, which amount to a draw or two
either way; a bar on that count kept the first coreset back for two to three
rounds of 100, a third of a run's full evaluations, where pricing an uncertain
rate lets the checks decide.
"""
function _is_coreset_choose(tr::_ISRows, rd::_ISRound, start::Vector{Int})
    nrows = tr.nrows
    rows = collect(1:nrows)
    focus = intersect(tr.focus, rows)
    newest = focus[max(1, length(focus) - rd.checks + 1):end]
    plain = setdiff(rows, newest)
    held = plain[[i % 4 == 0 for i in 1:length(plain)]]
    fitrows = setdiff(plain, held)
    top = _is_top(tr)
    rel = _is_relevance(tr.y, top)
    none(why; kw...) = (chosen=false, why=why, kw...)
    length(fitrows) < 10 && return none(:held; relheld=_is_effective(rel[held]))
    N = tr.N
    fullround = _is_full_round(rd)
    score = _is_score_rows(tr, top)
    # Every fourth scored draw is held out of the fit, as every fourth value
    # row is, for the size's force check (`_is_force_error`).
    zheld = score === nothing ? Set{Int}() : _is_held_draws(tr.zdraw)
    zin = [!(d in zheld) for d in tr.zdraw]
    scorefit = score === nothing ? nothing : (Z=score.Z[zin, :], z=score.z[zin])
    A, b, mx, my = _is_rows(tr, fitrows, rel[fitrows], scorefit)
    mheld = _is_effective(rel[held])
    mnewest = _is_effective(rel[newest])
    failrate(S, w, c, idx, m) = begin
        isempty(idx) && return NaN
        pred = [c + sum(w[i] * tr.L[S[i], j] for i in eachindex(S); init=0.0) for j in idx]
        (m * _is_fail_rate(pred .- tr.y[idx], rel[idx], pred, rd) + 1) / (m + 2)
    end
    nextsize(M) = M + max(1, M ÷ 10)
    sizes = Int[]
    M = max(1, length(start))
    while M <= N && _is_round_cost(M, N, rd, 0.0) < fullround
        push!(sizes, M)
        M = nextsize(M)
    end
    isempty(sizes) && return none(:cost; fullround=fullround)
    # The best size whatever it costs, for the record.
    priced = Ref((M=0, cost=Inf, fail=NaN))
    value = (units, w) -> begin
        Mu = length(units)
        bound = fullround - _is_round_cost(nextsize(Mu), N, rd, 0.0)
        # Weights that all come out zero leave a constant and no unit to
        # evaluate, which costs nothing by the formula: priced, such a
        # coreset was chosen, and the engine was then handed an empty unit
        # list.
        any(>(0), w) || return (0.0, bound, nothing)
        # Whether a coreset can follow the draws is the data's to say, not the
        # parameters': where the units depend on the parameters alike, a few
        # carry a sum over many parameters, and where each moves in its own
        # directions no fewer than those can (Campbell 2024 proves the bound
        # M > d under exactly that, the units' gradients independent). So the
        # size's gradient is held to the full sum's at held-out draws. Without
        # it the choice priced every size alike on thin evidence and took the
        # smallest, which failed: 1-16 units on DAIS's 33-parameter shell at
        # 1,000 units, and one unit a choice on a 266-parameter fit with 200.
        _is_force_error(score, tr, zheld, units, w, top) > _IS_CORESET_TOL.force &&
            return (0.0, bound, nothing)
        # Not `c`: the function's own `c` below would be this one, boxed.
        level = my - sum(mx[units] .* w)
        fail = maximum(filter(!isnan, [failrate(units, w, level, held, mheld),
            failrate(units, w, level, newest, mnewest)]); init=-Inf)
        cost = _is_round_cost(count(>(0), w), N, rd, fail)
        cost < priced[].cost && (priced[] = (M=Mu, cost=cost, fail=fail))
        (fullround - cost, bound, (M=Mu, cost=cost, fail=fail))
    end
    found = _coreset_select(A, b, ones(N), Inf, sizes, value; refit=true, start=start)
    best = priced[]
    found === nothing && return none(:cost; fullround=fullround,
        bestM=best.M, bestcost=best.cost, bestfail=best.fail)
    # Refitted on every row, at the size priced.
    Aall, ball, mxall, myall = _is_rows(tr, rows, rel, score)
    final = _coreset_select(Aall, ball, ones(N), Inf, [found.extra.M],
        (u, w) -> (1.0, 0.0, nothing); refit=true, start=start)
    final === nothing && return none(:empty)
    S = final.units
    w = final.weights
    any(>(0), w) || return none(:empty)
    c = myall - sum(mxall[S] .* w)
    # The coreset is the units it was priced by, those with weight: refits on
    # later rows re-weight them, and add none. Keeping the zero-weight units of
    # the path let refits give them weight, and coresets priced at 28 units of
    # 1,000 grew to 46.
    keep = w .> 0
    return (chosen=true, S=S[keep], w=w[keep], c=c, cost=found.extra.cost,
        fail=found.extra.fail, relheld=mheld)
end

"""A choice, made again only once the sampler has added at least as many
rows as the last coreset had units, one for each weight it fits: with fewer,
it refits to nearly the coreset that just failed. Rows, not seconds, so the
choice does not depend on the machine's load. `true` where one was chosen."""
function _is_choose!(cs::CTSEMISCoreset, start::Vector{Int})
    cs.has_cf && cs.rows.nrows - cs.rows_at_choose < _is_active_n(cs) && return false
    t = time()
    r = _is_coreset_choose(cs.rows, _is_round(cs), start)
    _is_tally!(cs, :choose, time() - t)
    cs.nchoose += 1
    cs.rows_at_choose = cs.rows.nrows
    code = Float64(something(findfirst(==(r.chosen ? :chosen : r.why),
        (:chosen, :held, :cost, :empty)), 0))
    get_(nt, k) = haskey(nt, k) ? Float64(nt[k]) : NaN
    push!(cs.trace, (Float64(cs.rows.nrows), r.chosen ? Float64(length(r.S)) : 0.0, code,
        r.chosen ? r.fail : get_(r, :bestfail), r.chosen ? r.cost : get_(r, :bestcost),
        get_(r, :fullround), get_(r, :bestM), get_(r, :relheld)))
    r.chosen || return false
    cs.has_cf = true
    cs.S = r.S
    cs.w = r.w
    cs.c = r.c
    cs.cf_cost = r.cost
    cs.cf_fail = r.fail
    push!(cs.sizes, _is_active_n(cs))
    cs.runs = 0
    cs.fails = 0
    cs.trusted = cs.last_passed
    return true
end

"""Whether the coreset in force is expected to cost more than full rounds."""
function _is_dearer(cs::CTSEMISCoreset)
    r = _is_round(cs)
    return _is_round_cost(_is_active_n(cs), cs.N, r, _is_failnow(cs)) >= _is_full_round(r)
end

"""A failed check: the coreset is chosen again from what the sampler has
taught it, and if none pays, full rounds until one does."""
function _is_failed!(cs::CTSEMISCoreset)
    cs.fails += 1
    cs.fails_total += 1
    cs.trusted = false
    cs.last_passed = false
    _is_choose!(cs, copy(cs.S)) || (cs.phase = :training)
    return nothing
end

"""A passed check: the weights refitted on every row; chosen again if they
all come out zero, and back to full rounds if its record makes it dearer."""
function _is_passed!(cs::CTSEMISCoreset)
    cs.trusted || (cs.proved += 1)
    cs.trusted = true
    cs.last_passed = true
    w, c = _is_weights(cs.rows, cs.S)
    if !any(>(0), w)
        _is_choose!(cs, copy(cs.S)) || (cs.phase = :training)
        return nothing
    end
    cs.w = w
    cs.c = c
    _is_dearer(cs) && (cs.phase = :training)
    return nothing
end

_is_cached(cs::CTSEMISCoreset, X) = !isempty(cs.lastgrad) &&
    all(haskey(cs.lastgrad, _is_key(view(X, :, j))) for j in 1:size(X, 2))

"""
A batch with no coreset to run on -- the first rounds, or the rest of one whose
check failed -- evaluated in full a half at a time, while the other half could
pay for a check on a coreset (one of the last coreset's size that never
fails). After each half a coreset is chosen from what the sampler has taught
so far, and the rest runs on it if one pays. Full rounds are what the sampler
pays without a coreset, so learning from them costs nothing; evaluated whole,
they and failed batches were a third to a half of a run's full evaluations.
"""
function _is_inparts!(cs::CTSEMISCoreset, X::AbstractMatrix, why::Symbol)
    n = size(X, 2)
    v = fill(NaN, n)
    done = 0
    allcached = _is_cached(cs, X)
    while done < n
        left = n - done
        half = left ÷ 2
        r0 = _is_round(cs)
        r = _ISRound(half, r0.grads, r0.checks)
        bc = cs.rows.nrows == 0 ? _is_batch_cost(0, 1, r, 0.0) :
            _is_batch_cost(cs.has_cf ? _is_active_n(cs) : 0, cs.N, r, 0.0)
        # Values a gradient batch already gave cost nothing in full.
        take = !allcached && half > 0 && bc.coreset < bc.full ? left - half : left
        idx = done .+ (1:take)
        t0 = _is_full!(cs, X[:, idx], why)
        v[idx] .= t0.value
        _is_learn!(cs.rows, t0)
        done += take
        if cs.rows.nrows > 0 && _is_choose!(cs, cs.has_cf ? copy(cs.S) : Int[])
            cs.phase = :coreset
            cs.switches += 1
            if done < n
                rest = (done + 1):n
                v[rest] .= _is_withcoreset!(cs, X[:, rest])
            end
            return v
        end
    end
    return v
end

"""A batch while there is a coreset: on it where that pays for the batch's
size, else in full."""
function _is_withcoreset!(cs::CTSEMISCoreset, X::AbstractMatrix)
    r0 = _is_round(cs)
    r = _ISRound(size(X, 2), r0.grads, r0.checks)
    bc = _is_batch_cost(_is_active_n(cs), cs.N, r, _is_failnow(cs))
    bc.coreset < bc.full && return _is_oncoreset!(cs, X)
    return _is_freecheck!(cs, X, :small)
end

"""A batch in full while there is a coreset -- one too small to pay for the
coreset's check, or whose values a gradient batch already gave -- so the check
is free: the coreset's prediction from the unit terms, at the draws a check
would read."""
function _is_freecheck!(cs::CTSEMISCoreset, X::AbstractMatrix, why::Symbol)
    n = size(X, 2)
    t0 = _is_full!(cs, X, why)
    pred = _is_predict(cs, t0.L)
    pick = _is_checkpick([isfinite(p) ? p : -Inf for p in pred])
    good = [isfinite(t0.total[j]) && isfinite(pred[j]) for j in pick]
    _is_learn!(cs.rows, t0, pick; focus=true)
    _is_learn!(cs.rows, t0, setdiff(1:n, pick))
    if any(good)
        # Against the units' sum, which the prediction stands in for; the log
        # posterior adds the prior, which a coreset's evaluation adds itself.
        err = (pred[pick] .- t0.total[pick])[good]
        append!(cs.checks, err)
        cs.runs += 1
        cs.free_total += 1
        rel = _is_relevance(t0.total[pick][good], _is_top(cs.rows))
        _is_passes(err, rel, _IS_CORESET_TOL.drift) ? _is_passed!(cs) : _is_failed!(cs)
    end
    return t0.value
end

"""A batch on the coreset, checked at a few of its draws evaluated in full."""
function _is_oncoreset!(cs::CTSEMISCoreset, X::AbstractMatrix)
    n = size(X, 2)
    v = _is_coreset_value!(cs, X)
    pick = _is_checkpick(v)
    t0 = _is_full!(cs, X[:, pick], :check)
    good = [isfinite(t0.value[i]) && isfinite(v[pick[i]]) for i in eachindex(pick)]
    if any(good)
        err = (v[pick] .- t0.value)[good]
        append!(cs.checks, err)
        v[pick[good]] .= t0.value[good]
        _is_learn!(cs.rows, t0; focus=true)
        rel = _is_relevance(t0.total[good], _is_top(cs.rows))
        cs.runs += 1
        cs.runs_total += 1
        if !_is_passes(err, rel, _IS_CORESET_TOL.drift)
            # Failed. None of the round's values comes from this coreset: the
            # sampler keeps every value it is handed, and a coreset that
            # overrated the draws it ranked highest by 20 log units (a first,
            # too small, coreset on the 250-subject strong shell) misled its
            # adaptation for thirty rounds and left 1,900 draws for the end
            # check to correct. The rest is evaluated as a batch with no
            # coreset is: its full halves are training data where the sampler
            # now is, and a coreset chosen from them takes what is left,
            # checked again.
            cs.fails += 1
            cs.fails_total += 1
            cs.trusted = false
            cs.last_passed = false
            cs.phase = :training
            others = setdiff(1:n, pick)
            if !isempty(others)
                v[others] .= _is_inparts!(cs, X[:, others], :failed_rest)
            elseif _is_choose!(cs, copy(cs.S))
                cs.phase = :coreset
            end
        else
            _is_passed!(cs)
        end
    end
    return v
end

"""
    ctsem_is_coreset_batch!(cs, X)

The density at the columns of `X`, as the sampler's value batch: in full, or on
the coreset where that pays, by the rules above. Unusable draws `-1e100`.
"""
function ctsem_is_coreset_batch!(cs::CTSEMISCoreset, X::AbstractMatrix)
    X = Matrix{Float64}(X)
    n = size(X, 2)
    # A batch's size, and the gradients asked for since the last one (none
    # before the first: until then a round's gradients are taken as 30% of
    # its draws, about what 'checked''s Newton components ask for).
    cs.values = _is_smooth(cs.values, n)
    cs.batches > 0 && (cs.grads = _is_smooth(cs.grads, cs.gradsnow))
    cs.batches += 1
    cs.gradsnow = 0.0
    _is_batch_start!(cs.rows)
    v = cs.phase != :coreset ? _is_inparts!(cs, X, :training) :
        _is_cached(cs, X) ? _is_freecheck!(cs, X, :small) : _is_withcoreset!(cs, X)
    return [isfinite(x) ? x : -1e100 for x in v]
end

"""
    ctsem_is_coreset_fullgrad!(cs, X)

A full gradient batch at the columns of `X` that teaches the coreset: the
batch's unit terms kept for a value batch at the same draws, and the scores of
a few of its draws learned. `(value, gradient)`, the gradient `npar x n`.
"""
function ctsem_is_coreset_fullgrad!(cs::CTSEMISCoreset, X::AbstractMatrix)
    X = Matrix{Float64}(X)
    t = time()
    n = size(X, 2)
    cs.cost += 2 * n
    res = ctsem_coreset_terms_gradient(cs.laplace, X, _IS_SCORE_DRAWS.per)
    pick = res.scored
    cs.lastgrad = Dict(_is_key(view(X, :, j)) => j for j in 1:n)
    cs.lastgrad_value = res.value
    cs.lastgrad_total = res.total
    cs.lastgrad_L = res.unit_loglik
    for j in 1:n
        isfinite(res.value[j]) && (cs.known[_is_key(view(X, :, j))] = res.value[j])
    end
    _is_learn_scores!(cs.rows, X[:, pick], res.scores, res.total[pick])
    _is_tally!(cs, :gradient_full, time() - t, n)
    value = [isfinite(x) ? x : -1e100 for x in res.value]
    return (value=value, gradient=res.gradient)
end

"""Each scored draw's rows: per parameter, the units' scores and their sum;
the newest `_IS_SCORE_DRAWS.kept` draws' rows are kept."""
function _is_learn_scores!(tr::_ISRows, Xs::AbstractMatrix, S::AbstractArray{<:Real,3},
        totals::AbstractVector)
    p = size(S, 2)
    for d in 1:size(Xs, 2)
        Sd = view(S, :, :, d)
        (isfinite(totals[d]) && all(isfinite, Sd)) || continue
        tr.scoredraws += 1
        tr.Z = vcat(tr.Z, Matrix(transpose(Sd)))
        append!(tr.zsum, vec(sum(Sd; dims=1)))
        append!(tr.zpar, 1:p)
        append!(tr.ztotal, fill(totals[d], p))
        append!(tr.zdraw, fill(tr.scoredraws, p))
        push!(tr.zdraws, Vector{Float64}(Xs[:, d]))
    end
    if tr.scoredraws > _IS_SCORE_DRAWS.kept
        keep = tr.zdraw .> tr.scoredraws - _IS_SCORE_DRAWS.kept
        tr.Z = tr.Z[keep, :]
        tr.zsum = tr.zsum[keep]
        tr.zpar = tr.zpar[keep]
        tr.ztotal = tr.ztotal[keep]
        tr.zdraw = tr.zdraw[keep]
        tr.zdraws = tr.zdraws[max(1, end - _IS_SCORE_DRAWS.kept + 1):end]
    end
    return nothing
end

"""
    ctsem_is_coreset_gradbatch!(cs, X; check=false)

The sampler's gradient batch: on the coreset where it is in force, and in full
otherwise. A coreset not yet trusted here has its gradient batch checked as a
batch of values is: its values at a check's draws evaluated in full, and the
batch's gradients taken in full if they fail. Where that is not expected to
cost less than the batch in full -- a few gradients, or a coreset that keeps
failing -- they are taken in full straight away. `check` checks every batch
so, trusted or not: for a caller with no value batches of its own to check
the coreset with, such as DAIS, whose every batch is values and gradients.
Without it a coreset chosen from DAIS's first draws, three units of 200,
served the rest of its run unchecked as it moved, and its answer was 0.9
posterior sds from the full density's.
"""
function ctsem_is_coreset_gradbatch!(cs::CTSEMISCoreset, X::AbstractMatrix;
    check::Bool=false)
    X = Matrix{Float64}(X)
    n = size(X, 2)
    cs.gradsnow += n
    _is_batch_start!(cs.rows)
    cs.phase != :coreset && return ctsem_is_coreset_fullgrad!(cs, X)
    units = cs.S[cs.w .> 0]
    weights = cs.w[cs.w .> 0]
    share = length(units) / cs.N
    k = min(8, n)
    checking = check || !cs.trusted
    checking && 2 * share * n + k + _is_failnow(cs) * 2 * n >= 2 * n &&
        return ctsem_is_coreset_fullgrad!(cs, X)
    t = time()
    out = ctsem_coreset_batch_gradient(cs.laplace, X, units, weights)
    _is_tally!(cs, :gradient_coreset, time() - t, n)
    cs.cost += 2 * n * share
    value = out[1, :] .+ cs.c
    if checking
        pick = _is_checkpick([isfinite(x) ? x : -Inf for x in value])
        t0 = _is_full!(cs, X[:, pick], :gradient_check)
        good = [isfinite(t0.value[i]) && isfinite(value[pick[i]]) for i in eachindex(pick)]
        if any(good)
            err = (value[pick] .- t0.value)[good]
            append!(cs.checks, err)
            _is_learn!(cs.rows, t0; focus=true)
            cs.runs += 1
            rel = _is_relevance(t0.total[good], _is_top(cs.rows))
            if !_is_passes(err, rel, _IS_CORESET_TOL.drift)
                _is_failed!(cs)
                return ctsem_is_coreset_fullgrad!(cs, X)
            end
            value[pick[good]] .= t0.value[good]
            _is_passed!(cs)
        end
    end
    return (value=[isfinite(x) ? x : -1e100 for x in value], gradient=out[2:end, :])
end

"""The full density at the columns of `X`: known values where a draw was
evaluated in full already, the rest evaluated now and counted for the end
check."""
function _is_fullvalues!(cs::CTSEMISCoreset, X::AbstractMatrix)
    n = size(X, 2)
    out = [get(cs.known, _is_key(view(X, :, j)), NaN) for j in 1:n]
    todo = findall(isnan, out)
    if !isempty(todo)
        t = time()
        out[todo] .= ctsem_evaluate_batch(cs.laplace, X[:, todo])
        cs.cost += length(todo)
        cs.why[:end] += length(todo)
        _is_tally!(cs, :full, time() - t, length(todo))
    end
    return out
end

"""
    ctsem_is_coreset_finalize!(cs, samples, logp, logw)

The end check. The full density where the weights sit; if the coreset was off
there by more than the tolerance either side of its typical offset (a constant
cancels), every draw carrying the weight is corrected to the full density, so
the weights handed back are exact where they matter. The corrected log
weights, or `nothing` where none needed correcting. `samples` is `npar x n`.
"""
function ctsem_is_coreset_finalize!(cs::CTSEMISCoreset, samples::AbstractMatrix,
        logp::AbstractVector, logw::AbstractVector)
    cs.has_cf || return nothing
    samples = Matrix{Float64}(samples)
    o = sortperm(-logw)
    top = o[1:min(20, length(o))]
    vt = _is_fullvalues!(cs, samples[:, top])
    delta = vt .- logp[top]
    ok = isfinite.(delta)
    offset = any(ok) ? median(delta[ok]) : 0.0
    worst = any(ok) ? maximum(abs.(delta[ok] .- offset)) : Inf
    cs.has_final = true
    cs.final_checked = length(top)
    cs.final_worst = worst
    cs.final_corrected = 0
    worst <= _IS_CORESET_TOL.final && return nothing
    mx = maximum(logw)
    cw = cumsum(exp.(logw[o] .- mx)) ./ sum(exp.(logw .- mx))
    need = union(o[1:something(findfirst(>=(0.999), cw), length(o))],
        findall(>(mx + log(1e-3)), logw))
    need = setdiff(need, top)
    full_all = fill(NaN, length(logw))
    full_all[top] .= vt
    isempty(need) || (full_all[need] .= _is_fullvalues!(cs, samples[:, need]))
    d = full_all .- logp
    d[.!isfinite.(d)] .= offset
    cs.final_corrected = length(need) + length(top)
    return logw .+ d
end

"""
    ctsem_is_coreset_record(cs)

What the coreset did, for the importance sampler's record of it in R (its
details' `importance_sampling`, `coreset`).
"""
function ctsem_is_coreset_record(cs::CTSEMISCoreset)
    order(d, ks) = (; (k => d[k] for k in ks)...)
    choices = isempty(cs.trace) ? zeros(0, 8) :
        [cs.trace[i][j] for i in eachindex(cs.trace), j in 1:8]
    return (used=cs.has_cf, units=cs.rows.nrows > 0 ? cs.N : 0,
        size=cs.has_cf ? _is_active_n(cs) : 0, sizes=cs.sizes, phase=String(cs.phase),
        switches=cs.switches, rounds_on_coreset=cs.runs_total, failed_checks=cs.fails_total,
        expected_round_cost=cs.has_cf ? cs.cf_cost : NaN, free_checks=cs.free_total,
        chosen=_is_chosen(cs), attempts=cs.nchoose, chosen_passed=cs.proved,
        seconds=order(cs.secs, (:full, :coreset, :gradient_full, :gradient_coreset, :choose)),
        draws=order(cs.draws, (:full, :coreset, :gradient_full, :gradient_coreset)),
        full_for=order(cs.why, (:training, :small, :check, :gradient_check, :failed_rest, :end)),
        choices=choices, training_draws=cs.rows.nrows, score_draws=cs.rows.scoredraws,
        reused_values=cs.reused,
        check_rms=isempty(cs.checks) ? NaN : sqrt(sum(abs2, cs.checks) / length(cs.checks)),
        check_max=isempty(cs.checks) ? NaN : maximum(abs, cs.checks),
        has_final=cs.has_final, final_checked=cs.final_checked, final_worst=cs.final_worst,
        final_corrected=cs.final_corrected, cost=cs.cost)
end
