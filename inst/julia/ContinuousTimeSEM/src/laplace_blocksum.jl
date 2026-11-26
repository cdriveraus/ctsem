# A nested unit's blocks combine their own maxima.
#
# A unit of persons in a study is one integral, and Laplace scores it at one
# maximum, `z*`. Where a person's effects have a second maximum the whole-unit
# mode sum (`_laplace_mode_sum`) kept a maximum of the whole unit for each such
# person: K maxima for K persons, each a pass over every member, and a valley
# between every pair, so an evaluation cost members times K squared (305 s on
# a study of 200 persons, 190 of them with two maxima). It also handed a
# fading person's share to whichever maximum ranked next, another person's.
#
# Here each block of the unit below the top keeps its own other maxima: its
# subtree (the block and the blocks under it) at another maximum of the
# subtree's objective `h_b`, its members' log likelihood and its coordinates'
# prior, with the rest of the unit held at `z*`. The block's configurations,
# its own at `z*` and those, combine as the mode sum combines a unit's
# maxima, and the unit's term is
#
#     T = T(z*) + sum_b D_b,   D_b = sum_c pi_bc a_bc,
#
# `pi` the shares `_laplace_mode_shares(f, a)`, `f_c` each configuration's fade
# (its subtree curvature's eigenvalues under `_LAPLACE_MODE_FADE`, its valleys
# towards the block's other configurations under `_LAPLACE_MODE_DEPTH`, along
# the subtree's segments), and `a_c` the configuration's value against `z*`'s
# (zero for the own one). Everything a block reads is its own members': a
# block whose own maximum is the best at full weight adds nothing, and costs a
# Newton step and a valley among its own members.
#
# The value. Held at `z*`, an alternative's conditional Laplace term
# `tau_c = h_b(y_c) - logdet(M_SS(y_c)) / 2` differs from the whole unit's at
# the maximum where the rest has moved with it, `z_c`: on persons in two
# studies by 0.003-0.075 nats at 20 persons a study, 0.003 median at 100 (the
# rest's gain from moving). Where a hand-over or a change of rank happens --
# configurations tied, or one fading -- the value must be the whole unit's, or
# the term jumps by that much when the inner solve moves `z*` to `z_c`. So
#
#     a_c = (1 - s_c) (tau_c - tau_own) + s_c (T(z_c) - T(z*)),
#     s_c = 1 - q_c prod_c' f_c',   q_c = step((|tau_c - tau_own| - d) / d),
#
# `d` `_LAPLACE_MODE_DEPTH` and `step` the fade's smoothstep: the whole unit's
# value wherever any configuration of the block fades or two are within `d`
# nats, the conditional one where they are all at full weight and apart, where
# only the best's value counts and is exact either way. `z_c` is solved, a pass
# over the unit, only where `s_c > 0`. The values ranked are the values
# combined, so the combination is continuous where ranks swap, as the mode sum
# is.
#
# The gradient takes each piece's: the conditional terms', the fades' and the
# valleys' through the subtree's own members, each configuration's position
# by the implicit function theorem on the subtree (a block solve), what that
# leaves on the rest of the unit through `z*` by one solve with the unit's
# curvature, and `T(z_c)`'s by the unit's own route at `z_c`.

"""
    _laplace_subtree(blocks, b)

Block `b`'s subtree in a unit laid out as `blocks`: its blocks in the unit's
order (`S`, the block last), their coordinates in the unit (`idx`), the blocks
relabelled as the subtree's own (`sub`, ancestors outside it dropped), where
each unit coordinate lies in it (`loc`, 0 outside), and its members
(`members`, block `b`'s).
"""
function _laplace_subtree(blocks::Vector{CTSEMLaplaceBlock}, b::Int)
    nb = length(blocks)
    S = Int[t for t in 1:nb if t == b || b in blocks[t].ancestors]
    pos = zeros(Int, nb)
    for (k, t) in enumerate(S)
        pos[t] = k
    end
    sub = Vector{CTSEMLaplaceBlock}(undef, length(S))
    idx = Int[]
    off = 0
    for (k, t) in enumerate(S)
        sub[k] = CTSEMLaplaceBlock(off, blocks[t].size, blocks[t].members,
            blocks[t].level, Int[pos[a] for a in blocks[t].ancestors if pos[a] > 0])
        append!(idx, (blocks[t].offset + 1):(blocks[t].offset + blocks[t].size))
        off += blocks[t].size
    end
    d = sum(k.size for k in blocks; init=0)
    loc = zeros(Int, d)
    for (k, i) in enumerate(idx)
        loc[i] = k
    end
    return (b=b, S=S, pos=pos, sub=sub, idx=idx, loc=loc, members=blocks[b].members)
end

"""The subtree's block of a unit's block curvature `M`, copied (the unit's
curvature is the slot's scratch) and symmetrised."""
function _laplace_subtree_matrix(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}, st) where {T}
    A = CTSEMBlockMatrix{T}(Matrix{T}[copy(M.diag[t]) for t in st.S],
        Vector{Matrix{T}}[Matrix{T}[copy(M.coupling[t][r])
            for (r, a) in enumerate(blocks[t].ancestors) if st.pos[a] > 0] for t in st.S])
    for D in A.diag
        _laplace_symmetrise!(D)
    end
    return A
end

"""
    _laplace_subtree_curvature(laplace, U, values, Ls, u, st)

The curvature of the subtree's objective `h_b` at the unit's coordinates `u`
(`_laplace_block_sum`), on the subtree's blocks: the identity less its members'
Hessian in its coordinates, by the members' sweeps over its columns alone
(`_laplace_columns_jacobian`).
"""
function _laplace_subtree_curvature(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    st) where {T}
    blocks = laplace.units.blocks[U]
    Jw = _laplace_columns_jacobian(laplace, U, values, Ls, u, copy(st.idx),
        copy(st.members), _laplace_hmode(laplace), nothing)
    diag = Matrix{T}[]
    coup = Vector{Matrix{T}}[]
    for (k, t) in enumerate(st.S)
        kt = blocks[t].size
        ot = blocks[t].offset
        ct = st.sub[k].offset
        D = Matrix{T}(undef, kt, kt)
        @inbounds for j in 1:kt, i in 1:kt
            D[i, j] = (i == j ? one(T) : zero(T)) - (Jw[ot + i, ct + j] + Jw[ot + j, ct + i]) / 2
        end
        push!(diag, D)
        row = Matrix{T}[]
        for a in blocks[t].ancestors
            st.pos[a] > 0 || continue
            ka = blocks[a].size
            oa = blocks[a].offset
            ca = st.sub[st.pos[a]].offset
            C = Matrix{T}(undef, kt, ka)
            @inbounds for j in 1:ka, i in 1:kt
                C[i, j] = -(Jw[ot + i, ca + j] + Jw[oa + j, ct + i]) / 2
            end
            push!(row, C)
        end
        push!(coup, row)
    end
    return CTSEMBlockMatrix{T}(diag, coup)
end

"""The subtree's objective `h_b` at the unit's coordinates `u`: its members' log
likelihood less its coordinates' prior."""
function _laplace_subtree_value(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    aws, st)
    v = _laplace_unit_loglik_value(laplace, U, values, Ls, u, aws, Float64[];
        positions=st.members).value
    isfinite(v) || return -Inf
    for i in st.idx
        v -= u[i]^2 / 2
    end
    return v
end

"""`h_b` and its gradient in the subtree's coordinates, at `u`."""
function _laplace_subtree_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    aws, st)
    r = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws, st.members)
    n = length(st.idx)
    isfinite(r.value) || return (value=-Inf, g=fill(NaN, n))
    g = Vector{Float64}(undef, n)
    v = Float64(r.value)
    @inbounds for (k, i) in enumerate(st.idx)
        g[k] = r.gradient[i] - u[i]
        v -= u[i]^2 / 2
    end
    return (value=v, g=g)
end

"""
    _laplace_block_newton(laplace, U, values, Ls, aws, base, st, y0)

Newton on the subtree's objective `h_b` from `y0`, the rest of the unit held at
`base`: the unit solve's rules (`_laplace_newton_unit_mode`) on the subtree's
members alone -- stopped on the exact step under `inner_tol`, at the value's
floor by whether a step halves the gradient (`_laplace_floor_step`), a curvature
that does not factor shifted until it does, a step that loses halved. Returns
the subtree's coordinates `y`, `h_b` there, whether it converged, and the
curvature there with its factorization (`M`, `F`).
"""
function _laplace_block_newton(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64}, st, y0::AbstractVector{Float64})
    idx = st.idx
    fail = (y=collect(Float64, y0), value=-Inf, converged=false, M=nothing, F=nothing)
    u = copy(base)
    u[idx] .= y0
    cur = _laplace_subtree_gradient(laplace, U, values, Ls, u, aws, st)
    isfinite(cur.value) || return fail
    for _ in 1:laplace.inner_maxiter
        M = _laplace_subtree_curvature(laplace, U, values, Ls, u, st)
        all(D -> all(isfinite, D), M.diag) || return fail
        F = _soft_shifted_factor(M, st.sub, 0.0)
        exact = F !== nothing
        if !exact
            scale = maximum((maximum(abs, D) for D in M.diag); init=1.0)
            shift = sqrt(eps(Float64)) * scale
            for _ in 1:30
                F = _soft_shifted_factor(M, st.sub, -shift)
                F === nothing || break
                shift *= 10
            end
            F === nothing && return fail
        end
        step = _soft_solve(F, st.sub, cur.g)
        if exact && maximum(abs, step; init=0.0) < laplace.inner_tol
            return (y=u[idx], value=cur.value, converged=true, M=M, F=F)
        end
        if exact && _laplace_floor_step(laplace, cur.g, step, cur.value)
            utry = copy(u)
            utry[idx] .+= step
            nxt = _laplace_subtree_gradient(laplace, U, values, Ls, utry, aws, st)
            if isfinite(nxt.value) && _laplace_floor_step_taken(cur.g, nxt.g)
                u, cur = utry, nxt
                continue
            end
            return (y=u[idx], value=cur.value, converged=true, M=M, F=F)
        end
        # A step too short to move the point moves it at no scale, and the
        # next iteration takes the same step: stuck, as in
        # `_laplace_newton_unit_mode`, where the line search accepted such a
        # step at every iteration to the cap. A loop, not `all` over a
        # closure: `u` is assigned again below, and a closure reading it
        # would box it.
        still = true
        for k in eachindex(idx)
            if u[idx[k]] + step[k] != u[idx[k]]
                still = false
                break
            end
        end
        if still
            exact && return (y=u[idx], value=cur.value, converged=true, M=M, F=F)
            return fail
        end
        t = 1.0
        moved = false
        for _ in 1:30
            utry = copy(u)
            utry[idx] .+= t .* step
            nxt = _laplace_subtree_gradient(laplace, U, values, Ls, utry, aws, st)
            if isfinite(nxt.value) && nxt.value >= cur.value
                u, cur = utry, nxt
                moved = true
                break
            end
            t /= 2
        end
        if !moved
            # Nothing the value resolves is left in the step's direction.
            exact && return (y=u[idx], value=cur.value, converged=true, M=M, F=F)
            return fail
        end
    end
    return fail
end

# `h_b` along the subtree's segment from `y1` to `y2`, the rest at `base`:
# `-Inf` where it is not finite.
function _laplace_block_segment(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64}, st, y1::Vector{Float64}, y2::Vector{Float64})
    w = copy(base)
    return t -> begin
        @inbounds for (k, i) in enumerate(st.idx)
            w[i] = y1[k] + t * (y2[k] - y1[k])
        end
        _laplace_subtree_value(laplace, U, values, Ls, w, aws, st)
    end
end

"""Whether two of a block's configurations are separate maxima of `h_b`: it
dips below the lower somewhere on the segment between them (seven interior
points), as `_laplace_modes_separated` asks of a unit's."""
function _laplace_block_separated(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64}, st, y1::Vector{Float64}, h1::Float64, y2::Vector{Float64},
    h2::Float64)
    lower = min(h1, h2)
    margin = 1e-8 * (1 + abs(lower))
    h = _laplace_block_segment(laplace, U, values, Ls, aws, base, st, y1, y2)
    # The middle first: two separate maxima usually dip there, and one point
    # below the lower settles it.
    for k in (4, 2, 6, 1, 3, 5, 7)
        g = h(k / 8)
        (isfinite(g) && g >= lower - margin) || return true
    end
    return false
end

"""
    _laplace_block_join!(configs, laplace, U, values, Ls, aws, base, st, r)

Block Newton's end point `r` into a block's configurations `configs` (the own
one first): a point the own one or another reaches without a dip is that one,
kept where it is the higher (never in place of the own one); one a valley
parts from all of them joins. Returns whether it joined.
"""
function _laplace_block_join!(configs::Vector{Any}, laplace::CTSEMLaplaceObjective,
    U::Integer, values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64}, st, r)
    for k in eachindex(configs)
        c = configs[k]
        _laplace_maxabsdiff(c.y, r.y) < 1e-4 && return false
        if !_laplace_block_separated(laplace, U, values, Ls, aws, base, st, c.y, c.h,
                r.y, r.value)
            k > 1 && r.value > c.h && (configs[k] = (y=r.y, h=r.value, M=r.M, F=r.F))
            return false
        end
    end
    push!(configs, (y=r.y, h=r.value, M=r.M, F=r.F))
    return true
end

"""
    _laplace_block_expand(blocks, base, e)

A block's other maximum `e` as a maximum of the whole unit laid out as `blocks`:
`base` with `e`'s subtree at `e.y`. How a unit's blocks' maxima go wherever the
unit's maxima go as vectors (a checkpoint, the fit's modes, the known modes);
`_laplace_block_classify` reads them back.
"""
function _laplace_block_expand(blocks::Vector{CTSEMLaplaceBlock}, base::Vector{Float64},
    e::CTSEMBlockMode)
    st = _laplace_subtree(blocks, e.b)
    v = copy(base)
    length(e.y) == length(st.idx) && (v[st.idx] .= e.y)
    return v
end

"""
    _laplace_block_classify(blocks, base, v; tol)

A unit's maximum `v`, read against its mode `base`: where it differs from
`base` by more than `tol` only within one block's subtree below the top, that
block's other maximum (`_laplace_block_expand`'s inverse), otherwise `nothing`
-- a maximum of the whole unit, or `base` itself. The subtree is the lowest
holding every block that differs (blocks come before their ancestors).
"""
function _laplace_block_classify(blocks::Vector{CTSEMLaplaceBlock},
    base::Vector{Float64}, v::AbstractVector{Float64}; tol::Float64=1e-6)
    length(v) == length(base) || return nothing
    differ = Int[]
    for (t, blk) in enumerate(blocks)
        o = blk.offset
        any(i -> abs(v[o + i] - base[o + i]) > tol, 1:blk.size) && push!(differ, t)
    end
    isempty(differ) && return nothing
    for t in eachindex(blocks)
        isempty(blocks[t].ancestors) && return nothing
        all(s -> s == t || t in blocks[s].ancestors, differ) || continue
        st = _laplace_subtree(blocks, t)
        return CTSEMBlockMode(t, collect(Float64, v[st.idx]))
    end
    return nothing
end

"""Whether two of a unit's blocks' other maxima are the same: one block, and
coordinates within `1e-4`, as the mode sum tells maxima apart."""
_laplace_block_same(e::CTSEMBlockMode, f::CTSEMBlockMode) = e.b == f.b &&
    length(e.y) == length(f.y) && _laplace_maxabsdiff(e.y, f.y) < 1e-4

"""
    _laplace_block_project!(blockset, laplace, U, values, Ls, aws, base, ru)

A nested unit's blocks' other maxima from a maximum `ru` of the whole unit
that a search reached from elsewhere (a known start, the top block's probes):
each block below the top whose own coordinates differ from `base`'s by a tenth
of a prior unit or more, its subtree solved from `ru`'s with the rest held at
`base` (`_laplace_block_newton`), kept in `blockset[b]` where a valley parts it
from `base`'s. Returns whether any did: `ru` is then persons at their other
maxima, which the blocks' combinations credit one by one, and kept as a
maximum of the whole unit as well it was credited twice -- ranked by its own
term beside the blocks' credits at the unit's mode, it took the lead as the
persons' preferences turned, and the term fell by the credits, 3 nats on
persons with sign-ambiguous loadings in studies of six.
"""
function _laplace_block_project!(blockset, laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64}, ru::AbstractVector{Float64})
    blocks = laplace.units.blocks[U]
    # A maximum whose rest stands apart from `base`'s (`_laplace_rest_shift`)
    # is the whole unit's, a study's effects in another basin with its persons
    # following: given the rest at `base` they fall back to their own. Solved
    # person by person, as until 2026-11-25, that was most of a search on the
    # affect pilot's study units (20 persons' subtree Newtons for a maximum
    # with the study 2.8 away, every one back at its own), and the maximum was
    # then dropped where a search reached it from a known start.
    _laplace_rest_shift(blocks, base, ru) >= _LAPLACE_REST_SHIFT && return false
    any_ = false
    for (b, blk) in enumerate(blocks)
        isempty(blk.ancestors) && continue
        o = blk.offset
        maximum((abs(ru[o + i] - base[o + i]) for i in 1:blk.size); init=0.0) >= 0.1 ||
            continue
        st = _laplace_subtree(blocks, b)
        r = _laplace_block_newton(laplace, U, values, Ls, aws, base, st, ru[st.idx])
        (r.converged && isfinite(r.value)) || continue
        ho = _laplace_subtree_value(laplace, U, values, Ls, base, aws, st)
        isfinite(ho) || continue
        _laplace_block_separated(laplace, U, values, Ls, aws, base, st, base[st.idx], ho,
            r.y, r.value) || continue
        push!(blockset[b], r.y)
        any_ = true
    end
    return any_
end

"""
How far, in prior units, a maximum's rest -- every block another lies under, a
study's own coordinates -- may stand from the unit's mode for it to be read
as persons at their other maxima (`_laplace_block_project!`): persons
switching move the study's effects by a little (0.21 for nine of twenty on the
affect pilot), a study's other basin by a lot (1.2-2.8 there).
"""
const _LAPLACE_REST_SHIFT = 0.5

"""The largest difference between `v` and `base` in a nested unit's rest, every
block another lies under (`_LAPLACE_REST_SHIFT`)."""
function _laplace_rest_shift(blocks::Vector{CTSEMLaplaceBlock}, base::AbstractVector,
    v::AbstractVector)
    under = fill(false, length(blocks))
    for block in blocks, a in block.ancestors
        under[a] = true
    end
    s = 0.0
    for (t, blk) in enumerate(blocks)
        under[t] || continue
        @inbounds for i in 1:blk.size
            s = max(s, abs(v[blk.offset + i] - base[blk.offset + i]))
        end
    end
    return s
end

"""
    _laplace_block_settle!(blockset, laplace, U, values, Ls, aws, best)

A nested unit's mode with each block at its configuration of best conditional
term (`_laplace_block_sum`'s `tau`), given the rest: every block with other
maxima in `blockset` takes the best of them where it beats its own, its own
kept beside it, and the unit's Newton then makes that the whole unit's maximum
(`_laplace_whole_mode_record`); repeated while a block changes, at most
`_LAPLACE_SETTLE_ROUNDS` times. Persons' choices move the study's effects, and
those move every other person's: on persons with sign-ambiguous loadings in
studies of six, three persons switching left two others whose other maximum
was then the better by 0.26 nats. Returns the unit's mode record.
"""
function _laplace_block_settle!(blockset, laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws, best)
    blocks = laplace.units.blocks[U]
    for _ in 1:_LAPLACE_SETTLE_ROUNDS
        base = collect(Float64, best.u)
        changed = false
        for (b, blk) in enumerate(blocks)
            (isempty(blk.ancestors) || isempty(blockset[b])) && continue
            st = _laplace_subtree(blocks, b)
            n = length(st.idx)
            ho = _laplace_subtree_value(laplace, U, values, Ls, base, aws, st)
            isfinite(ho) || continue
            Fo = _soft_shifted_factor(_laplace_subtree_curvature(laplace, U, values, Ls,
                base, st), st.sub, 0.0)
            Fo === nothing && continue
            top = ho - Fo.logdet / 2
            pick = 0
            picked = nothing
            for (k, y) in enumerate(blockset[b])
                length(y) == n || continue
                r = _laplace_block_newton(laplace, U, values, Ls, aws, base, st, y)
                (r.converged && isfinite(r.value) && r.F !== nothing) || continue
                tau = r.value - r.F.logdet / 2
                if tau > top + 1e-8
                    top, pick, picked = tau, k, r.y
                end
            end
            pick == 0 && continue
            own = base[st.idx]
            blockset[b] = vcat(Vector{Float64}[y for (k, y) in enumerate(blockset[b])
                if k != pick], [own])
            base[st.idx] .= picked
            changed = true
        end
        changed || break
        rn = _laplace_whole_mode_record(laplace, U, values, Ls, aws, base, best)
        rn === nothing && break
        best = rn
    end
    return best
end

"""How many rounds `_laplace_block_settle!` takes at most: each a pass over the
blocks with other maxima and the unit's Newton."""
const _LAPLACE_SETTLE_ROUNDS = 5

"""
    _laplace_block_partners!(blockset, laplace, U, values, Ls, aws, base)

`_laplace_partner_modes!` for a nested unit's blocks: past each block's
configuration (at `base`, and those in `blockset[b]`) whose subtree curvature
has an eigenvalue under `_LAPLACE_MODE_FADE`, the configuration the subtree's
Newton reaches 0.25 to 2 prior units either way along that eigenvector, given
the rest at `base`, joined to `blockset[b]` where a valley parts it from the
block's others. Past the fold the inner solve falls into one of them, and the
block's combination hands over to it only if it holds it.
"""
function _laplace_block_partners!(blockset, laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    base::Vector{Float64})
    blocks = laplace.units.blocks[U]
    # The base's curvature once, a pass over the unit, its blocks each block's
    # own configuration's.
    Mb = _laplace_copy_blockmatrix(_laplace_unit_curvature(laplace, U, values, Ls,
        copy(base)))
    for (b, blk) in enumerate(blocks)
        isempty(blk.ancestors) && continue
        st = _laplace_subtree(blocks, b)
        n = length(st.idx)
        Mo = _laplace_subtree_matrix(Mb, blocks, st)
        lowo = _laplace_mode_lowest(Mo, st.sub, n)
        # Nothing near a fold here and no other configuration: nothing to do.
        (isempty(blockset[b]) && (lowo === nothing || isempty(lowo.values))) && continue
        ho = _laplace_subtree_value(laplace, U, values, Ls, base, aws, st)
        isfinite(ho) || continue
        configs = Any[(y=base[st.idx], h=ho, M=Mo, F=nothing)]
        for y in blockset[b]
            length(y) == n || continue
            u = copy(base)
            u[st.idx] .= y
            h = _laplace_subtree_value(laplace, U, values, Ls, u, aws, st)
            isfinite(h) && push!(configs, (y=y, h=h, M=nothing, F=nothing))
        end
        for i in 1:length(configs)
            c = configs[i]
            u = copy(base)
            u[st.idx] .= c.y
            low = if i == 1
                lowo
            else
                Mc = c.M === nothing ?
                    _laplace_subtree_curvature(laplace, U, values, Ls, u, st) : c.M
                _laplace_mode_lowest(Mc, st.sub, n)
            end
            (low === nothing || isempty(low.values)) && continue
            v = low.vectors[:, 1]
            for s in (0.25, 0.5, 1.0, 2.0), sgn in (1.0, -1.0)
                r = _laplace_block_newton(laplace, U, values, Ls, aws, base, st,
                    c.y .+ (sgn * s) .* v)
                (r.converged && isfinite(r.value)) || continue
                _laplace_block_join!(configs, laplace, U, values, Ls, aws, base, st, r) &&
                    push!(blockset[b], r.y)
            end
        end
    end
    return blockset
end

"""
    _laplace_block_sums(laplace, U, theta, Ls, aws, zstar, tstar, Mstar)

Nested unit `U`'s blocks' combinations (`_laplace_block_sum`) at `theta`, the
unit at its maximum `zstar`, whose term is `tstar` and curvature `Mstar` (a
copy). Each block's stored other maxima (`block_modes[U]`) are followed from
where they stood, and those kept are written back. Returns the sum of the
blocks' `D_b` and what the gradient needs of the blocks that move it
(`blocks`), or `nothing` where the unit has none.
"""
function _laplace_block_sums(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, aws, zstar::Vector{Float64},
    tstar::Float64, Mstar::CTSEMBlockMatrix{Float64})
    entries = laplace.block_modes[U]
    isempty(entries) && return nothing
    blocks = laplace.units.blocks[U]
    # Where the unit's leaves integrate their own soft directions
    # (laplace_blocksoft.jl), each leaf configuration's value carries its own
    # correction, weighted as the unit's term weights them, `1 - w`.
    softw = nothing
    if _laplace_block_soft_on(laplace, U)
        rest = _laplace_rest_lowest(Mstar, blocks, _laplace_leaves(blocks))
        softw = 1 - (rest === nothing ? 0.0 :
            _laplace_smoothstep_weight(rest.lambda, laplace.soft_blend_lo,
                laplace.soft_blend_hi))
    end
    kept = CTSEMBlockMode[]
    parts = Any[]
    total = 0.0
    for b in sort!(unique(Int[e.b for e in entries]))
        (1 <= b <= length(blocks) && !isempty(blocks[b].ancestors)) || continue
        st = _laplace_subtree(blocks, b)
        n = length(st.idx)
        alts = Vector{Float64}[e.y for e in entries if e.b == b && length(e.y) == n]
        r = _laplace_block_sum(laplace, U, theta, Ls, aws, zstar, tstar, Mstar, st, alts;
            softw=(softw !== nothing && length(st.S) == 1) ? softw : nothing)
        if r === nothing
            # Not scored here (its own curvature does not factor): its maxima
            # are kept for the next evaluation, and it adds nothing.
            append!(kept, CTSEMBlockMode[CTSEMBlockMode(b, y) for y in alts])
            continue
        end
        append!(kept, CTSEMBlockMode[CTSEMBlockMode(b, y) for y in r.kept])
        total += r.value
        r.active && push!(parts, r)
    end
    laplace.block_modes[U] = kept
    return (value=total, blocks=parts)
end

"""
    _laplace_block_sum(laplace, U, theta, Ls, aws, zstar, tstar, Mstar, st, alts)

One block's `D_b` (the file's header): its other maxima `alts` followed by
`_laplace_block_newton` given `zstar`, joined to its own configuration
(`_laplace_block_join!`), scored, faded and combined. Returns the value, the
maxima kept, whether its gradient is not zero (`active`), and the seeds of that
gradient (`_laplace_block_sum_gradient!`); `nothing` where the own
configuration's curvature does not factor or its value is not finite.
"""
function _laplace_block_sum(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, aws, zstar::Vector{Float64},
    tstar::Float64, Mstar::CTSEMBlockMatrix{Float64}, st, alts::Vector{Vector{Float64}};
    softw::Union{Nothing,Float64}=nothing)
    blocks = laplace.units.blocks[U]
    idx = st.idx
    n = length(idx)
    Mo = _laplace_subtree_matrix(Mstar, blocks, st)
    Fo = _soft_shifted_factor(Mo, st.sub, 0.0)
    Fo === nothing && return nothing
    ho = _laplace_subtree_value(laplace, U, theta, Ls, zstar, aws, st)
    isfinite(ho) || return nothing
    tauo = ho - Fo.logdet / 2
    configs = Any[(y=zstar[idx], h=ho, M=Mo, F=Fo)]
    for y0 in alts
        r = _laplace_block_newton(laplace, U, theta, Ls, aws, zstar, st, y0)
        (r.converged && isfinite(r.value) && r.F !== nothing) || continue
        # A maximum more than the mode sum's span below the own one is out of
        # reach of any share (`_LAPLACE_MODE_SPAN`).
        r.value - r.F.logdet / 2 >= tauo - _LAPLACE_MODE_SPAN || continue
        _laplace_block_join!(configs, laplace, U, theta, Ls, aws, zstar, st, r)
    end
    K = length(configs)
    kept = Vector{Float64}[configs[c].y for c in 2:K]
    K == 1 && return (value=0.0, kept=kept, active=false)
    tau = Float64[c.h - c.F.logdet / 2 for c in configs]
    dtau = tau .- tau[1]
    # A leaf whose soft directions it integrates itself: each configuration's
    # correction (`_laplace_leaf_soft`), the own one's already in the unit's
    # term, so a configuration's value against it moves by the difference.
    Sv = zeros(Float64, K)
    if softw !== nothing
        point = y -> (u = copy(zstar); u[idx] .= y; u)
        for c in 1:K
            _laplace_leaf_soft_gate(laplace, configs[c].M, 1,
                laplace.soft_blend_hi) || continue
            sc = _laplace_leaf_soft(laplace, U, theta, Ls, aws, point(configs[c].y), st,
                configs[c].M)
            sc === nothing && return nothing
            Sv[c] = sc.S
        end
        dtau .+= softw .* (Sv .- Sv[1])
    end
    low = Any[_laplace_mode_lowest(c.M, st.sub, n) for c in configs]
    any(isnothing, low) && return nothing
    fade, dscale = _LAPLACE_MODE_FADE, _LAPLACE_MODE_DEPTH
    valley = fill(Inf, K, K)
    vt = zeros(Float64, K, K)
    for k in 1:K, j in (k + 1):K
        h = _laplace_block_segment(laplace, U, theta, Ls, aws, zstar, st,
            configs[k].y, configs[j].y)
        # The valley is the segment's lowest point, so where the middle lies
        # `_LAPLACE_MODE_DEPTH` below both ends neither fades by it, and the
        # middle stands in for it exactly: one pass rather than seven.
        mid = h(0.5)
        v = (isfinite(mid) && mid <= min(configs[k].h, configs[j].h) - dscale) ?
            (value=mid, t=0.5) : _laplace_segment_valley(h, configs[k].h, configs[j].h)
        valley[k, j] = valley[j, k] = v.value
        vt[k, j] = v.t
        vt[j, k] = 1 - v.t
    end
    x = Vector{Float64}[low[c].values ./ fade for c in 1:K]
    yv = Vector{Float64}[Float64[j == c ? Inf :
        max(configs[c].h - valley[c, j], 0.0) / dscale for j in 1:K] for c in 1:K]
    f = Float64[prod(_laplace_fade_weight, x[c]; init=1.0) *
        prod(_laplace_fade_weight, yv[c]; init=1.0) for c in 1:K]
    H = prod(f)
    q = Float64[c == 1 ? 0.0 : _laplace_fade_weight((abs(dtau[c]) - dscale) / dscale)
        for c in 1:K]
    s = Float64[c == 1 ? 0.0 : 1 - q[c] * H for c in 1:K]
    # The whole unit's value where it counts: `z_c`, the unit's maximum with
    # the block at configuration `c`, from there by the unit's Newton, kept
    # where it stays in that configuration's basin.
    J = copy(dtau)
    wholes = Any[nothing for _ in 1:K]
    for c in 2:K
        s[c] > 0 || continue
        u = copy(zstar)
        u[idx] .= configs[c].y
        rn = _laplace_newton_unit_mode(laplace, U, theta, Ls, aws, u)
        (rn.converged && isfinite(rn.value)) || continue
        _laplace_left_basin(laplace, U, theta, Ls, aws, u, rn) && continue
        zc = copy(rn.u)
        p = _laplace_unit_primal(laplace, U, theta, Ls, copy(zc), aws; given=rn.final)
        isfinite(p.term) || continue
        # The unit's curvature is the slot's scratch: kept, it is a copy.
        wholes[c] = (z=zc, term=p.term, M=_laplace_copy_blockmatrix(p.M), fac=p.fac,
            soft_quadrature=p.soft_quadrature, floored=p.floored)
        J[c] = p.term - tstar
    end
    a = Float64[c == 1 ? 0.0 : (wholes[c] === nothing ? dtau[c] :
        (1 - s[c]) * dtau[c] + s[c] * J[c]) for c in 1:K]
    shares = _laplace_mode_shares(f, a)
    pi = shares.pi
    value = sum(pi[c] * a[c] for c in 1:K if pi[c] > 0; init=0.0)
    # dD = sum_c pi_c da_c + sum_c C_c dlog f_c, the shares depending on the
    # values only through their ranking (`_laplace_mode_fade_coefficients`);
    # da_c = (1 - s_c) d dtau_c + s_c (dT(z_c) - dT(z*)) + (J_c - dtau_c) ds_c,
    # ds_c = -q_c H dlog H - H dq_c, and dlog H the sum of the fades' dlog f.
    C = _laplace_mode_fade_coefficients(shares, f, a, value)
    hcoef = 0.0
    ktau = zeros(Float64, K)
    kt = zeros(Float64, K)
    for c in 2:K
        pi[c] == 0 && continue
        if wholes[c] === nothing
            ktau[c] = pi[c]
            continue
        end
        e = (abs(dtau[c]) - dscale) / dscale
        dq = _laplace_fade_slope(e) * sign(dtau[c]) / dscale
        ktau[c] = pi[c] * ((1 - s[c]) - (J[c] - dtau[c]) * H * dq)
        kt[c] = pi[c] * s[c]
        hcoef -= pi[c] * (J[c] - dtau[c]) * q[c] * H
    end
    ktau[1] = -sum(ktau[2:K])
    # Each configuration's seeds: the coefficient on `h_b` there (`alpha`), on
    # `logdet(M_SS) / 2` (`gam`, as `-1/2 M^-1`'s scale), its fading
    # eigenvalues' `v v'` (`w`, `v`), and the valleys' points (`pts`).
    alpha = copy(ktau)
    gam = -ktau ./ 2
    lw = [Float64[] for _ in 1:K]
    lv = [Vector{Float64}[] for _ in 1:K]
    pts = Any[]
    for c in 1:K
        f[c] > 0 || continue
        cf = C[c] + hcoef
        cf == 0 && continue
        for i in eachindex(x[c])
            sx = _laplace_fade_slope(x[c][i])
            sx == 0 && continue
            push!(lw[c], cf * sx / (_laplace_fade_weight(x[c][i]) * fade))
            push!(lv[c], low[c].vectors[:, i])
        end
        for j in 1:K
            sy = _laplace_fade_slope(yv[c][j])
            sy == 0 && continue
            e = cf * sy / (_laplace_fade_weight(yv[c][j]) * dscale)
            alpha[c] += e
            t = vt[c, j]
            push!(pts, (e=e, c=c, j=j, t=t,
                y=configs[c].y .+ t .* (configs[j].y .- configs[c].y)))
        end
    end
    # The corrections' seeds: `d dtau_c` carries `softw (dS_c - dS_own)` and
    # `-(S_c - S_own) dw`, so each configuration's `dS_c` takes `softw
    # ktau_c` and the rest's weight `-sum_c ktau_c S_c`.
    kS = softw === nothing ? zeros(Float64, K) : softw .* ktau
    dwcoef = softw === nothing ? 0.0 : -sum(ktau[c] * Sv[c] for c in 1:K)
    active = any(!iszero, alpha) || any(!isempty, lw) || !isempty(pts) ||
        any(!iszero, kt) || any(c -> kS[c] != 0 && Sv[c] != 0, 1:K)
    return (value=value, kept=kept, active=active, st=st, configs=configs,
        alpha=alpha, gam=gam, lw=lw, lv=lv, pts=pts, kt=kt, kstar=-sum(kt),
        wholes=wholes, pi=pi, f=f, s=s, dtau=dtau, J=J, a=a, S=Sv, kS=kS,
        dwcoef=dwcoef)
end

"""
    _laplace_block_sum_gradient!(thetabar, zbar, laplace, U, theta, Ls, dL,
        positions, r)

One block's `dD_b` (`_laplace_block_sum`'s `r`) but for the whole unit's
routes: into `thetabar` with respect to the raw vector, and into `zbar` what it
leaves on the unit's coordinates at `z*`, which the caller carries back
through `z*`'s own movement. The valleys' points and each configuration's
seeds are passes over the subtree's members (`_laplace_soft_quadrature_vjp!`
on its objective, `_SoftQuadratureThetaGradient`'s `group` and `prior`); each
other configuration's movement is a block solve and one more pass. `false`
where a pass fails.
"""
function _laplace_block_sum_gradient!(thetabar::Vector{Float64}, zbar::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dL, positions, zstar::Vector{Float64}, r)
    st = r.st
    idx = st.idx
    d = length(zbar)
    G = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dL, positions, nothing;
        group=collect(Int, st.members), prior=collect(Int, idx))
    K = length(r.configs)
    point(y) = (u = copy(zstar); u[idx] .= y; u)
    zS = [zeros(Float64, length(idx)) for _ in 1:K]
    rest = zeros(Float64, d)
    for q in r.pts
        g = G(point(q.y))
        g === nothing && return false
        thetabar .-= q.e .* g.theta
        @inbounds for i in 1:d
            k = st.loc[i]
            if k == 0
                rest[i] -= q.e * g.z[i]
            else
                zS[q.c][k] -= q.e * (1 - q.t) * g.z[i]
                zS[q.j][k] -= q.e * q.t * g.z[i]
            end
        end
    end
    for c in 1:K
        (r.alpha[c] != 0 || r.gam[c] != 0 || !isempty(r.lw[c])) || continue
        Gamma = _SubtreeMemberGamma(_SoftGamma(r.gam[c], r.configs[c].F, st.sub,
            r.lw[c], r.lv[c], r.lv[c]), st.loc)
        zb = _laplace_soft_quadrature_vjp!(thetabar, G, point(r.configs[c].y),
            r.alpha[c], nothing, Gamma)
        zb === nothing && return false
        @inbounds for i in 1:d
            k = st.loc[i]
            k == 0 ? (rest[i] += zb[i]) : (zS[c][k] += zb[i])
        end
    end
    # The leaf's corrections, each through its configuration's own maximum
    # (`_laplace_leaf_soft_gradient!`), what they leave on the rest to `z*`.
    for c in 1:K
        (r.kS[c] != 0 && r.S[c] != 0) || continue
        zb = zeros(Float64, d)
        _laplace_leaf_soft_gradient!(thetabar, zb, laplace, U, theta, Ls, dL, positions,
            _laplace_workspace!(laplace, Float64, length(theta)), point(r.configs[c].y),
            st, r.configs[c].M, r.kS[c]) === nothing && return false
        rest .+= zb
    end
    # Each other configuration's own movement, given `z*`: `dy = M_SS^-1 d(grad
    # h_b)`, along its subtree's members.
    for c in 2:K
        any(!iszero, zS[c]) || continue
        beta = _soft_solve(r.configs[c].F, st.sub, zS[c])
        B = zeros(Float64, d)
        B[idx] .= beta
        g = _soft_quadrature_directional(G, point(r.configs[c].y), B)
        g === nothing && return false
        thetabar .+= g.theta
        @inbounds for i in 1:d
            st.loc[i] == 0 && (rest[i] += g.z[i])
        end
    end
    @inbounds for (k, i) in enumerate(idx)
        zbar[i] += zS[1][k]
    end
    zbar .+= rest
    return true
end

"""
    _laplace_unit_combo_gradient!(row, laplace, U, theta, Ls, dLlevels, dLsolved,
        combo)

Unit `U`'s gradient where its blocks combine their own maxima
(`_laplace_store_unit_primal!`'s `combo`: `top`, the whole-unit mode sum's or
`nothing`; `blocks`, `_laplace_block_sum`'s; `lead`, `z*` with its primal):
the unit's own gradient (the mode sum's, or its route at `z*`), the routes at
`z*` and each `z_c` the blocks weight, the blocks' own passes, and what those
leave on `z*`'s coordinates by its movement, `dz* = M*^-1 d(grad g)` (one solve
with the unit's curvature and one pass). `false` where a route or a pass fails.
"""
function _laplace_unit_combo_gradient!(row::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64}, Ls, dLlevels,
    dLsolved, combo)
    npar = length(theta)
    lead = combo.lead
    if combo.top !== nothing
        _laplace_mode_sum_gradient!(row, laplace, U, theta, Ls, dLlevels, dLsolved,
            combo.top) || return false
    else
        good, _ = _laplace_route_at!(row, laplace, U, theta, Ls, dLlevels, dLsolved,
            lead.z, lead)
        good || return false
    end
    part = zeros(Float64, npar)
    ks = sum(r.kstar for r in combo.blocks; init=0.0)
    if ks != 0
        good, _ = _laplace_route_at!(part, laplace, U, theta, Ls, dLlevels, dLsolved,
            lead.z, lead)
        good || return false
        row .+= ks .* part
    end
    for r in combo.blocks, c in eachindex(r.kt)
        r.kt[c] == 0 && continue
        w = r.wholes[c]
        good, _ = _laplace_route_at!(part, laplace, U, theta, Ls, dLlevels, dLsolved,
            w.z, w)
        good || return false
        row .+= r.kt[c] .* part
    end
    spec = laplace.spec
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    thetabar = zeros(Float64, npar)
    zbar = zeros(Float64, length(lead.z))
    for r in combo.blocks
        _laplace_block_sum_gradient!(thetabar, zbar, laplace, U, theta, Ls, dLlevels,
            positions, lead.z, r) || return false
    end
    # The leaves' corrections move with the rest's weight
    # (`_laplace_block_soft_term`), through its eigenvalue at `z*`.
    wcoef = sum(r.dwcoef for r in combo.blocks; init=0.0)
    if wcoef != 0
        blocks = laplace.units.blocks[U]
        rest = _laplace_rest_lowest(lead.M, blocks, _laplace_leaves(blocks))
        if rest !== nothing
            _, wp = _laplace_smoothstep_weight_slope(rest.lambda, laplace.soft_blend_lo,
                laplace.soft_blend_hi)
            if wp != 0
                glam = _laplace_rest_lambda_gradient(laplace, U, theta, Ls, dLlevels,
                    positions, lead.z, lead.M, rest.v)
                glam === nothing && return false
                thetabar .+= (wcoef * wp) .* glam
            end
        end
    end
    if any(!iszero, zbar)
        lead.F === nothing && return false
        blocks = laplace.units.blocks[U]
        beta = _soft_solve(lead.F, blocks, zbar)
        Gw = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dLlevels, positions,
            nothing)
        g = _soft_quadrature_directional(Gw, lead.z, beta)
        g === nothing && return false
        thetabar .+= g.theta
    end
    row .+= thetabar
    return all(isfinite, row)
end

"""
    _laplace_route_at!(row, laplace, U, theta, Ls, dLlevels, dLsolved, z, p)

The unit's gradient at its maximum `z` by its own route
(`_laplace_unit_gradient_route!`), `p` that maximum's primal (its curvature
`M`, factorization `fac`, and how its term was scored), the unit's records set
to that maximum's for the call and put back after it.
"""
function _laplace_route_at!(row::Vector{Float64}, laplace::CTSEMLaplaceObjective,
    U::Integer, theta::Vector{Float64}, Ls, dLlevels, dLsolved, z::Vector{Float64}, p)
    saved = (laplace.modes[U], laplace.logdet_floored[U], laplace.mode_repaired[U])
    laplace.modes[U] = z
    laplace.logdet_floored[U] = p.floored
    laplace.mode_repaired[U] = p.fac.repaired
    try
        return _laplace_unit_gradient_route!(row, laplace, U, theta, Ls, dLlevels,
            dLsolved, p.soft_quadrature, p.M, (p.fac.factors, p.fac.coupling))
    finally
        laplace.modes[U], laplace.logdet_floored[U], laplace.mode_repaired[U] = saved
    end
end
