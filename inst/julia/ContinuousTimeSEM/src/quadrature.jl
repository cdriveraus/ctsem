"""
Adaptive Gauss-Hermite quadrature over the same inner integral the Laplace
approximation approximates -- as a reference, as a diagnostic, and as a
correction to the estimate.

# Why this exists

`laplace.jl` replaces each subject's integral

    I_i = int exp(ll_i(theta, z)) phi(z) dz

by a single Gaussian fitted at the mode. That is exact when the integrand is
Gaussian in `z`, which it is when the varying parameters have identity
transforms and enter the state mean linearly, and it is *not* exact otherwise --
most visibly for a random effect on DRIFT, whose transform is `-log1p_exp`.

The error is not a constant. Profiling one 40-subject dataset along the
population sd of DRIFT, holding everything else at the Laplace estimate, the
gap between the Laplace value and a nine-node quadrature runs

    pop sd    0.36   0.50   0.67   0.89   1.14   1.43   1.75   2.09
    gap       0.89   1.63   2.60   3.62   4.60   5.62   6.66   7.64

-- monotone in the scale. A term that grows with the parameter tilts the
profile, and the tilt is the whole of the shrinkage: on that dataset the
Laplace profile peaks at a population sd of 0.67 and the quadrature profile at
0.89, against a generating 1.00. Two thirds of the shortfall is the
approximation rather than the estimand or the data.

# What is here

  * `ctsem_laplace_quadrature` -- the log marginal likelihood by adaptive
    Gauss-Hermite quadrature, using the mode and curvature the Laplace fit has
    already computed as the rule's location and scale -- the curvature
    clipped from below at the prior's, `M~ = V max(Lambda, 1) V'`, see
    `_quadrature_clipped_scale`. `nodes = 1` reproduces the Laplace value
    exactly wherever `M >= I`, which is what `test_quadrature.jl` asserts; on a
    unit with an eigenvalue below one it gives `g(uhat) - sum log max(lambda, 1)
    / 2` instead, the eigenwise-floored term, because that is the one-point rule
    at the clipped scale.
  * `ctsem_laplace_correction` -- the first-order correction to the estimate,
    `delta = (-H)^-1 grad(quadrature - laplace)`, at the cost of `2 * npar`
    quadrature evaluations and no refit. Cheap enough to report beside every
    fit; large entries say the estimate is meaningfully approximation-limited.
  * `ctsem_laplace_refine` -- the estimate that actually maximises the
    quadrature objective, for when the correction says it matters.

# How a hierarchy is handled

A unit spanning several subjects does not factor over them, so a product rule
over the whole unit latent vector would need `m^dim(u)` nodes and `dim(u)` grows
with the study. But the integral does factor *conditionally* -- given the study
effect the subjects are independent -- and that is exactly the tree the block
structure already encodes. So the rule recurses over it:

    value(b, u) = logsumexp_j [ w_j - z_j'z_j/2
                                + ( sum over b's children of value(child, u) ,
                                    or b's own members' log likelihoods when b
                                    is a leaf ) ]

with `u[b] = z_j` set before descending. Cost is `m^k` per block times the
number of blocks, so it is *linear* in the number of groups at each level and
exponential only in that level's own number of random effects -- 5 nodes on a
one-effect study level over 40 subjects with two subject effects is
`5 * 40 * 25` process evaluations, not `5^41`.

Two different rules, for two different reasons. A leaf block re-solves its
conditional mode at every outer node, because that mode moves a long way with
the effect above it and the accuracy of the whole thing rests on it. An outer
block keeps the joint mode, and takes its scale from the *eliminated* diagonal
the block factorization already produced: the curvature of that block after
every block beneath it has been integrated out, which is the conditional
precision of the block given its ancestors under the Gaussian approximation --
exactly the distribution the recursion is standing in when it gets there. The
rule does not need the exact conditional mode to be *valid*, only to be
efficient, and that centre and scale have already been paid for.

With one level the tree is a single leaf with no ancestors, the recursion is one
call, and the result is identical to the flat rule this file started as -- which
is what `test_quadrature.jl` pins by checking `nodes = 1` against the Laplace
value at one, two and three levels (on concave fixtures; see the clipping
above). With the clip, `nodes = 1` at more than one level clips each block's
eliminated precision separately, which is not the same as clipping the unit's
eigenvalues: the two agree only where nothing is clipped.
"""

using LinearAlgebra

export ctsem_laplace_quadrature, ctsem_laplace_correction, ctsem_laplace_refine

"""
    _gauss_hermite(m)

Nodes and weights of the `m`-point Gauss-Hermite rule for weight `exp(-x^2)`,
by Golub-Welsch: the nodes are the eigenvalues of the Jacobi matrix and the
weights come from the first components of its eigenvectors.

Computed rather than tabulated so that any `m` is available, and cached because
the same `m` is asked for once per subject per evaluation.
"""
function _gauss_hermite(m::Integer)
    m >= 1 || throw(ArgumentError("need at least one quadrature node"))
    key = Int(m)
    # Locked, following `_TRANSFORM_CACHE` in `r_interface.jl`. Only the
    # quadrature path warms this cache first; the binary measurement kernels
    # (`binary_measurement.jl`, `kalman_filters.jl`, `kalman_trace.jl`) call it
    # once per row from every threaded loop in the engine, so an ordinary fit
    # with binary indicators can reach a cold cache from several threads at
    # once, and a concurrent `setindex!` during a rehash corrupts a `Dict`. A
    # node count is built once per session and an uncontended lock is tens of
    # nanoseconds.
    lock(_GH_CACHE_LOCK) do
        cached = get(_GH_CACHE, key, nothing)
        cached === nothing || return cached
        if key == 1
            built = ([0.0], [sqrt(pi)])
        else
            offdiag = [sqrt(i / 2) for i in 1:(key - 1)]
            decomposition = eigen(SymTridiagonal(zeros(key), offdiag))
            built = (decomposition.values,
                sqrt(pi) .* (decomposition.vectors[1, :] .^ 2))
        end
        _GH_CACHE[key] = built
        return built
    end
end

const _GH_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const _GH_CACHE_LOCK = ReentrantLock()

"""
    _gh_grid(k, m)

The `k`-dimensional product rule as `(points, logweights)`, with the
`exp(x^2)` factor that turns the Gauss-Hermite weight into the one an adaptive
rule wants already folded into the log weight.
"""
function _gh_grid(k::Integer, m::Integer)
    key = (Int(k), Int(m))
    cached = get(_GH_GRID_CACHE, key, nothing)
    cached === nothing || return cached
    x, w = _gauss_hermite(m)
    npoints = Int(m)^Int(k)
    points = Vector{Vector{Float64}}(undef, npoints)
    logweights = Vector{Float64}(undef, npoints)
    for (position, index) in enumerate(Iterators.product(ntuple(_ -> 1:Int(m), Int(k))...))
        points[position] = [x[i] for i in index]
        logweights[position] = sum(log(w[i]) + x[i]^2 for i in index; init=0.0)
    end
    built = (points, logweights)
    _GH_GRID_CACHE[key] = built
    return built
end

# Unlocked, unlike `_GH_CACHE`: the only reader is `_quadrature_block`, which
# runs strictly after `_quadrature_warm_caches` has asked for every
# `(block.size, nodes)` key it will use, so the cache is read-only for the whole
# threaded region. Any new caller reached from a thread needs the warm too, or
# this needs a lock.
const _GH_GRID_CACHE = Dict{Tuple{Int,Int},Tuple{Vector{Vector{Float64}},Vector{Float64}}}()

"""
    _quadrature_warm_caches(laplace, nodes)

Populate the Gauss-Hermite caches for every size this evaluation will ask for.

Called from single-threaded code before the chunk loop spawns, so that
`_gauss_hermite` and `_gh_grid` are pure reads inside the threads. The set is
small and known in advance: one node count for the whole call, and one block
dimension per distinct block size in the unit tree.
"""
function _quadrature_warm_caches(laplace::CTSEMLaplaceObjective, nodes::Integer)
    m = Int(nodes)
    _gauss_hermite(m)
    # `_CTSEM_BINARY_NODES` is what the measurement kernels ask for, and it is
    # a different count from the quadrature's own.
    _gauss_hermite(_CTSEM_BINARY_NODES[])
    seen = Set{Int}()
    for blocks in laplace.units.blocks, block in blocks
        k = block.size
        k >= 1 && !(k in seen) && (push!(seen, k); _gh_grid(k, m))
    end
    return nothing
end

"""
    _quadrature_children(blocks)

For each block, the blocks immediately beneath it, and the blocks with nothing
above them.

`ancestors` is ordered outward from the immediate parent, so `c` sits directly
under `b` exactly when `c.ancestors[1] == b`. A block with no ancestors is a
root, and a unit normally has exactly one, since the unit *is* its outermost
group; several appear only when the outermost level carries no random effects
and the tree is a forest of the level below.
"""
function _quadrature_children(blocks::Vector{CTSEMLaplaceBlock})
    children = [Int[] for _ in blocks]
    roots = Int[]
    for (b, block) in enumerate(blocks)
        if isempty(block.ancestors)
            push!(roots, b)
        else
            push!(children[block.ancestors[1]], b)
        end
    end
    return (children=children, roots=roots)
end

"""
    _quadrature_clipped_scale(F, k)

The rule's scale and log scale from a block's precision `D = F.L F.L'`, with
`D`'s eigenvalues clipped from below at one: `S S' = inv(D~)`,
`D~ = V max(Lambda, 1) V'`.

The precision is the Gaussian the rule is centred on, and scaling the nodes by
`D^(-1/2)` places them where that Gaussian puts its mass. Where the log
likelihood is concave in `u` every eigenvalue is at least one and that is right.
Where it has gone convex it is not: an eigenvalue of 1e-3 scatters the nodes
some thirty prior standard deviations out, where the integrand is negligible,
and the rule reports mostly its own extrapolation. Measured on 40-subject
weak-data fits, a 5-point rule scaled by `D` was 1.8 nats high on average on
such units (2.9 at worst) against an exact reference, and 0.06 when clipped.
The clip never makes a node wider than the prior's own spread, and the rule is
a change of variables whatever `S` is, so it stays a quadrature of the same
integral; only where its nodes go changes.

`D - I` positive definite -- every eigenvalue above one -- returns the Cholesky
scale exactly as before, bit for bit. Otherwise `D~` is built by the engine's
own symmetric eigensolver and factored the same way, so the scale is continuous
in `D` across the switch (`D~ -> D` as its smallest eigenvalue rises to one).
"""
function _quadrature_clipped_scale(F::CTSEMCholesky, k::Integer)
    Lf = Matrix(F.L)
    D = Lf * transpose(Lf)
    shifted = D - Matrix{Float64}(LinearAlgebra.I, k, k)
    issuccess(_ctsem_cholesky(shifted, k)) &&
        return (scale=_ctsem_cholesky_uinv(F), logdetscale=-logdet(F) / 2,
            clipped=false)
    E = _ctsem_symeig(D)
    clipped = E.vectors * Diagonal(max.(E.values, 1.0)) * transpose(E.vectors)
    G = _ctsem_cholesky(Matrix(_laplace_symmetrise(clipped)), k)
    return (scale=_ctsem_cholesky_uinv(G), logdetscale=-logdet(G) / 2,
        clipped=true)
end

"""
    _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws)

The conditional mode and scale for a leaf block, with its ancestors held at
whatever `u` currently says.

The leaf mode is where the accuracy is. It moves a long way with the effect
above it -- a subject whose study effect has been pushed two standard deviations
has a different best guess about its own -- so re-solving it at every outer node
is the point of the recursion rather than an optimisation within it.

Newton on this block's coordinates only. The ancestors are data here, and a
member outside this block contributes nothing to its curvature, so the problem
is `k x k` however large the unit is.

Returns `(ok, centre, scale, logdetscale)` with `scale * scale' = M~^-1`, the
curvature clipped at the prior's; see `_quadrature_clipped_scale`.
"""
function _quadrature_leaf_rule!(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, aws)
    block = laplace.units.blocks[U][b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    members = block.members
    failure = (ok=false, centre=Float64[], scale=zeros(Float64, 0, 0), logdetscale=0.0)

    inner_gradient = function (z)
        work = copy(u)
        @inbounds for (t, c) in enumerate(columns); work[c] = z[t]; end
        result = _laplace_unit_loglik_gradient(laplace, U, theta, Ls, work, aws, members)
        isfinite(result.value) || return fill(NaN, k)
        return [result.gradient[c] for c in columns] .- z
    end
    # The log likelihood's gradient in this block alone, as a function of this
    # block alone -- differentiating it gives the block's curvature.
    loglik_gradient = function (z)
        S = eltype(z)
        ws = _laplace_workspace!(laplace, S, length(theta))
        work = convert(Vector{S}, u)
        @inbounds for (t, c) in enumerate(columns); work[c] = z[t]; end
        result = _laplace_unit_loglik_gradient(laplace, U, convert(Vector{S}, theta),
            [convert(Matrix{S}, L) for L in Ls], work, ws, members)
        return [result.gradient[c] for c in columns]
    end
    curvature_at = function (z)
        A = ForwardDiff.jacobian(loglik_gradient, z)
        M = Matrix{Float64}(LinearAlgebra.I, k, k) .- _laplace_symmetrise(A)
        # The engine's factorization rather than LAPACK's: once per inner Newton
        # step per subject, inside the subject loop. See small_linalg.jl.
        return _ctsem_cholesky(M, k)
    end

    z = Float64[u[c] for c in columns]
    gradient = inner_gradient(z)
    all(isfinite, gradient) || return failure
    for _ in 1:laplace.inner_maxiter
        maximum(abs, gradient) < laplace.inner_tol && break
        factorization = curvature_at(z)
        issuccess(factorization) || return failure
        candidate = z .+ (factorization \ gradient)
        trial = inner_gradient(candidate)
        all(isfinite, trial) || break
        z = candidate
        gradient = trial
    end
    factorization = curvature_at(z)
    issuccess(factorization) || return failure
    rule = _quadrature_clipped_scale(factorization, k)
    return (ok=true, centre=z, scale=rule.scale, logdetscale=rule.logdetscale)
end

"""
    _quadrature_block(laplace, U, theta, Ls, b, u, context, aws)

Block `b`'s contribution to its unit's log marginal: its own coordinates
integrated, and whatever sits beneath it recursed into.

`u` carries the values already chosen for `b`'s ancestors and is written in
place as the recursion descends. Entries belonging to other branches are stale
and unread, because a member reads only the blocks on its own path.
"""
function _quadrature_block(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, context, aws)
    blocks = laplace.units.blocks[U]
    block = blocks[b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    children = context.children[b]
    leaf = isempty(children)

    rule = if leaf
        _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws)
    else
        # An outer block keeps the joint mode, and takes its scale from the
        # *eliminated* diagonal the block factorization already produced: the
        # curvature of this block after every block beneath it has been
        # integrated out. That is the conditional precision of `u_b` given its
        # ancestors under the Gaussian approximation, which is exactly the
        # distribution the recursion is standing in when it gets here.
        #
        # The marginal covariance -- this block's diagonal of `inv(M)` -- is the
        # wrong object and only looks right at two levels, where a block with
        # nothing above it has marginal and conditional coincide. At three it
        # loses the coupling to the level above, and the `nodes = 1` identity
        # (which forces `sum_b logdet(scale_b) = -logdet(M)/2`) fails, because
        # the elimination is what makes that sum telescope.
        # Clipped the same way as a leaf's; the eliminated diagonal of an
        # `M >= I` is itself at least `I` (its inverse is a diagonal block of
        # `inv(M) <= I`), so this changes nothing where the unit is concave.
        clippedrule = _quadrature_clipped_scale(context.factors[b], k)
        (ok=true, centre=Float64[context.mode[c] for c in columns],
         scale=clippedrule.scale, logdetscale=clippedrule.logdetscale)
    end
    rule.ok || return NaN

    points, logweights = _gh_grid(k, context.nodes)
    terms = Vector{Float64}(undef, length(points))
    @inbounds for j in eachindex(points)
        z = rule.centre .+ sqrt(2) .* (rule.scale * points[j])
        for (t, c) in enumerate(columns); u[c] = z[t]; end
        inner = 0.0
        if leaf
            for m in block.members
                shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                    laplace.units.offsets[U][m])
                inner += laplace.objective.subject_objectives[
                    laplace.units.members[U][m]](shifted)
            end
        else
            for c in children
                inner += _quadrature_block(laplace, U, theta, Ls, c, u, context,
                    aws)
            end
        end
        terms[j] = isfinite(inner) ? inner - dot(z, z) / 2 + logweights[j] : -Inf
    end
    peak = maximum(terms)
    isfinite(peak) || return NaN
    accumulated = 0.0
    @inbounds for j in eachindex(terms)
        accumulated += exp(terms[j] - peak)
    end
    # `(2 pi)^(-k/2)` from the standard normal density and `2^(k/2)` from the
    # change of variable `z = centre + sqrt(2) * scale * x`. The Laplace value's
    # own constants cancel against each other instead -- see `_laplace_unit_term`
    # -- so they are written out here rather than shared.
    return peak + log(accumulated) + rule.logdetscale + k * (log(2) - log(2 * pi)) / 2
end

"""
    ctsem_laplace_quadrature(laplace, values; nodes=5, contributions=false)

The log marginal likelihood by adaptive Gauss-Hermite quadrature, plus the
prior -- the same quantity `ctsem_laplace_evaluate` approximates, computed with
`nodes` points per random effect instead of one.

The rule is *adaptive* in the standard sense: centred at an inner mode and
scaled by the inverse of the curvature there, both of which the Laplace
machinery already produces -- with the curvature's eigenvalues clipped from
below at one, the prior's own, so that no node is placed wider than the prior
spreads (`_quadrature_clipped_scale`). Where every eigenvalue is at least one
that is the unclipped curvature, bit for bit, and with `nodes = 1` the rule has
a single point at the mode and the result is the Laplace value to machine
precision, the cheapest available check that the two agree about what they are
integrating. Where one is below one, `nodes = 1` gives the eigenwise-floored
term `g(uhat) - sum log max(lambda, 1) / 2` rather than the fit's own, so the
gap `ctLaplaceCheck` reports is no longer zero at one node on such a unit.

Cost is `nodes^k` process log likelihoods per block. At one level that is
`nodes^k` per subject -- 5 or 25 for the common one or two random effects. With
a hierarchy it recurses over the block tree rather than taking a product over
it, so the cost stays linear in the number of groups at each level; see the
module docstring.

`contributions` returns the per-unit terms, which sum to the value less the
prior. They are per *unit*, not per subject: with a group level the integral
does not decompose over a group's members, which is the same reason
`ctsem_laplace_evaluate` spreads its own per-unit term evenly.
"""
function ctsem_laplace_quadrature(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, contributions::Bool=false)
    nodes >= 1 || throw(ArgumentError("need at least one quadrature node"))
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    nunits = length(laplace.units.members)
    unit_term = zeros(Float64, nunits)
    Ls = _laplace_popchols(theta, laplace.spec)

    nchunks = _ctsem_nchunks(nunits)
    # The pool sizes the scratch store, not this function, and its workers carry
    # their own slot -- so the unit functions below need nothing passed to them,
    # and a nested region gets a sub-budget rather than the whole pool again.
    _laplace_ensure_pool!(laplace)
    ranges = _ctsem_chunk_assignment(_laplace_unit_weights(laplace), nchunks)
    failed = fill(false, nchunks)

    # Fill the quadrature caches here, before anything is spawned.
    #
    # `_gauss_hermite` and `_gh_grid` memoise into plain global `Dict`s, and the
    # chunks below call both -- `_gh_grid` once per block, `_gauss_hermite` once
    # per row through the measurement kernels. Concurrent `setindex!` during a
    # rehash corrupts a `Dict`; `laplace.jl` already carries the note about
    # finding exactly that on the workspace store, where Julia caught it with
    # "Multiple concurrent writes to Dict detected!" rather than returning wrong
    # numbers. There is nothing here that makes these two safer, only rarer:
    # they are hit by every thread on the first row it touches.
    #
    # `_gauss_hermite` now takes a lock as well, because its other callers are
    # on threaded paths that never reach this warm. Warming still matters here:
    # after this loop `_GH_GRID_CACHE` is read-only for the rest of the call,
    # which is what lets it stay unlocked on the per-block path.
    _quadrature_warm_caches(laplace, nodes)

    # A trial point can be invalid in ways that *throw* rather than return a
    # non-finite number -- a parameter vector an optimizer wandered into can
    # make the matrix exponential's own scaling step take `ceil(Int, NaN)`. The
    # engine's contract is that an evaluation reports NaN, and a throw inside a
    # spawned task escapes as a `TaskFailedException` that kills whatever loop
    # is above it, so it is caught here rather than left to every caller.

    run = function (c)
        try
            _quadrature_chunk!(laplace, theta, Ls, ranges[c], unit_term, nodes)
        catch err
            # Numerical failures become NaN; bugs do not.
            #
            # This catch is here for a trial point the model cannot evaluate --
            # a parameter vector whose matrix exponential scaling step takes
            # `ceil(Int, NaN)` and throws rather than returning a non-finite
            # number. Converting *that* to NaN is the engine's contract.
            #
            # It caught a `MethodError` too, and that is a different thing
            # entirely. When the worker pool removed the `slot` argument from
            # the unit functions, two call sites here kept passing one; every
            # chunk raised `MethodError`, every chunk was marked failed, and
            # `ctLaplaceCheck` reported a NaN gap that read exactly like a
            # quadrature that could not be computed for this model. Nothing
            # errored, nothing warned, and the only reason it was found is that
            # a test asserted the result was finite.
            #
            # So the errors that mean "this code is wrong" are rethrown. A
            # `TaskFailedException` from a nested region is unwrapped first,
            # because the pool may have spawned inside the chunk.
            _ctsem_must_propagate(err) && rethrow()
            failed[c] = true
        end
        return nothing
    end
    _laplace_parallel(laplace, 1:nchunks) do c
        run(c)
        return true
    end
    if any(failed)
        return contributions ? (value=NaN, subject_loglik=unit_term) :
            (value=NaN, subject_loglik=Float64[])
    end
    total = sum(unit_term) + _ctsem_log_prior(laplace.objective, theta)
    return contributions ? (value=total, subject_loglik=unit_term) :
        (value=total, subject_loglik=Float64[])
end

"""
    _quadrature_chunk!(laplace, theta, Ls, units, unit_term, nodes)

One chunk of the unit loop, writing each unit's term into `unit_term`.

Throws rather than flagging: the caller wraps it, so that a failure inside a
spawned task becomes a NaN evaluation rather than an escaping exception.
"""
function _quadrature_chunk!(laplace::CTSEMLaplaceObjective, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, units, unit_term::Vector{Float64},
    nodes::Integer)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    for U in units
        blocks = laplace.units.blocks[U]
        if isempty(blocks)
            # No random effects anywhere in this unit: there is no integral,
            # and the term is the plain log likelihood.
            unit_term[U] = _laplace_unit_objective_gradient(laplace, U, theta,
                Ls, Float64[], aws).value
            continue
        end
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        mode = copy(laplace.modes[U])
        M = _laplace_unit_curvature(laplace, U, theta, Ls, mode)
        factorization = _laplace_factor_repaired!(M, blocks)
        factorization.ok ||
            error("the inner curvature could not be factorized")
        tree = _quadrature_children(blocks)
        context = (children=tree.children, nodes=Int(nodes), mode=mode,
            factors=factorization.factors)
        u = copy(mode)
        total = 0.0
        for root in tree.roots
            total += _quadrature_block(laplace, U, theta, Ls, root, u, context,
                aws)
        end
        isfinite(total) || error("quadrature produced a non-finite unit term")
        unit_term[U] = total
    end
    return nothing
end

"""
    _quadrature_gap_gradient(laplace, values, nodes, step)

Central difference of `quadrature - laplace` in every coordinate.

The *difference* is differenced, not the two objectives separately, because the
difference is the small, smooth quantity: both objectives carry the same
several-hundred-log-unit level and the same sharp curvature, and subtracting
first keeps the cancellation out of the finite difference.
"""
function _quadrature_gap_gradient(laplace::CTSEMLaplaceObjective,
    values::Vector{Float64}, nodes::Integer, step::Real)
    npar = length(values)
    gradient = zeros(Float64, npar)
    gap = function (x)
        result = try
            quad = ctsem_laplace_quadrature(laplace, x; nodes=nodes).value
            lap = ctsem_laplace_evaluate(laplace, x; gradient=false).value
            quad - lap
        catch err
            err isa InterruptException && rethrow()
            NaN
        end
        # A step that lands outside the support carries no information about the
        # gap; treating it as zero difference is the only honest fallback, and
        # it leaves the corresponding entry of the correction at zero rather
        # than at NaN, which would poison the whole Newton step.
        return isfinite(result) ? result : 0.0
    end
    for j in 1:npar
        h = step * max(1.0, abs(values[j]))
        plus = copy(values); plus[j] += h
        minus = copy(values); minus[j] -= h
        gradient[j] = (gap(plus) - gap(minus)) / (2h)
    end
    return gradient
end

"""
    ctsem_laplace_correction(laplace, values, hessian; nodes=5, step=1e-3)

How far the Laplace estimate sits from the one the exact integral would give,
to first order, without refitting.

At the Laplace estimate the Laplace gradient is zero, so the quadrature
objective's gradient there is exactly the *gap* gradient `grad(Q - T)`, and one
Newton step against the curvature already computed for the standard errors
gives

    delta = (-H)^-1 grad(Q - T)

Costs `2 * npar` quadrature evaluations. Report it beside the estimate: an
entry small against that parameter's standard error says the approximation is
not what is limiting the answer, and a large one says it is. On a 40-subject
model with a random DRIFT the population scale's entry is about two thirds of
the distance to the generating value.

`hessian` is the outer Hessian at `values`, as `ctsem_laplace_hessian` returns
it; it is taken as an argument rather than recomputed because the caller
computing standard errors has already paid for it.
"""
function ctsem_laplace_correction(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, hessian::AbstractMatrix; nodes::Integer=5,
    step::Real=1e-3)
    theta = collect(Float64, values)
    gap = _quadrature_gap_gradient(laplace, theta, nodes, step)
    H = Symmetric((hessian .+ transpose(hessian)) ./ 2)
    delta, dropped = _correction_step(H, gap)
    quadrature = ctsem_laplace_quadrature(laplace, theta; nodes=nodes).value
    laplacevalue = ctsem_laplace_evaluate(laplace, theta; gradient=false).value
    return (delta=delta, corrected=theta .+ delta, gap_gradient=gap,
        quadrature=quadrature, laplace=laplacevalue,
        gap=quadrature - laplacevalue, nodes=Int(nodes),
        dropped_directions=dropped)
end

"""
    _correction_step(H, gap) -> (delta, dropped)

Solve `-H \\ gap` along the directions `H` actually identifies, and report how
many it did not.

A plain solve here is wrong in a way that produces a number rather than an
error. The gap gradient is a central difference of the quadrature value, so on a
model where Laplace is exact it sits at floating-point noise -- measured at
7e-12 and 2e-11 on two runs of the same three-level fixture. If `H` is also
near-singular, and it is whenever a population scale is weakly identified (six
studies, in that fixture), `H \\ gap` multiplies that noise by 1/lambda_min.
The two runs returned corrections of 3.3e+08 and 3.5e-10 standard errors from
the same fit: one absurd, one fine, both meaningless, and nothing distinguished
them.

Truncating at `sqrt(eps)` relative to the largest eigenvalue is the standard
answer and the right one here: a first-order correction along a direction the
data does not identify is not estimable, and reporting zero for it is honest
where reporting 3e+08 is not. `dropped_directions` on the result says when it
happened rather than leaving it to be inferred.
"""
function _correction_step(H::Symmetric, gap::AbstractVector)
    n = length(gap)
    n == 0 && return (Float64[], 0)
    decomposition = try
        eigen(H)
    catch err
        err isa InterruptException && rethrow()
        return (fill(NaN, n), 0)
    end
    lambda = decomposition.values
    scale = maximum(abs, lambda)
    (!isfinite(scale) || scale == 0) && return (zeros(n), n)
    tolerance = sqrt(eps(Float64)) * scale
    projected = transpose(decomposition.vectors) * gap
    dropped = 0
    for i in eachindex(lambda)
        if abs(lambda[i]) <= tolerance
            projected[i] = 0.0
            dropped += 1
        else
            projected[i] /= lambda[i]
        end
    end
    return (-(decomposition.vectors * projected), dropped)
end

"""
    ctsem_laplace_refine(laplace, values; nodes=5, maxiter=50, step=1e-4)

Maximise the quadrature objective, started from the Laplace estimate.

For when `ctsem_laplace_correction` says the first-order step is large enough to
matter and a linear correction is not enough. The gradient is a central
difference of the quadrature value -- exact-to-quadrature but `2 * npar`
evaluations each, so this is minutes rather than the seconds a Laplace fit
takes, and it is deliberately not the default.

Starting at the Laplace estimate is what makes it affordable: the two optima
are close, the inner modes are warm, and L-BFGS has a good initial metric.
"""
function ctsem_laplace_refine(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, maxiter::Integer=50,
    step::Real=1e-4, g_tol::Real=1e-5, verbose::Bool=false)
    start = collect(Float64, values)
    invalid = floatmax(Float64) / 1e8
    objective = function (x)
        value = try
            ctsem_laplace_quadrature(laplace, x; nodes=nodes).value
        catch err
            err isa InterruptException && rethrow()
            NaN
        end
        return isfinite(value) ? -value : invalid
    end
    gradient! = function (G, x)
        for j in eachindex(x)
            h = step * max(1.0, abs(x[j]))
            plus = copy(x); plus[j] += h
            minus = copy(x); minus[j] -= h
            G[j] = (objective(plus) - objective(minus)) / (2h)
        end
        return G
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol,
        show_trace=verbose, store_trace=false)
    result = Optim.optimize(objective, gradient!, start, Optim.LBFGS(), options)
    minimizer = collect(Optim.minimizer(result))
    return (minimizer=minimizer,
        maximum_loglik=-Optim.minimum(result),
        laplace_start=start,
        shift=minimizer .- start,
        iterations=Optim.iterations(result),
        converged=Optim.converged(result),
        nodes=Int(nodes))
end

################################################################################
# Reference implementation of the gated floor's value
################################################################################
#
# The `:gated` floor (`_laplace_gated_term` in laplace.jl) is generic in its
# element type because its gradient differentiates it. These two functions
# compute the same rule in plain Float64, written separately and first, and
# `test_laplace.jl` holds the floor's value to them. They are not on any
# fitting path.

"""
    _laplace_soft_rule_unit(laplace, values, U; nodes=5, ndirs=1, tau=0.0,
                            recenter=true, newton_steps=0)

One unit's log marginal likelihood by a Gauss-Hermite rule along the unit's
softest directions and Laplace in the rest. Value only: the reference the
gated floor is tested against; see
`CT-SEM/review/LAPLACE-eigenwise-floor-2026-09-23.md`, second and third
addenda. Reads the modes of the last `ctsem_laplace_evaluate` at `values`.

In the eigenbasis `M = V Lambda V'` at the mode:

  * the `ndirs` smallest directions, and any other with `lambda < tau`, get
    `nodes` points at the *clipped* scale `1/sqrt(max(lambda, 1))` -- never
    wider than the prior, which is what `ctsem_laplace_quadrature`'s
    `M^(-1/2)` scaling gets wrong at a near-singular unit;
  * at each such node the remaining coordinates are moved to their conditional
    mode (`recenter`, Newton) and integrated by Laplace with their conditional
    curvature, clipped the same way. Without the recentring the rule misses the
    ridge a nonlinear effect bends the integrand along, and does worse than one
    node; without the clip, a conditional curvature that collapses at a far node
    produces a spike of several nats.

`newton_steps = k > 0` replaces the converged recentring by `k` Newton steps,
each with the node's own conditional curvature, and takes the conditional
logdet from the last curvature evaluated: `k` gradients, `k` curvatures and one
value per node. `k = 1` at 3 nodes is the cheap rule of the third addendum.
The middle node of an odd rule sits at the mode and costs nothing.

`nodes = 1`, `ndirs = 0`, `tau = 1` is `g(uhat) - sum log max(lambda, 1) / 2`.
Uses `_ctsem_symeig` and `_ctsem_cholesky` throughout (no LAPACK); a rule of
`n` nodes costs about `n` conditional Newton solves of a few steps, each a
curvature of the unit, so `n` times the primal work of one unit.
"""
function _laplace_soft_rule_unit(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, U::Integer; nodes::Integer=5, ndirs::Integer=1,
    tau::Real=0.0, recenter::Bool=true, newton_steps::Integer=0)
    theta = collect(Float64, values)
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    uhat = copy(laplace.modes[U])
    d = length(uhat)
    blocks = laplace.units.blocks[U]
    gof(u) = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
    Mof(u) = _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta, Ls, u),
        blocks, d)
    d == 0 && return gof(uhat).value
    E = _ctsem_symeig(Mof(uhat))
    soft = sort(unique(vcat(collect(1:min(Int(ndirs), d)),
        findall(<(tau), E.values))))
    stiff = setdiff(1:d, soft)
    ns, nh = length(soft), length(stiff)
    Vs, Vh = E.vectors[:, soft], E.vectors[:, stiff]
    lsoft = max.(E.values[soft], 1.0)
    xs, ws = _gauss_hermite(nodes)

    # log int exp(g(ubase + Vh z)) dz by Laplace at the conditional mode, with
    # the conditional curvature clipped at the prior's.
    stiff_part = function (ubase::Vector{Float64})
        local z, r, gz, Mz, F, ez, lam, uu
        nh == 0 && return gof(ubase).value
        z = zeros(nh)
        if recenter && newton_steps > 0
            lam = max.(E.values[stiff], 1.0)
            maximum(abs, ubase .- uhat) < 1e-12 &&
                return gof(uhat).value - sum(log, lam; init=0.0) / 2
            for _ in 1:Int(newton_steps)
                r = gof(ubase .+ Vh * z)
                ez = _ctsem_symeig(transpose(Vh) * Mof(ubase .+ Vh * z) * Vh)
                lam = max.(ez.values, 1.0)
                z .+= ez.vectors * ((transpose(ez.vectors) *
                    (transpose(Vh) * r.gradient)) ./ lam)
            end
            return gof(ubase .+ Vh * z).value - sum(log, lam; init=0.0) / 2
        end
        if recenter
            for _ in 1:50
                r = gof(ubase .+ Vh * z)
                gz = transpose(Vh) * r.gradient
                maximum(abs, gz) < laplace.inner_tol && break
                Mz = transpose(Vh) * Mof(ubase .+ Vh * z) * Vh
                F = _ctsem_cholesky(Matrix(_laplace_symmetrise(Mz)), nh)
                issuccess(F) || break
                z .+= F \ gz
            end
        end
        uu = ubase .+ Vh * z
        ez = _ctsem_symeig(transpose(Vh) * Mof(uu) * Vh)
        lam = recenter ? max.(ez.values, 1.0) : E.values[stiff]
        return gof(uu).value - sum(log, lam; init=0.0) / 2
    end
    terms = Float64[]
    for idx in Iterators.product(ntuple(_ -> 1:Int(nodes), ns)...)
        x = Float64[xs[i] for i in idx]
        v = stiff_part(uhat .+ Vs * (sqrt(2) .* x ./ sqrt.(lsoft)))
        push!(terms, isfinite(v) ?
            v + sum(log(ws[i]) + xs[i]^2 for i in idx; init=0.0) : -Inf)
    end
    peak = maximum(terms)
    isfinite(peak) || return NaN
    # Per soft direction sqrt(2/lambda~) from the change of variable against the
    # sqrt(pi) the weights carry and the (2 pi)^(-1/2) of the density.
    return peak + log(sum(exp.(terms .- peak))) - sum(log, lsoft; init=0.0) / 2 -
        ns * log(pi) / 2
end


"""
    _laplace_soft_weight(lambda, lo, hi)

The C1 hand-off from the soft rule to the Laplace term: one below `lo`, zero
above `hi`, and `1 - (3 s^2 - 2 s^3)` between, with `s = (lambda - lo)/(hi - lo)`.
"""
function _laplace_soft_weight(lambda::Real, lo::Real, hi::Real)
    s = clamp((lambda - lo) / (hi - lo), 0.0, 1.0)
    return 1 - (3 * s^2 - 2 * s^3)
end

"""
    _laplace_gated_unit_reference(laplace, values, U; lo=0.2, hi=0.7, nodes=3,
                                  newton_steps=1, solves=5)

`T = T_total + w(lambda_min) (T_soft - T_total)` for one unit, value only, at the
modes of the last `ctsem_laplace_evaluate` at `values` under the total floor.
The Float64 reference for `_laplace_gated_term`; see the section comment.

The gate is exact and costs no likelihood sweep: `M - hi I` is block-factored
(the same elimination as `M`'s own, no fill-in), and a unit it accepts returns
`T_total` untouched. Only a flagged unit places `lambda_min`: by the engine's
own symmetric eigensolver for a unit up to `_LAPLACE_EIGEN_MAXDIM`, else by
`solves` steps of inverse iteration on the factor of `M` from a fixed start,
which converges at the eigengap's rate. Only a unit with `w > 0` evaluates the
soft rule.

Returns `(value, weight, lambda_min, flagged)`.
"""
function _laplace_gated_unit_reference(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, U::Integer; lo::Real=0.2, hi::Real=0.7,
    nodes::Integer=3, newton_steps::Integer=1, solves::Integer=5)
    theta = collect(Float64, values)
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    u = laplace.modes[U]
    blocks = laplace.units.blocks[U]
    inner = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
    isempty(u) && return (value=inner.value, weight=0.0, lambda_min=Inf, flagged=false)
    M = _laplace_unit_curvature(laplace, U, theta, Ls, u)
    ok, logdetM, _, _ = _laplace_block_factor(M, blocks)
    ok || return (value=NaN, weight=0.0, lambda_min=NaN, flagged=false)
    total = inner.value - max(logdetM, 0.0) / 2
    _laplace_exceeds_identity(M, blocks, hi) &&
        return (value=total, weight=0.0, lambda_min=Inf, flagged=false)
    d = length(u)
    dense = _laplace_block_dense(M, blocks, d)
    lambda = if d <= _LAPLACE_EIGEN_MAXDIM[]
        # A small unit is decomposed outright: no likelihood sweep either way,
        # and inverse iteration converges only as fast as the eigengap allows
        # -- five solves placed every unit of the weak-data fixture, where the
        # eigenvalue ratio is under 0.4, and were 4% off on a unit whose two
        # eigenvalues are 0.956 and 1.175.
        _ctsem_symeig(dense).values[1]
    else
        F = _ctsem_cholesky(copy(dense), d)
        x = fill(1 / sqrt(d), d)
        for _ in 1:Int(solves)
            y = F \ x
            x = y ./ sqrt(sum(abs2, y))
        end
        sum(x .* (dense * x))
    end
    w = _laplace_soft_weight(lambda, lo, hi)
    w == 0 && return (value=total, weight=0.0, lambda_min=lambda, flagged=true)
    soft = _laplace_soft_rule_unit(laplace, theta, U; nodes=nodes,
        ndirs=1, tau=0.0, recenter=true, newton_steps=newton_steps)
    return (value=total + w * (soft - total), weight=w, lambda_min=lambda,
        flagged=true)
end

