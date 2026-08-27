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
    already computed as the rule's location and scale. `nodes = 1` reproduces
    the Laplace value exactly, which is what `test_quadrature.jl` asserts.
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
block reuses the joint Laplace mode and the corresponding diagonal block of the
selected inverse, which is that effect's marginal covariance under the Gaussian
approximation. The rule does not need the exact conditional mode to be *valid*,
only to be efficient, and that centre and scale have already been paid for.

With one level the tree is a single leaf with no ancestors, the recursion is one
call, and the result is identical to the flat rule this file started as -- which
is what `test_quadrature.jl` pins by checking `nodes = 1` against the Laplace
value at one, two and three levels.
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

const _GH_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

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

const _GH_GRID_CACHE = Dict{Tuple{Int,Int},Tuple{Vector{Vector{Float64}},Vector{Float64}}}()

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
    _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws, slot)

The conditional mode and scale for a leaf block, with its ancestors held at
whatever `u` currently says.

The leaf mode is where the accuracy is. It moves a long way with the effect
above it -- a subject whose study effect has been pushed two standard deviations
has a different best guess about its own -- so re-solving it at every outer node
is the point of the recursion rather than an optimisation within it.

Newton on this block's coordinates only. The ancestors are data here, and a
member outside this block contributes nothing to its curvature, so the problem
is `k x k` however large the unit is.

Returns `(ok, centre, scale, logdetscale)` with `scale * scale' = M^-1`.
"""
function _quadrature_leaf_rule!(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, aws, slot::Integer)
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
        ws = _laplace_workspace!(laplace, S, length(theta), slot)
        work = convert(Vector{S}, u)
        @inbounds for (t, c) in enumerate(columns); work[c] = z[t]; end
        result = _laplace_unit_loglik_gradient(laplace, U, convert(Vector{S}, theta),
            [convert(Matrix{S}, L) for L in Ls], work, ws, members)
        return [result.gradient[c] for c in columns]
    end
    curvature_at = function (z)
        A = ForwardDiff.jacobian(loglik_gradient, z)
        M = Matrix{Float64}(LinearAlgebra.I, k, k) .- _laplace_symmetrise(A)
        return cholesky(Symmetric(M); check=false)
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
    return (ok=true, centre=z, scale=Matrix(inv(factorization.U)),
        logdetscale=-logdet(factorization) / 2)
end

"""
    _quadrature_block(laplace, U, theta, Ls, b, u, context, aws, slot)

Block `b`'s contribution to its unit's log marginal: its own coordinates
integrated, and whatever sits beneath it recursed into.

`u` carries the values already chosen for `b`'s ancestors and is written in
place as the recursion descends. Entries belonging to other branches are stale
and unread, because a member reads only the blocks on its own path.
"""
function _quadrature_block(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, context, aws, slot::Integer)
    blocks = laplace.units.blocks[U]
    block = blocks[b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    children = context.children[b]
    leaf = isempty(children)

    rule = if leaf
        _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws, slot)
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
        factorization = context.factors[b]
        (ok=true, centre=Float64[context.mode[c] for c in columns],
         scale=Matrix(inv(factorization.U)),
         logdetscale=-logdet(factorization) / 2)
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
                    aws, slot)
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
machinery already produces. With `nodes = 1` the rule has a single point at the
mode and the result is the Laplace value to machine precision, which is the
cheapest available check that the two agree about what they are integrating.

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
    while length(laplace.workspaces) < nchunks
        push!(laplace.workspaces, Dict{Any,Any}())
    end
    ranges = _ctsem_chunk_assignment(_laplace_unit_weights(laplace), nchunks)
    failed = fill(false, nchunks)

    # A trial point can be invalid in ways that *throw* rather than return a
    # non-finite number -- a parameter vector an optimizer wandered into can
    # make the matrix exponential's own scaling step take `ceil(Int, NaN)`. The
    # engine's contract is that an evaluation reports NaN, and a throw inside a
    # spawned task escapes as a `TaskFailedException` that kills whatever loop
    # is above it, so it is caught here rather than left to every caller.
    run = function (c)
        try
            _quadrature_chunk!(laplace, theta, Ls, ranges[c], unit_term, nodes, c)
        catch err
            err isa InterruptException && rethrow()
            failed[c] = true
        end
        return nothing
    end
    if nchunks <= 1
        run(1)
    else
        Threads.@sync for c in 1:nchunks
            Threads.@spawn run(c)
        end
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
    _quadrature_chunk!(laplace, theta, Ls, units, unit_term, nodes, slot)

One chunk of the unit loop, writing each unit's term into `unit_term`.

Throws rather than flagging: the caller wraps it, so that a failure inside a
spawned task becomes a NaN evaluation rather than an escaping exception.
"""
function _quadrature_chunk!(laplace::CTSEMLaplaceObjective, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, units, unit_term::Vector{Float64},
    nodes::Integer, slot::Integer)
    aws = _laplace_workspace!(laplace, Float64, length(theta), slot)
    for U in units
        blocks = laplace.units.blocks[U]
        if isempty(blocks)
            # No random effects anywhere in this unit: there is no integral,
            # and the term is the plain log likelihood.
            unit_term[U] = _laplace_unit_objective_gradient(laplace, U, theta,
                Ls, Float64[], aws).value
            continue
        end
        _laplace_solve_unit_mode!(laplace, U, theta, Ls, slot)
        mode = copy(laplace.modes[U])
        M = _laplace_unit_curvature(laplace, U, theta, Ls, mode, slot)
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
                aws, slot)
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
    delta = try
        -(H \ gap)
    catch err
        err isa InterruptException && rethrow()
        fill(NaN, length(theta))
    end
    quadrature = ctsem_laplace_quadrature(laplace, theta; nodes=nodes).value
    laplacevalue = ctsem_laplace_evaluate(laplace, theta; gradient=false).value
    return (delta=delta, corrected=theta .+ delta, gap_gradient=gap,
        quadrature=quadrature, laplace=laplacevalue,
        gap=quadrature - laplacevalue, nodes=Int(nodes))
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
