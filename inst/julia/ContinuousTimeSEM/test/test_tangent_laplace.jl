# The hand-coded Laplace sweeps (`tangent_laplace.jl`) against forward mode
# through the filter at the duals and jets they replaced, their oracles
# (`dual_oracles.jl`): the curvature's dense and blocked routes, one to
# three levels, effects on the means and on the dynamics and noise, both at
# once by the split route, a curved mean, a drift evaluated by row, and
# categorical indicators.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))
isdefined(@__MODULE__, :_adjoint_test_dataframe) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))
isdefined(@__MODULE__, :_oracle_hessian) ||
    include(joinpath(@__DIR__, "dual_oracles.jl"))

# One latent read by an ordinal item in three categories, a count and an item
# censored at -0.25 and 0.3 (`tangent_categorical.jl`): the thresholds, the
# count's dispersion and intercept and the censored item's sd and intercept
# free, and random effects on the initial mean, the intercept and the drift. A
# fourth item, ordinal on the first's thresholds, reads no state, so its
# predictor has no variance: the filter adds its likelihood at its mean, which
# shares the censored item's intercept, and moves nothing.
const _TAN_CATEGORICAL_OBJECTIVE = let
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.3;;],
        lambda=[1.0; 0.6; 0.8; 0.0;;], jy=[1.0; 0.6; 0.8; 0.0;;], manifestmeans=zeros(4, 1),
        manifestvar=Matrix(0.3I, 4, 4), t0var=[0.5;;], t0means=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "param[1]"),
            (:CINT, 1, 1) => (2, "param[2]"),
            (:DRIFT, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:JAx, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:MANIFESTMEANS, 2, 1) => (5, "param[5]"),
            (:THRESHOLDS, 1, 1) => (6, "param[6]"),
            (:THRESHOLDS, 1, 2) => (7, "log1p_exp(param[7])"),
            (:MANIFESTVAR, 2, 2) => (8, "log1p_exp(param[8])"),
            (:MANIFESTVAR, 3, 3) => (9, "log1p_exp(param[9])"),
            (:MANIFESTMEANS, 3, 1) => (10, "param[10]"),
            (:MANIFESTMEANS, 4, 1) => (10, "0.5 * param[10] - 0.2"),
            (:THRESHOLDS, 4, 1) => (6, "param[6]"),
            (:THRESHOLDS, 4, 2) => (7, "log1p_exp(param[7])"),
        ),
        extra=[:THRESHOLDS => zeros(4, 2)])
    sp = ekf_from_data_frame(df; manifesttype=[2, 3, 4, 2], ncategories=[3, 0, 0, 3],
        censormin=[-Inf, -Inf, -0.25, -Inf], censormax=[Inf, Inf, 0.3, Inf])
    starts, times, data = _laplace_test_data(6, 4; seed=4)
    for c in axes(data, 2)
        v = data[:, c]
        data[1, c] = 1.0 + (v[1] > -0.1) + (v[1] > 0.15)
        data[2, c] = round(2 + 4 * v[2])
        data[3, c] = clamp(v[3], -0.25, 0.3)
        data[4, c] = 1.0 + (v[4] > -0.2) + (v[4] > 0.1)
        c % 5 == 2 && (data[2, c] = NaN)
        c % 7 == 3 && (data[1, c] = NaN)
    end
    ctsem_objective(sp, starts, times, data)
end
_fresh_categorical() = (ctsem_laplace_objective(_TAN_CATEGORICAL_OBJECTIVE, [1, 2, 3],
    [11, 12, 13], [14, 15, 16], [1.0, 1.0, 1.0]),
    [0.1, -0.1, 0.2, -0.3, 0.2, -0.2, 0.3, -0.8, -1.0, 0.05, -0.3, -0.4, -0.5, 0.1,
        -0.1, 0.2])
# One effect, the intercept's: a unit's moments along its one lane
# (`_tan_cat_moments_lanes!`).
_fresh_categorical1() = (ctsem_laplace_objective(_TAN_CATEGORICAL_OBJECTIVE, [2], [11],
    Int[], [1.0]), [0.1, -0.1, 0.2, -0.3, 0.2, -0.2, 0.3, -0.8, -1.0, 0.05, -0.3])

# The linear model with its initial mean's transform curved: an effect there
# moves the innovations nonlinearly, so Gauss-Newton is not its curvature, and
# the objective says so as R's would, `affine` naming only the intercept.
const _TAN_CURVED_MEAN_OBJECTIVE = let
    df = _laplace_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.4;;], t0means=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "2 * log1p_exp(param[1])"),
            (:CINT, 1, 1) => (2, "param[2]"),
            (:DRIFT, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:JAx, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:MANIFESTMEANS, 1, 1) => (5, "param[5]"),
        ),
    )
    subject_starts, times, data = _laplace_test_data(6, 1)
    ctsem_objective(ekf_from_data_frame(df), subject_starts, times, data)
end

@testset "hand-coded Laplace curvature against the dual" begin
    C = ContinuousTimeSEM
    function curvatures(laplace, values, seed)
        Ls = C._laplace_popchols(values, laplace.spec)
        rng = Random.MersenneTwister(seed)
        out = Tuple{Matrix{Float64},Matrix{Float64}}[]
        for U in eachindex(laplace.units.members)
            d = length(laplace.units.blocks[U]) == 0 ? 0 :
                sum(b -> b.size, laplace.units.blocks[U])
            d == 0 && continue
            u = 0.4 .* randn(rng, d)
            hand = C._laplace_block_dense(C._laplace_unit_curvature(laplace, U, values, Ls,
                u), laplace.units.blocks[U], d)
            # The curvature is minus the unit objective's Hessian.
            push!(out, (hand, -_oracle_unit_hessian(laplace, U, values, Ls, u)))
        end
        return out
    end
    function served(laplace, values)
        ctx = get(laplace.workspaces[C._laplace_slot()],
            (:tan_curvature, length(values), C._HYBRID_PROFILE), nothing)
        return ctx !== nothing
    end
    cases = [
        ("means", () -> _fresh_linear()),
        # Every effect in the means affinely: Gauss-Newton for the whole member.
        ("means, quadratic", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]; quadratic=true), copy(_LAPLACE_LINEAR_VALUES))),
        # Two effects on the dynamics and noise beside three in the means: the
        # split, one full sweep and one forward (`_tan_laplace_split_member!`).
        ("means beside dynamics, split", () -> (ctsem_laplace_objective(
            _LAPLACE_LINEAR_OBJECTIVE, [1, 2, 5, 3, 4], collect(6:10), collect(11:20),
            ones(5); affine_level=[1, 1, 1], affine_effect=[1, 2, 3]),
            vcat(_LAPLACE_LINEAR_VALUES[1:5], [-0.3, -0.15, -0.2, -0.4, -0.35],
                0.1 .* sin.(1:10)))),
        # An intercept beside an initial mean whose transform is curved: the
        # objective's `affine` keeps the member off Gauss-Newton.
        ("curved mean", () -> (ctsem_laplace_objective(_TAN_CURVED_MEAN_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]; affine_level=[1], affine_effect=[2]),
            copy(_LAPLACE_LINEAR_VALUES))),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
        ("two levels", () -> _fresh_twolevel()),
        ("three levels", () -> _fresh_threelevel()),
        ("categorical", () -> _fresh_categorical()),
        ("categorical, one effect", () -> _fresh_categorical1()),
    ]
    threshold = C._LAPLACE_BLOCK_THRESHOLD[]
    try
        for (name, make) in cases, blocked in (false, true)
            C.ctsem_set_block_threshold!(blocked ? 1 : threshold)
            laplace, values = make()
            pairs = curvatures(laplace, values, 3)
            @testset "$name, $(blocked ? "blocked" : "dense")" begin
                @test served(laplace, values)
                @test !isempty(pairs)
                for (hand, dual) in pairs
                    @test all(isfinite, hand)
                    @test maximum(abs, hand .- dual) <= 1e-9 * max(1.0, maximum(abs, dual))
                end
            end
        end
    finally
        C.ctsem_set_block_threshold!(threshold)
    end
end

@testset "hand-coded seeded sweeps against the dual" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("drift by row", () -> _fresh_nonlinear()),
        ("three levels", () -> _fresh_threelevel()),
        ("categorical", () -> _fresh_categorical()),
    ]
    for (name, make) in cases
        laplace, values = make()
        Ls = C._laplace_popchols(values, laplace.spec)
        rng = Random.MersenneTwister(7)
        @testset "$name" begin
            for U in eachindex(laplace.units.members)
                d = sum(b -> b.size, laplace.units.blocks[U]; init=0)
                u = 0.4 .* randn(rng, d)
                nm = length(laplace.units.members[U])
                for d1 in (randn(rng, length(values)), zeros(length(values)))
                    hand = C._laplace_unit_seeded_gradient(laplace, U, values, Ls, u,
                        1:nm, d1)
                    h0, h1 = copy(hand.d0), copy(hand.d1c)
                    dual = _oracle_seeded_gradient(laplace, U, values, Ls, u, 1:nm, d1)
                    @test hand.ok && dual.ok
                    @test maximum(abs, h0 .- dual.d0) <= 1e-9 * max(1.0, maximum(abs, dual.d0))
                    @test maximum(abs, h1 .- dual.d1c) <= 1e-9 * max(1.0, maximum(abs, dual.d1c))
                end
            end
        end
    end
end

@testset "hand-coded sequential predictions' passes against the dual" begin
    # Each subject's seeded plug-in pass (`_tan_laplace_trace!`), which the
    # sequential predictions build from, against the filter at the curvature's
    # dual on the same duals: the moments and row likelihoods it records, values
    # and lanes.
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
        ("categorical", () -> _fresh_categorical()),
    ]
    for (name, make) in cases
        laplace, values = make()
        out = C.ctsem_laplace_sequential(laplace, values)
        served = get(laplace.workspaces[C._laplace_slot()],
            (:tan_curvature, length(values), C._HYBRID_OFF), nothing) !== nothing
        S = C._laplace_curvature_dual(laplace)
        W = ForwardDiff.npartials(S)
        theta = C._laplace_extend(laplace, collect(Float64, values))
        Ls = C._laplace_popchols(theta, laplace.spec)
        sp = laplace.objective.params
        ews = C._laplace_workspace!(laplace, Float64, length(theta)).ekf_ws
        rng = Random.MersenneTwister(5)
        vals(A) = map(ForwardDiff.value, A)
        lane(A, w) = map(z -> ForwardDiff.partials(z, w), A)
        @testset "$name" begin
            @test served
            @test any(isfinite, out.llrow)
            for U in eachindex(laplace.units.members)
                i = laplace.units.members[U][1]
                u = 0.3 .* randn(rng, laplace.units.dims[U])
                shifted = C._laplace_member_values!(zeros(length(theta)), theta,
                    laplace.spec, Ls, u, laplace.units.offsets[U][1])
                dirs = randn(rng, length(theta), W)
                x = S[S(shifted[t], ForwardDiff.Partials(ntuple(w -> dirs[t, w], Val(W))))
                    for t in eachindex(theta)]
                T = size(laplace.objective.subject_objectives[i].data, 2)
                function trace()
                    tr = C.CTSEMKalmanTrace(S, ews.state_dim, ews.manifest_dim, T, 1,
                        length(sp.mutables))
                    tr.offset = 0
                    tr.current_subject = 1
                    return tr
                end
                th, to = trace(), trace()
                okh = C._tan_laplace_trace!(th, laplace, theta, i, x, W)
                @test okh == _oracle_trace!(to, laplace, i, x)
                okh || continue
                for f in (:eta, :etacov, :y, :ycov, :llrow)
                    a, b = getfield(th, f), getfield(to, f)
                    for get in (vals, (A -> lane(A, 1)), (A -> lane(A, W)))
                        ga, gb = get(a), get(b)
                        @test isequal(isnan.(ga), isnan.(gb))
                        ok = .!isnan.(gb)
                        @test maximum(abs, ga[ok] .- gb[ok]; init=0.0) <=
                            1e-9 * max(1.0, maximum(abs, gb[ok]; init=0.0))
                    end
                end
            end
        end
    end
end

@testset "hand-coded mode cross derivatives and gated passes against the dual" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
        ("two levels", () -> _fresh_twolevel()),
        ("categorical", () -> _fresh_categorical()),
    ]
    rel(a, b) = maximum(abs, a .- b) / max(1.0, maximum(abs, b))
    for (name, make) in cases
        laplace, values = make()
        C.ctsem_laplace_evaluate(laplace, values; gradient=false)
        theta = C._laplace_extend(laplace, collect(Float64, values))
        spec = laplace.spec
        Ls = C._laplace_popchols(theta, spec)
        dL = C._laplace_level_chol_derivatives(theta, spec)
        positions = Vector{Int}[C._laplace_level_positions(spec, l)
            for l in eachindex(spec.levels)]
        rng = Random.MersenneTwister(9)
        @testset "$name" begin
            # The mode Jacobian runs whole on the cross derivatives.
            @test all(J -> all(isfinite, J), C.ctsem_laplace_mode_jacobian(laplace, values))
            for U in eachindex(laplace.modes)
                z = copy(laplace.modes[U])
                isempty(z) && continue
                # The mode's cross derivatives (`_laplace_mode_cross`).
                H = C._laplace_mode_cross(laplace, U, theta, Ls, dL, positions, z, nothing)
                Hd = _oracle_mode_cross(laplace, U, theta, Ls, dL, positions, z, nothing)
                @test all(isfinite, H)
                @test rel(H, Hd) <= 1e-9
                # The gated rule's first-order pass and trace jet pass.
                G = C._GatedThetaGradient(laplace, U, theta, Ls, dL, positions, nothing)
                beta = randn(rng, length(z))
                r = C._gated_directional(G, z, beta)
                rd = _oracle_gated_directional(G, z, beta)
                @test rel(r.theta, rd.theta) <= 1e-8
                @test rel(r.z, rd.z) <= 1e-8
                A = randn(rng, length(z), length(z))
                E = C._ctsem_symeig(A .+ A')
                for chunk in ([1], collect(1:length(z)))
                    W = C._ctsem_jet_width(length(chunk))
                    tb, zb = zeros(length(theta)), zeros(length(z))
                    @test C._laplace_gated_jet_pass!(tb, zb, G, z, E, chunk, 1.0, Val(W))
                    jd = _oracle_gated_jet(G, z, E, chunk, W)
                    @test rel(-tb, jd.theta) <= 1e-8
                    @test rel(-zb, jd.z) <= 1e-8
                end
            end
            # The gated adjoint runs whole on those passes.
            @test any(!isnothing, [C._laplace_gated_unit_adjoint(laplace, U, values)
                for U in eachindex(laplace.modes)])
        end
    end
end

@testset "hand-coded trace jets against the jet" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("means beside dynamics", () -> (ctsem_laplace_objective(
            _LAPLACE_LINEAR_OBJECTIVE, [1, 2, 5, 3, 4], collect(6:10), collect(11:20),
            ones(5)), vcat(_LAPLACE_LINEAR_VALUES[1:5], [-0.3, -0.15, -0.2, -0.4, -0.35],
            0.1 .* sin.(1:10)))),
        ("curved mean", () -> (ctsem_laplace_objective(_TAN_CURVED_MEAN_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("two levels", () -> _fresh_twolevel()),
        ("three levels", () -> _fresh_threelevel()),
        # Cells by row: the engine's own code at the jet, a chunk at a time.
        ("drift by row", () -> _fresh_nonlinear()),
        ("categorical", () -> _fresh_categorical()),
        ("categorical, one effect", () -> _fresh_categorical1()),
    ]
    try
        for (name, make) in cases
            laplace, values = make()
            Ls = C._laplace_popchols(values, laplace.spec)
            rng = Random.MersenneTwister(11)
            npar = length(values)
            @testset "$name" begin
                for U in eachindex(laplace.units.members)
                    d = sum(b -> b.size, laplace.units.blocks[U]; init=0)
                    u = 0.4 .* randn(rng, d)
                    nm = length(laplace.units.members[U])
                    # One direction (the narrow jet), three, and seven (two
                    # chunks of the sweep's width).
                    for nd in (1, 3, 7)
                        D = zeros(npar, C._CTSEM_JET_WIDTH)
                        D[:, 1:nd] .= randn(rng, npar, nd)
                        hand = C._laplace_unit_jet_gradient(laplace, U, values, Ls, u,
                            1:nm, D, nd)
                        h0, hS = copy(hand.d0), copy(hand.dS)
                        jet = _oracle_jet_gradient(laplace, U, values, Ls, u, 1:nm, D, nd,
                            C._ctsem_jet_width(nd))
                        @test hand.ok && jet.ok
                        @test maximum(abs, h0 .- jet.d0) <=
                            1e-9 * max(1.0, maximum(abs, jet.d0))
                        @test maximum(abs, hS .- jet.dS) <=
                            1e-9 * max(1.0, maximum(abs, jet.dS))
                    end
                end
                served = get(laplace.workspaces[C._laplace_slot()],
                    (:tan_jet, npar, C.CTSEMTraceJet{C._CTSEM_JET_WIDTH}, C._HYBRID_PROFILE),
                    nothing) !== nothing
                @test served
            end
        end
    finally
    end
end

@testset "hand-coded Laplace sweeps on the series' tables against the dual and the jet" begin
    # Every interval on the series (`ctsem_set_series_discretization_below!`):
    # the curvature, the seeded sweeps and the trace jets take it by its tables
    # (`TanSeriesTable`), whose cotangents each member's sweep releases at its
    # flush; the oracles take the filter's own series at the dual and the jet.
    C = ContinuousTimeSEM
    below = C._CTSEM_SERIES_BELOW[]
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("categorical", () -> _fresh_categorical()),
    ]
    tol(a, b) = maximum(abs, a .- b) <= 1e-9 * max(1.0, maximum(abs, b))
    try
        C.ctsem_set_series_discretization_below!(Inf)
        counts = C._CTSEM_OPCOUNT
        for (name, make) in cases
            laplace, values = make()
            Ls = C._laplace_popchols(values, laplace.spec)
            rng = Random.MersenneTwister(5)
            npar = length(values)
            @testset "$name" begin
                for U in eachindex(laplace.units.members)
                    d = sum(b -> b.size, laplace.units.blocks[U]; init=0)
                    d == 0 && continue
                    u = 0.4 .* randn(rng, d)
                    nm = length(laplace.units.members[U])
                    hand = C._laplace_block_dense(C._laplace_unit_curvature(laplace, U,
                        values, Ls, u), laplace.units.blocks[U], d)
                    before = (counts.series_intercept[], counts.series_noise[])
                    dual = -_oracle_unit_hessian(laplace, U, values, Ls, u)
                    # The filter took both series.
                    @test counts.series_intercept[] > before[1]
                    @test counts.series_noise[] > before[2]
                    @test tol(hand, dual)
                    d1 = randn(rng, npar)
                    hs = C._laplace_unit_seeded_gradient(laplace, U, values, Ls, u, 1:nm, d1)
                    h0, h1 = copy(hs.d0), copy(hs.d1c)
                    ds = _oracle_seeded_gradient(laplace, U, values, Ls, u, 1:nm, d1)
                    @test hs.ok && ds.ok
                    @test tol(h0, ds.d0)
                    @test tol(h1, ds.d1c)
                    for nd in (1, 7)
                        D = zeros(npar, C._CTSEM_JET_WIDTH)
                        D[:, 1:nd] .= randn(rng, npar, nd)
                        hj = C._laplace_unit_jet_gradient(laplace, U, values, Ls, u, 1:nm,
                            D, nd)
                        j0, jS = copy(hj.d0), copy(hj.dS)
                        jet = _oracle_jet_gradient(laplace, U, values, Ls, u, 1:nm, D, nd,
                            C._ctsem_jet_width(nd))
                        @test hj.ok && jet.ok
                        @test tol(j0, jet.d0)
                        @test tol(jS, jet.dS)
                    end
                end
            end
        end
    finally
        C.ctsem_set_series_discretization_below!(below)
    end
end
