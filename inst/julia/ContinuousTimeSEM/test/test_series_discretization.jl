using LinearAlgebra
using Random
using ForwardDiff

# The series route of the discretisation (`series_discretization.jl`): the
# intercept and process-noise integrals taken directly where the closed forms
# would divide by a pivot that is small against the interval.
#
# The kernels are checked against an independent answer -- the same integrals
# summed in 256-bit arithmetic with no scaling -- rather than against the closed
# forms, which are exactly what fails in the cases that matter here. The
# likelihood is checked through a drift that is exactly singular, where the
# closed forms return NaN, and the adjoint against forward mode there to first
# and second order.

isdefined(@__MODULE__, :_adjoint_test_dataframe) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))

# int_0^h e^{J s} ds and int_0^h e^{J s} Q e^{J' s} ds, summed term by term in
# BigFloat: h sum_j (J h)^j / (j+1)! and sum_j h^(j+1) / (j+1)! L^j(Q), with
# L(W) = J W + W J'. Four hundred terms is far past convergence for every
# case below.
function _series_reference_phi(J, h; terms=400)
    Jb = BigFloat.(J)
    hb = BigFloat(h)
    term = Matrix{BigFloat}(I, size(J)...) * hb
    total = copy(term)
    for j in 1:terms
        term = Jb * term * hb / (j + 1)
        total += term
    end
    return total
end

function _series_reference_noise(J, Q, h; terms=400)
    Jb = BigFloat.(J)
    hb = BigFloat(h)
    L = BigFloat.(Q)
    coefficient = hb
    total = coefficient * L
    for j in 1:terms
        L = Jb * L + L * Jb'
        coefficient = coefficient * hb / (j + 1)
        total += coefficient * L
    end
    return total
end

_series_relerr(a, b) = maximum(abs.(Float64.(a) .- Float64.(b))) /
    max(maximum(abs.(Float64.(b))), floatmin(Float64))

@testset "series kernels match a 256-bit reference" begin
    rng = MersenneTwister(1)
    cases = [
        ("stable", [-0.8 0.3 0.1; 0.2 -0.5 0.05; -0.1 0.2 -1.2], 0.7),
        ("singular", [-1.0 0.5; 0.5 -0.25], 0.8),
        ("saddle", [0.0 1.0; 1.0 0.0], 0.6),
        ("stiff with a near-zero eigenvalue", [-40.0 0.0; 3.0 -1e-9], 1.3),
        ("short interval", [-0.5 0.2; 0.1 -0.3], 1e-6),
        ("long interval", [-0.5 0.2; 0.1 -0.3], 25.0),
        ("explosive", [0.4 0.1; 0.0 0.2], 3.0),
    ]
    for (label, J, h) in cases
        p = size(J, 1)
        # One wider than needed, so the kernels' leading-block indexing is
        # exercised as the filter uses it.
        buf = ContinuousTimeSEM.SeriesKernelBuffer(Float64, p + 1)
        Phi = copy(ContinuousTimeSEM._series_intercept!(buf, J, h, p))
        @test _series_relerr(Phi, _series_reference_phi(J, h)) < 1e-13
        B = randn(rng, p, p)
        Q = B * B' + 0.1I
        V = copy(ContinuousTimeSEM._series_noise!(buf, J, Q, h, p))
        @test _series_relerr(V, _series_reference_noise(J, Q, h)) < 1e-13
        @test V == V'
    end
end

@testset "series kernels agree with the closed forms where those are well conditioned" begin
    J = [-0.8 0.3 0.1; 0.2 -0.5 0.05; -0.1 0.2 -1.2]
    h = 0.7
    Q = [0.3 0.05 0.0; 0.05 0.2 0.01; 0.0 0.01 0.4]
    buf = ContinuousTimeSEM.SeriesKernelBuffer(Float64, 3)
    E = exp(J * h)
    @test ContinuousTimeSEM._series_intercept!(buf, J, h, 3) ≈ J \ (E - I) rtol = 1e-13
    X = lyap(J, Q)
    @test ContinuousTimeSEM._series_noise!(buf, J, Q, h, 3) ≈ X - E * X * E' rtol = 1e-12
end

@testset "series pullbacks match forward mode" begin
    rng = MersenneTwister(2)
    phi_of(Jv, h, C) = begin
        p = size(C, 1)
        buf = ContinuousTimeSEM.SeriesKernelBuffer(eltype(Jv), p)
        sum(ContinuousTimeSEM._series_intercept!(buf, reshape(Jv, p, p), h, p) .* C)
    end
    noise_of(v, h, C) = begin
        p = size(C, 1)
        Qh = reshape(v[p*p+1:end], p, p)
        buf = ContinuousTimeSEM.SeriesKernelBuffer(eltype(v), p)
        sum(ContinuousTimeSEM._series_noise!(buf, reshape(v[1:p*p], p, p),
            (Qh + Qh') / 2, h, p) .* C)
    end
    for (J, h) in [([-1.0 0.5; 0.5 -0.25], 0.8),
                   ([-0.8 0.3 0.1; 0.2 -0.5 0.05; -0.1 0.2 -1.2], 2.7),
                   ([0.0 1.0; 1.0 0.0], 0.6), ([-0.5 0.2; 0.1 -0.3], 1e-3)]
        p = size(J, 1)
        C = randn(rng, p, p)
        buf = ContinuousTimeSEM.SeriesKernelBuffer(Float64, p)

        forward = ForwardDiff.gradient(v -> phi_of(v, h, C), vec(J))
        ContinuousTimeSEM._series_intercept!(buf, J, h, p)
        Jbar = zeros(p, p)
        ContinuousTimeSEM._series_intercept_pullback!(Jbar, buf, p, C)
        @test vec(Jbar) ≈ forward rtol = 1e-12 atol = 1e-14

        B = randn(rng, p, p)
        Q = B * B' + 0.1I
        Cs = (C + C') / 2
        forward = ForwardDiff.gradient(v -> noise_of(v, h, Cs), vcat(vec(J), vec(Q)))
        ContinuousTimeSEM._series_noise!(buf, J, Q, h, p)
        Jbar = zeros(p, p)
        Qbar = zeros(p, p)
        ContinuousTimeSEM._series_noise_pullback!(Jbar, Qbar, buf, p, Cs)
        @test vec(Jbar) ≈ forward[1:p*p] rtol = 1e-12 atol = 1e-14
        # Forward mode differentiates through (Q + Q') / 2.
        @test vec((Qbar + Qbar') / 2) ≈ forward[p*p+1:end] rtol = 1e-12 atol = 1e-14
    end
end

# The tables (`SeriesTable`) sum the same truncated series from powers built
# once per drift, and their pullback releases a run of intervals' cotangents in
# one pass through the power chains: against the reference, against Horner's
# kernels, and against the sum of Horner's pullbacks over the run.
@testset "series tables reproduce the kernels and their pullbacks" begin
    CT = ContinuousTimeSEM
    rng = MersenneTwister(3)
    cases = [
        ("stable", [-0.8 0.3 0.1; 0.2 -0.5 0.05; -0.1 0.2 -1.2], (1e-3, 0.3, 0.7, 2.7)),
        ("singular", [-1.0 0.5; 0.5 -0.25], (0.8, 4.0)),
        ("stiff with a near-zero eigenvalue", [-40.0 0.0; 3.0 -1e-9], (1e-4, 1.3)),
        ("long interval", [-0.5 0.2; 0.1 -0.3], (25.0, 0.01)),
        ("explosive", [0.4 0.1; 0.0 0.2], (3.0,)),
    ]
    for (label, J, hs) in cases
        p = size(J, 1)
        n = p + 1
        Jn = zeros(n, n)
        Jn[1:p, 1:p] = J
        B = randn(rng, p, p)
        Q = B * B' + 0.1I
        Qn = zeros(n, n)
        Qn[1:p, 1:p] = Q
        noise = CT.SeriesTable(Float64, n, true)
        intercept = CT.SeriesTable(Float64, n, false)
        @test !CT._series_table_key!(noise, Jn, Qn, p)
        @test !CT._series_table_key!(intercept, Jn, nothing, p)
        @test CT._series_table_key!(noise, Jn, Qn, p)
        noise_bar = CT.SeriesTableBar(Float64)
        intercept_bar = CT.SeriesTableBar(Float64)
        Jb_sum = zeros(n, n)
        Qb_sum = zeros(n, n)
        Ji_sum = zeros(n, n)
        for h in hs
            horner = CT.SeriesKernelBuffer(Float64, n)
            table = CT.SeriesKernelBuffer(Float64, n)
            V = copy(CT._series_noise!(horner, Jn, Qn, h, p))
            Vt = copy(CT._series_noise_table!(table, noise, h, p))
            @test _series_relerr(Vt, _series_reference_noise(J, Q, h)) < 1e-13
            @test Vt ≈ V rtol = 1e-13
            @test Vt == Vt'
            C = randn(rng, p, p)
            CT._series_noise_pullback!(Jb_sum, Qb_sum, horner, p, C)
            CT._series_noise_table_pullback!(noise_bar, table, noise, p, C)

            Phi = copy(CT._series_intercept!(horner, Jn, h, p))
            Phit = copy(CT._series_intercept_table!(table, intercept, h, p))
            @test _series_relerr(Phit, _series_reference_phi(J, h)) < 1e-13
            @test Phit ≈ Phi rtol = 1e-13
            C = randn(rng, p, p)
            CT._series_intercept_pullback!(Ji_sum, horner, p, C)
            CT._series_intercept_table_pullback!(intercept_bar, table, intercept, p, C)
        end
        Jb = zeros(n, n)
        Qb = zeros(n, n)
        CT._series_table_pullback!(Jb, Qb, noise, noise_bar)
        @test Jb ≈ Jb_sum rtol = 1e-11 atol = 1e-13
        @test Qb ≈ Qb_sum rtol = 1e-11 atol = 1e-13
        Ji = zeros(n, n)
        CT._series_table_pullback!(Ji, nothing, intercept, intercept_bar)
        @test Ji ≈ Ji_sum rtol = 1e-11 atol = 1e-13
        # Released: a second release adds nothing, and the next run starts
        # from zero however long the last one was.
        @test !noise_bar.pending
        CT._series_table_pullback!(Jb, Qb, noise, noise_bar)
        @test Jb ≈ Jb_sum rtol = 1e-11 atol = 1e-13
        table = CT.SeriesKernelBuffer(Float64, n)
        CT._series_noise_table!(table, noise, hs[1], p)
        C = randn(rng, p, p)
        CT._series_noise_table_pullback!(noise_bar, table, noise, p, C)
        horner = CT.SeriesKernelBuffer(Float64, n)
        CT._series_noise!(horner, Jn, Qn, hs[1], p)
        J1 = zeros(n, n)
        Q1 = zeros(n, n)
        CT._series_noise_pullback!(J1, Q1, horner, p, C)
        J2 = zeros(n, n)
        Q2 = zeros(n, n)
        CT._series_table_pullback!(J2, Q2, noise, noise_bar)
        @test J2 ≈ J1 rtol = 1e-11 atol = 1e-13
        @test Q2 ≈ Q1 rtol = 1e-11 atol = 1e-13
    end
end

# Two states whose cross effect is one parameter in both cells, as in the R
# fixture that found this: det(DRIFT) = a1 a2 - c^2, which is exactly zero at
# a1 = -1, a2 = -0.25, c = 0.5 -- an eigenvalue of zero, so both the intercept
# solve and the Lyapunov operator are singular there, not merely ill
# conditioned. CINT and DIFFUSION are free so both integrals reach the
# gradient.
function _series_singular_parameters()
    df = _adjoint_test_dataframe(
        drift=[-1.0 0.5; 0.5 -0.25], jax=[-1.0 0.5; 0.5 -0.25],
        cint=[0.0; 0.0;;], diffusion=[0.2 0.0; 0.0 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
        t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "param[1]"), (:JAx, 1, 1) => (1, "param[1]"),
            (:DRIFT, 2, 2) => (2, "param[2]"), (:JAx, 2, 2) => (2, "param[2]"),
            (:DRIFT, 1, 2) => (3, "param[3]"), (:JAx, 1, 2) => (3, "param[3]"),
            (:DRIFT, 2, 1) => (3, "param[3]"), (:JAx, 2, 1) => (3, "param[3]"),
            (:CINT, 1, 1) => (4, "param[4]"), (:CINT, 2, 1) => (5, "param[5]"),
            (:DIFFUSION, 1, 1) => (6, "log1p_exp(param[6])"),
            (:DIFFUSION, 2, 1) => (7, "param[7]"),
        ),
    )
    ekf_from_data_frame(df)
end

_SERIES_SP = _series_singular_parameters()
_SERIES_DATA = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.25], 2, :)
_SERIES_TIMES = [0.0, 0.5, 1.2, 2.0]
_series_values(c) = [-1.0, -0.25, c, 0.3, -0.2, 0.1, 0.05]

@testset "a two-state likelihood through a singular drift" begin
    objective = ContinuousTimeSEM.ctsem_objective(_SERIES_SP, [1], _SERIES_TIMES,
        _SERIES_DATA)
    for c in (0.3, sqrt(0.25 - 1e-4), sqrt(0.25 - 1e-8), 0.5, sqrt(0.25 + 1e-8), 0.6)
        values = _series_values(c)
        @test isfinite(objective(values))
        result = ContinuousTimeSEM.ctsem_validate_forward_gradient(objective, values)
        @test result.relative_error < 1e-6
        @test result.adjoint_relative_error < 1e-9
    end

    # Smooth through the singular point: second differences on a 1e-6 grid are
    # at rounding level. The closed forms give NaN at its centre and second
    # differences of 2e-9 around it, and at det = 1e-8 a gradient with no
    # correct digit (0.9 relative error against finite differences).
    values = [objective(_series_values(0.5 + i * 1e-6)) for i in -5:5]
    @test maximum(abs.(diff(diff(values)))) < 1e-11

    # Which route ran, from the counters: none of the series at a
    # well-conditioned drift, both kernels at the singular one.
    ContinuousTimeSEM.ctsem_reset_opcounts!()
    objective(_series_values(0.3))
    counts = ContinuousTimeSEM.ctsem_opcounts()
    @test counts.series_intercept == 0
    @test counts.series_noise == 0
    objective(_series_values(0.5))
    counts = ContinuousTimeSEM.ctsem_opcounts()
    @test counts.series_intercept > 0
    @test counts.series_noise > 0
end

@testset "forward over reverse through the series" begin
    objective = ContinuousTimeSEM.ctsem_objective(_SERIES_SP, [1], _SERIES_TIMES,
        _SERIES_DATA)
    values = _series_values(0.5)
    H = ContinuousTimeSEM.ctsem_hessian(objective, values)
    @test H ≈ ForwardDiff.hessian(objective, values) rtol = 1e-9
end

@testset "the series everywhere reproduces the closed forms" begin
    values = _series_values(0.3)
    closed = ContinuousTimeSEM.ctsem_objective(_SERIES_SP, [1], _SERIES_TIMES,
        _SERIES_DATA)
    value = closed(values)
    gradient = ContinuousTimeSEM.ctsem_adjoint_gradient(closed, values).gradient
    old = ContinuousTimeSEM._CTSEM_SERIES_BELOW[]
    try
        ContinuousTimeSEM.ctsem_set_series_discretization_below!(Inf)
        series = ContinuousTimeSEM.ctsem_objective(_SERIES_SP, [1], _SERIES_TIMES,
            _SERIES_DATA)
        @test series(values) ≈ value rtol = 1e-12
        @test ContinuousTimeSEM.ctsem_adjoint_gradient(series, values).gradient ≈
            gradient rtol = 1e-10
        # Substeps, so one row's tape holds several series records.
        substepped = ContinuousTimeSEM.ctsem_objective(_SERIES_SP, [1], _SERIES_TIMES,
            _SERIES_DATA, zeros(0, 4), zeros(1, 0), 0.2)
        result = ContinuousTimeSEM.ctsem_validate_forward_gradient(substepped, values)
        @test result.adjoint_relative_error < 1e-9
    finally
        ContinuousTimeSEM.ctsem_set_series_discretization_below!(old)
    end
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_set_series_discretization_below!(-1)
end
