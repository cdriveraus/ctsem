using DataFrames, LinearAlgebra, Random

# `ctsem_evaluate_batch`: one bridge call in place of one per proposal draw.
#
# The point of the test is not the arithmetic -- `ctsem_evaluate_batch` is a
# loop over `_ctsem_probe_value`, which is already exercised elsewhere -- it is
# that the loop reproduces the column-by-column answer exactly, for each of the
# three objective shapes the R side hands it (`CTSEMObjective`,
# `CTSEMLaplaceObjective`, `CTSEMJointObjective`), including the `-Inf` on a
# column the model cannot evaluate.
#
# Fixtures are shared rather than rebuilt, the way test_quadrature.jl reuses
# test_laplace.jl's: `_adjoint_linear_1d_parameters` (adjoint_fixtures.jl) has
# two free parameters, so a batch actually varies across columns, and
# `_LAPLACE_LINEAR_OBJECTIVE` (laplace_fixtures.jl) is a real Laplace fit with
# an inner solve to fail. Included explicitly so this file stays independent
# of load order under both runtests.jl and runtests_parallel.jl.
isdefined(@__MODULE__, :_adjoint_linear_1d_parameters) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))
isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

function _evalbatch_marginal_objective()
    sp = _adjoint_linear_1d_parameters()
    subject_starts = [1, 4]
    times = [0.0, 0.5, 1.0, 0.0, 0.7, 1.3]
    data = reshape([0.1, -0.2, 0.05, 0.0, 0.2, -0.1], 1, :)
    ContinuousTimeSEM.ctsem_objective(sp, subject_starts, times, data)
end

@testset "batch matches column-by-column for the marginal objective, good and bad columns" begin
    objective = _evalbatch_marginal_objective()
    rng = MersenneTwister(401)
    X = 0.3 .* randn(rng, 2, 6)
    # One column the model cannot evaluate, by the same mechanism the R-side
    # guard tests use: a NaN reaches the arithmetic and every finite check
    # downstream fails, rather than relying on a particular exception.
    X[1, 3] = NaN
    expected = [ContinuousTimeSEM._ctsem_probe_value(objective, view(X, :, j))
                for j in 1:size(X, 2)]
    got = ContinuousTimeSEM.ctsem_evaluate_batch(objective, X)
    @test got == expected
    @test isfinite(expected[1]) && isfinite(expected[2])
    @test expected[3] == -Inf
    @test all(isfinite, got[[1, 2, 4, 5, 6]])

    # And through the generic `ctsem_evaluate` directly, for the columns it can
    # evaluate at all -- the batch is not a different computation from the one
    # `ctsem_evaluate_batch`'s docstring promises it dispatches through.
    for j in (1, 2, 4, 5, 6)
        @test got[j] == ContinuousTimeSEM.ctsem_evaluate(objective, X[:, j]; gradient=false).value
    end
end

@testset "batch matches column-by-column for the joint objective" begin
    objective = _evalbatch_marginal_objective()
    npar = 2
    joint = ContinuousTimeSEM.ctsem_joint_objective(objective, npar)
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    rng = MersenneTwister(402)
    X = 0.3 .* randn(rng, ndim, 5)
    expected = [ContinuousTimeSEM._ctsem_probe_value(joint, view(X, :, j))
                for j in 1:size(X, 2)]
    got = ContinuousTimeSEM.ctsem_evaluate_batch(joint, X)
    @test got == expected
    @test all(isfinite, got)
end

@testset "batch matches column-by-column for the laplace objective, including a non-converged column" begin
    objective = _LAPLACE_LINEAR_OBJECTIVE
    good = _LAPLACE_LINEAR_VALUES
    npar = length(good)
    rng = MersenneTwister(403)
    X = good .+ 0.02 .* randn(rng, npar, 4)
    X = hcat(X, good)
    # A column the model cannot evaluate, the same mechanism as the marginal
    # and joint cases above. `_ctsem_probe_value`'s laplace override has its
    # own further reason to return -Inf -- an inner solve that does not
    # converge, at a finite raw vector -- but that failure mode depends on the
    # fixture's own curvature and is not reliably reproduced by picking a raw
    # magnitude by hand; the NaN case already proves the batch route carries
    # whatever `_ctsem_probe_value` decides, for whichever reason.
    bad = copy(good); bad[1] = NaN
    X = hcat(X, bad)
    expected = [ContinuousTimeSEM._ctsem_probe_value(objective, view(X, :, j))
                for j in 1:size(X, 2)]
    got = ContinuousTimeSEM.ctsem_evaluate_batch(objective, X)
    @test got == expected
    @test all(isfinite, got[1:(end - 1)])
    @test got[end] == -Inf

    # The value at the fixture's own point is not merely finite but the number
    # `ctsem_evaluate` reports for it directly -- so the batch route has not
    # quietly picked up a different evaluation of the same objective.
    direct = ContinuousTimeSEM.ctsem_evaluate(objective, good; gradient=false).value
    @test got[5] == direct
end

@testset "an empty matrix is refused rather than crossing the bridge" begin
    objective = _evalbatch_marginal_objective()
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_evaluate_batch(objective, zeros(2, 0))
end

@testset "column order is not special" begin
    objective = _evalbatch_marginal_objective()
    joint = ContinuousTimeSEM.ctsem_joint_objective(objective, 2)
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    rng = MersenneTwister(404)
    X = 0.3 .* randn(rng, ndim, 5)
    forward = ContinuousTimeSEM.ctsem_evaluate_batch(joint, X)
    reversed = ContinuousTimeSEM.ctsem_evaluate_batch(joint, X[:, end:-1:1])
    @test forward == reverse(reversed)
end
