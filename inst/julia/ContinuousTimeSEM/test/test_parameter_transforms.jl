using ComponentArrays
using DataFrames
using ForwardDiff

# Materialization combines free transformed parameters with fixed values. This
# is the bridge between optimizer vectors and structured EKF matrices.
@testset "Parameter materialization" begin
    layout = DataFrame(
        matrix = [:A, :A, :A, :A],
        row = [1, 2, 1, 2],
        col = [1, 1, 2, 2],
    )
    axis = retrieve_axes(layout)
    sp = ContinuousTimeSEM.EKFParameters(
        BitVector([true, false, true, false]),
        [true, false, true, false],
        [false, false, false, false],
        [false, false, false, false],
        Function[p -> p[1] + 1, p -> 2 * p[2]],
        Function[],
        Function[],
        [1, 2],
        axis,
        [false, true, false, true],
        AbstractFloat[10.0, 20.0],
    )

    @test sp.regular_transforms isa Tuple
    @test sp.fixed_values isa Vector{Float64}

    all_params = zeros(4)
    ContinuousTimeSEM._materialize_all_params!(all_params, [3.0, 5.0], sp)
    @test all_params == [4.0, 10.0, 10.0, 20.0]
    @test_throws DimensionMismatch ContinuousTimeSEM._materialize_all_params!(zeros(3), [3.0, 5.0], sp)

    f(x) = begin
        dual_params = Vector{typeof(x)}(undef, 4)
        ContinuousTimeSEM._materialize_all_params!(dual_params, [x, x], sp)
        sum(dual_params)
    end
    @test ForwardDiff.derivative(f, 3.0) ≈ 3.0
end

# Complex transforms depend on the current state and PARS block, so test both
# correct index assignment and mismatch validation.
@testset "Complex parameter transforms" begin
    all_params = zeros(4)
    indices = [2, 4]
    transforms = Function[
        (state, pars) -> state[1] + pars[1],
        (state, pars) -> state[2] * pars[2],
    ]

    ContinuousTimeSEM.apply_complex_transforms_at_indices!(
        all_params,
        indices,
        transforms,
        [1.0, 2.0],
        [10.0, 5.0],
    )

    @test all_params == [0.0, 11.0, 0.0, 10.0]

    tuple_params = zeros(4)
    ContinuousTimeSEM.apply_complex_transforms_at_indices!(
        tuple_params,
        indices,
        Tuple(transforms),
        [1.0, 2.0],
        [10.0, 5.0],
    )
    @test tuple_params == all_params

    @test_throws DimensionMismatch ContinuousTimeSEM.apply_complex_transforms_at_indices!(
        all_params,
        [1],
        transforms,
        [1.0, 2.0],
        [10.0, 5.0],
    )
end
