using ComponentArrays
using DataFrames

# Transform strings are parsed from the R-side parameter table; this checks that
# the generated lambdas have the expected signatures and values.
@testset "Transform string helpers" begin
    tf = eval(Meta.parse(ContinuousTimeSEM.generate_transform_string("param[1] + 2")))
    ctf = eval(Meta.parse(ContinuousTimeSEM.generate_complex_transform_string("state[1] + PARS[1]")))

    @test ContinuousTimeSEM.log1p_exp(0.0) ≈ log(2.0)
    @test tf([3.0]) == 5.0
    context = ContinuousTimeSEM.CTSEMRowContext([1.0], (PARS=[4.0],),
        Float64[], Float64[], 0.0, 0.0, 1, 1)
    @test ctf(context) == 5.0
end

# Fixed-value mapping fills the non-optimized positions in the flattened
# parameter vector, so both mutating APIs should preserve indexing.
@testset "Fixed-value mapping" begin
    target = zeros(4)
    positions = [true, false, true, false]

    @test ContinuousTimeSEM.map_fixed_values!(target, positions, [10.0, 20.0]) === nothing
    @test target == [10.0, 0.0, 20.0, 0.0]

    returned = ContinuousTimeSEM.map_fixed_values(target, [false, true, false, true], [30.0, 40.0])
    @test returned === target
    @test target == [10.0, 30.0, 20.0, 40.0]

    expanded = ContinuousTimeSEM.params_to_all_params([1.0, 2.0], BitVector([true, false, true]))
    @test expanded[[1, 3]] == [1.0, 2.0]
end

# Component axes determine how the flat parameter vector is viewed as named
# matrices; incorrect axes would corrupt every downstream parameter lookup.
@testset "Component axes from data frames" begin
    layout = DataFrame(
        matrix = [:A, :A, :A, :A, :B, :B],
        row = [1, 2, 1, 2, 1, 1],
        col = [1, 1, 2, 2, 1, 2],
    )

    names, nrows, ncols = retrieve_names_and_dims(layout)
    @test names == [:A, :B]
    @test nrows == [2, 1]
    @test ncols == [2, 2]

    axis = retrieve_axes(layout)
    values = ComponentVector(zeros(6), axis)
    @test size(values.A) == (2, 2)
    @test size(values.B) == (1, 2)
end

# This exercises the full DataFrame-to-EKFParameters parser, including fixed
# values, mutable masks, and regular/predict/update transforms.
@testset "EKF parameter parsing from data frame" begin
    df = DataFrame(
        matrix = [:A, :A, :A, :A],
        row = [1, 2, 1, 2],
        col = [1, 1, 2, 2],
        parnumber = Union{Missing, Int}[1, missing, 2, missing],
        value = Union{Missing, Float64}[missing, 4.0, missing, 8.0],
        transform = Union{Missing, String}["param[1]", missing, "param[2]^2", missing],
        predicttransform = Union{Missing, String}[missing, "state[1] + PARS[1]", missing, missing],
        updatetransform = Union{Missing, String}[missing, missing, missing, "state[1] - PARS[1]"],
    )

    sp = ekf_from_data_frame(df)

    @test sp.regular_transforms isa Tuple
    @test sp.predict_transforms isa Tuple
    @test sp.update_transforms isa Tuple
    @test sp.fixed_values isa Vector{Float64}
    @test sp.mutables == BitVector([true, false, true, false])
    @test sp.fixed_indices == [false, true, false, true]
    @test sp.fixed_values == [4.0, 8.0]
    @test sp.diffusion_state_indices == Int[]
    @test findall(sp.predict_transforms_indices) == [2]
    @test findall(sp.update_transforms_indices) == [4]
    @test sp.regular_transforms[1]([2.0, 3.0]) == 2.0
    @test sp.regular_transforms[2]([2.0, 3.0]) == 9.0
    context = ContinuousTimeSEM.CTSEMRowContext([2.0], (PARS=[5.0],),
        Float64[], Float64[], 0.0, 0.0, 1, 1)
    @test sp.predict_transforms[1](context) == 7.0
    @test sp.update_transforms[1](context) == -3.0

    indexed = ekf_from_data_frame(df, DataFrame(
        parameter=Int[], predictor=Int[], coefficient=Int[]), [1])
    @test indexed.diffusion_state_indices == [1]
end
