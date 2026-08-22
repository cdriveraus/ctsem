using LinearAlgebra
using Random

"""
    constraincorsqrt1(mat, epsilon=1e-5)

Transform an unconstrained lower-triangular correlation square-root parameter
matrix into a valid correlation Cholesky-like factor.

The returned matrix has rows scaled so that `o * o'` is positive definite, with
`epsilon` used as a small numerical margin.
"""
function constraincorsqrt1(mat::AbstractArray{T}, epsilon = 1e-5) where T
    d = size(mat, 1)
    o = zeros(T, d, d)
    ss = zeros(T, d)
    s = zeros(T, d)
    
    # Step 1: Compute sums of squares (ss) and sums (s) for each row
    # TODO: Rewrite this for performance
    for i in 1:d
        for j in 1:d
            if j > i
                ss[i] += mat[j, i]^2
                s[i] += mat[j, i]
            elseif j < i
                ss[i] += mat[i, j]^2
                s[i] += mat[i, j]
            end
        end
        s[i] += epsilon
        ss[i] += epsilon
    end
    
    # Step 2: Compute orthogonalized values for correlation matrix
    # TODO: Rewrite this for performance
    for i in 1:d
        o[i, i] = 0
        r1 = sqrt(ss[i])
        r3 = abs(s[i]) / r1 - 1
        r4 = sqrt(log1p(exp(2 * (abs(s[i]) - s[i] - 1) - 4)))
        r = (r4 * r3 + 1) * r4 + 1
        r = sqrt(ss[i] + r)
        
        for j in 1:d
            if j > i
                o[i, j] = mat[j, i] / r
            elseif j < i
                o[i, j] = mat[i, j] / r
            end
        end
        o[i, i] = sqrt(1 - sum(o[i, :].^2) + epsilon)
    end
    
    return o
end

"""
    sdcovsqrt2cov(mat, choleskymats)

Convert standard-deviation/correlation square-root parameters to a covariance.

The diagonal of `mat` supplies standard deviations, and the lower triangle
supplies unconstrained correlation parameters. `choleskymats` is currently
accepted for compatibility with the R-side interface.
"""
function sdcovsqrt2cov(mat, choleskymats)
    # TODO: Rewrite this for performance
    # if size(mat, 1) == 0
    #     return Symmetric(mat, :L) 
    # elseif choleskymats < 1

        # TODO: 
        diag_vals = Diagonal(mat)
        # corr_mat = constraincorsqrt1(mat)
        corr_mat = constraincorsqrt1_vec(Symmetric(mat, :L))
        return Symmetric(diag_vals * corr_mat * corr_mat' * diag_vals, :L)
    # else
    #     # TODO: Implement this as a specialization of the function
    #     return Symmetric(mat * mat', :L)
    # end
end

"""
    sdcovsqrt2cov!(buffer, mat, choleskymats)

In-place buffered version of `sdcovsqrt2cov`.

The covariance is written to `buffer.out`; other fields of `buffer` are used as
scratch storage. `choleskymats` is currently accepted for compatibility.
"""
function sdcovsqrt2cov!(buffer, mat, choleskymats)
    return sdcovsqrt2cov!(buffer, mat, choleskymats, Val(size(mat, 1)))
end

function sdcovsqrt2cov!(buffer, mat, choleskymats, dim::Val{d}) where {d}
    # TODO: Rewrite this for performance
    # if size(mat, 1) == 0
    #     # return Symmetric(mat, :L) 
    #     copyto!(buffer.out, mat)
    #     return nothing
    # elseif choleskymats < 1
        constraincorsqrt1_vec!(buffer, mat, 1e-5, dim)

        @inbounds for j in 1:d, i in 1:d
            buffer.intermediate[i, j] = mat[i, i] * buffer.out[i, j]
        end

        _mul_right_transpose!(buffer.out, buffer.intermediate, buffer.intermediate, dim, dim, dim)
        return nothing
    # else
    #     # TODO: Implement this as a specialization of the function
    #     # return Symmetric(mat * mat', :L)
    #     mul!(buffer.out, mat, mat')
    #     return nothing
    # end
end

"""
    is_positive_definite(mat)

Return whether all eigenvalues of `mat` are strictly positive.
"""
function is_positive_definite(mat::AbstractMatrix{<:AbstractFloat})
    return all(eigvals(mat) .> 0)
end

"""
    run_tests()

Run simple interactive checks for the correlation and covariance transforms.

This helper prints diagnostic output and throws an error if a generated
correlation matrix is not positive definite.
"""
function run_tests()
    Random.seed!(123)
    test_matrices = [
        randn(2, 2),
        randn(3, 3),
        randn(4, 4),
        diagm([1.0, 2.0, 3.0, 4.0]),  # Diagonal case
        zeros(4, 4)                   # Edge case
    ]
    # Ensure upper triangles are zero
    for mat in test_matrices
        mat[triu(ones(size(mat)), 1) .== 1] .= 0.0
    end

    for (i, test_mat) in enumerate(test_matrices)
        println("\nTest case $i:")
        println("Input matrix:")
        println(test_mat)

        # Constrain to correlation matrix
        constrained_mat = constraincorsqrt1(test_mat)
        corr_mat = constrained_mat * constrained_mat'

        println("Constrained Cholesky factor:")
        println(constrained_mat)
        println("Resulting correlation matrix:")
        println(corr_mat)

        if !is_positive_definite(corr_mat)
            error("Test case $i failed: Not positive definite.")
        else
            println("Test case $i passed: Positive definite.")
        end

        # Test sdcovsqrt2cov
        cov_mat = sdcovsqrt2cov(test_mat, 0)
        println("Resulting covariance matrix:")
        println(cov_mat)

        if !is_positive_definite(cov_mat)
            println("Covariance test case $i failed: Not positive definite.")
        else
            println("Covariance test case $i passed: Positive definite.")
        end
    end
end

# Run tests
# run_tests()

# Tests passed and the results equal the ones in R, so the functions are assumed to be correct

"""
    first_part(mat, epsilon=1e-5)

Compute row-wise off-diagonal sums of squares and sums for a symmetric matrix.

This is the reference loop implementation used by the correlation constraint
experiments in this file.
"""
function first_part(mat::AbstractArray{T}, epsilon = 1e-5) where T
    d = size(mat, 1)
    ss = zeros(T, d)
    s = zeros(T, d)
    
    # Step 1: Compute sums of squares (ss) and sums (s) for each row
    for i in 1:d
        for j in 1:d
            if j > i
                ss[i] += mat[j, i]^2
                s[i] += mat[j, i]
            elseif j < i
                ss[i] += mat[i, j]^2
                s[i] += mat[i, j]
            end
        end
        s[i] += epsilon
        ss[i] += epsilon
    end
    return ss, s
end

"""
    first_part_modified(mat, epsilon=1e-5)

Vectorized variant of `first_part` using `Symmetric(mat, :L)`.
"""
function first_part_modified(mat::AbstractArray{T}, epsilon = 1e-5) where T
    sym_mat = Symmetric(mat, :L)
    s = sum(sym_mat, dims=1)' .- diag(sym_mat) .+ epsilon
    ss = sum(abs2, sym_mat, dims=1)' .- diag(sym_mat).^2 .+ epsilon
    return ss, s
end

"""
    first_part_transposed(mat, epsilon=1e-5)

Vectorized variant of `first_part` based on a transposed upper-triangle view.
"""
function first_part_transposed(mat::AbstractArray{T}, epsilon = 1e-5) where T
    sym_mat = Symmetric(mat', :U)
    s = sum(sym_mat, dims=2) .- diag(sym_mat) .+ epsilon
    ss = sum(abs2, sym_mat, dims=2) .- diag(sym_mat).^2 .+ epsilon
    return ss, s
end

"""
    first_part_diag(mat, epsilon=1e-5)

Vectorized variant of `first_part` that materializes the diagonal once.
"""
function first_part_diag(mat::AbstractArray{T}, epsilon = 1e-5) where T
    sym_mat = Symmetric(mat, :L)
    mat_diag = collect(diag(sym_mat))
    s = sum(sym_mat, dims=2) .- mat_diag .+ epsilon
    ss = sum(abs2, sym_mat, dims=2) .- mat_diag.^2 .+ epsilon
    return ss, s
end

"""
    second_part(mat, ss, s, epsilon=1e-5)

Build the constrained correlation square-root factor from precomputed `ss` and
`s` vectors.

This is the reference loop implementation corresponding to `first_part`.
"""
function second_part(mat::AbstractArray{T}, ss::AbstractArray{T}, s::AbstractArray{T}, epsilon = 1e-5) where T
    d = length(ss)
    o = zeros(T, d, d)
    for i in 1:d
        r1 = sqrt(ss[i])
        r3 = abs(s[i]) / r1 - 1
        r4 = sqrt(log1p(exp(2 * (abs(s[i]) - s[i] - 1) - 4)))
        r = (r4 * r3 + 1) * r4 + 1
        r = sqrt(ss[i] + r)
        
        for j in 1:d
            if j > i
                o[i, j] = mat[j, i] / r
            elseif j < i
                o[i, j] = mat[i, j] / r
            end
        end
        o[i, i] = sqrt(1 - sum(o[i, :].^2) + epsilon)
    end
    
    return o
end

"""
    second_part_improved(mat, ss, s, epsilon=1e-5)

Vectorized scaling variant of `second_part`.
"""
function second_part_improved(mat::AbstractArray{T}, ss::AbstractArray{T}, s::AbstractArray{T}, epsilon = 1e-5) where T
    d = length(ss)
    o = zeros(T, d, d)
    r1 = sqrt.(ss)
    r3 = @. abs(s) / r1 - 1
    r4 = @. sqrt(log1p(exp(2 * (abs(s) - s - 1) - 4)))
    r = @. sqrt(ss + (r4 * r3 + 1) * r4 + 1)
    for i in 1:d
        for j in 1:d
            if j > i
                o[i, j] = mat[j, i] / r[i]
            elseif j < i
                o[i, j] = mat[i, j] / r[i]
            end
        end
        o[i, i] = sqrt(1 - sum(o[i, :].^2) + epsilon)
    end
    
    return o
end

"""
    second_part_improved_2(mat, ss, s, epsilon=1e-5)

Variant of `second_part` that reads off-diagonal entries through a symmetric
view of the lower triangle.
"""
function second_part_improved_2(mat::AbstractArray{T}, ss::AbstractArray{T}, s::AbstractArray{T}, epsilon = 1e-5) where T
    d = length(ss)
    o = zeros(T, d, d)
    r3 = @. abs(s) / sqrt(ss) - 1
    r4 = @. sqrt(log1p(exp(2 * (abs(s) - s - 1) - 4)))
    r = @. sqrt(ss + (r4 * r3 + 1) * r4 + 1)
    sym_mat = Symmetric(mat, :L)
    for i in 1:d
        for j in 1:d
            if j != i
                o[i, j] = sym_mat[i, j] / r[i]
            end
        end
        o[i, i] = sqrt(1 - sum(o[i, :].^2) + epsilon)
    end
    
    return o
end

"""
    second_part_improved_3(mat, ss, s, epsilon=1e-5)

Broadcast-oriented variant of `second_part` using a symmetric lower-triangle
view and row-wise scaling.
"""
function second_part_improved_3(mat::AbstractArray{T}, ss::AbstractArray{T}, s::AbstractArray{T}, epsilon = 1e-5) where T
    tmp = @. sqrt(log1p(exp(2 * (abs(s) - s - 1) - 4)))
    r = @. sqrt(ss + (tmp * (abs(s) / sqrt(ss) - 1) + 1) * tmp + 1)
    sym_mat = Symmetric(mat, :L)
    o = sym_mat ./ r
    o[diagind(o)] .= sqrt.(1 .- sum(x -> x^2, o, dims = 2) .+ diag(o).^2 .+ epsilon)
    return o
end

# inp = rand(4, 4)

# ss, s = first_part(inp)
# first_part_modified(inp)
# first_part(inp) .≈ first_part_modified(inp)
# first_part(inp) .≈ first_part_transposed(inp)
# first_part(inp) .≈ first_part_diag(inp)
# first_part_diag(inp)

# println("Original:")
# @benchmark first_part($inp)
# println("Modified:")
# @benchmark first_part_modified($inp)
# println("Transposed:")
# @benchmark first_part_transposed($inp)
# println("Diagonal:")
# @benchmark first_part_diag($inp)


# second_part(inp, first_part(inp)...) ≈ constraincorsqrt1(inp)
# second_part(inp, first_part(inp)...) ≈ second_part_improved(inp, ss, s)
# second_part(inp, first_part(inp)...) ≈ second_part_improved_2(inp, ss, s)
# second_part(inp, first_part(inp)...) ≈ second_part_improved_3(inp, ss, s)

# println("Original:")
# @benchmark second_part($inp, $ss, $s)
# println("Improved:")
# @benchmark second_part_improved($inp, $ss, $s)
# println("Improved 2:")
# @benchmark second_part_improved_2($inp, $ss, $s)
# println("Improved 3:")
# @benchmark second_part_improved_3($inp, $ss, $s)

"""
    constraincorsqrt1_vec(sym_mat, epsilon=1e-5)

Vectorized correlation square-root constraint for a symmetric matrix view.

The input is treated as symmetric, the off-diagonal entries are scaled row-wise,
and the diagonal is adjusted so rows have valid correlation-factor norms.
"""
function constraincorsqrt1_vec(sym_mat, epsilon = 1e-5)
    d = collect(diag(sym_mat))
    s = sum(sym_mat, dims=1)' .- d .+ epsilon
    ss = sum(abs2, sym_mat, dims=1)' .- d.^2 .+ epsilon
    tmp = @. sqrt(log1p(exp(2 * (abs(s) - s - 1) - 4)))
    r = @. sqrt(ss + (tmp * (abs(s) / sqrt(ss) - 1) + 1) * tmp + 1)
    o = sym_mat ./ r
	p = o - diagm(diag(o))
	dd = reshape(sqrt.(1 .- sum(abs2, p, dims = 2) .+ epsilon), :) 
    return p + diagm(dd)
end

# buffer requires:
# s = Vector{S}(undef, d)
# ss = Vector{S}(undef, d)
# r = Vector{S}(undef, d)
# row_sq = Vector{S}(undef, d)
# out = Matrix{S}(undef, d, d)

"""
    _sym_lower_get(mat, i, j)

Return `mat[i, j]` from a matrix whose lower triangle stores symmetric entries.
"""
@inline function _sym_lower_get(mat::AbstractMatrix, i::Int, j::Int)
    return i >= j ? mat[i, j] : mat[j, i]
end

"""
    constraincorsqrt1_vec!(buffer, sym_mat, epsilon=1e-5)

In-place buffered version of `constraincorsqrt1_vec`.

The constrained square-root factor is written to `buffer.out`; the vectors in
`buffer` hold intermediate row sums, scales, and row norms.
"""
function constraincorsqrt1_vec!(buffer, sym_mat, epsilon = 1e-5)
    return constraincorsqrt1_vec!(buffer, sym_mat, epsilon, Val(size(sym_mat, 1)))
end

function constraincorsqrt1_vec!(buffer, sym_mat, epsilon, ::Val{d}) where {d}
    S = eltype(sym_mat)
    e = convert(S, epsilon)

    @inbounds for i in 1:d
        si = e
        ssi = e
        for j in 1:d
            if j != i
                v = _sym_lower_get(sym_mat, i, j)
                si += v
                ssi += v * v
            end
        end
        buffer.s[i] = si
        buffer.ss[i] = ssi
    end

    @inbounds for i in 1:d
        si = buffer.s[i]
        ssi = buffer.ss[i]
        abs_si = abs(si)
        tmp = sqrt(log1p(exp(2 * (abs_si - si - one(S)) - 4)))
        buffer.r[i] = sqrt(ssi + (tmp * (abs_si / sqrt(ssi) - one(S)) + one(S)) * tmp + one(S))
    end

    @inbounds for i in 1:d
        inv_ri = inv(buffer.r[i])
        sq = zero(S)
        for j in 1:d
            if j == i
                buffer.out[i, j] = zero(S)
            else
                v = _sym_lower_get(sym_mat, i, j) * inv_ri
                buffer.out[i, j] = v
                sq += v * v
            end
        end
        buffer.row_sq[i] = sq
    end

    @inbounds for i in 1:d
        buffer.out[i, i] = sqrt(one(S) - buffer.row_sq[i] + e)
    end
end

# constraincorsqrt1(inp) ≈ constraincorsqrt1_vec(inp)
# constraincorsqrt1(inp)
# constraincorsqrt1_vec(inp)

# inp = rand(4, 4)
# println("Original:")
# @benchmark constraincorsqrt1($inp)
# println("Combined:")
# @benchmark constraincorsqrt1_vec($inp)
