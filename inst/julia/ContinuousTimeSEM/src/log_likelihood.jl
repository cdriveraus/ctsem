using LinearAlgebra

"""
    _kalman_loglikelihood_cholesky!(ll_buffer, S, ỹ, log2π_const)

Return the Gaussian Kalman innovation log-likelihood.

`S` must be a Cholesky factorization of the innovation covariance. The vector
`ll_buffer` is overwritten with `S \\ ỹ` and reused to avoid allocations.
"""
function _kalman_loglikelihood_cholesky!(ll_buffer::AbstractVector, S, ỹ, log2π_const::Real)
    # Deliberately *not* `Val(length(...))`. The observed-variable count
    # varies row to row with the missingness pattern, so a `Val` built from
    # it is a runtime value: it forces a dynamic dispatch (and an
    # allocation) on every row, and specializes the callee afresh for every
    # distinct count. The loop below is over the Cholesky diagonal, so a
    # compile-time bound buys nothing here anyway.
    return _kalman_loglikelihood_cholesky!(ll_buffer, S, ỹ, log2π_const, length((ỹ)))
end

"""Compile-time-length form, kept for callers that already hold a `Val`."""
function _kalman_loglikelihood_cholesky!(ll_buffer::AbstractVector, S, ỹ, log2π_const::Real, ::Val{d}) where {d}
    return _kalman_loglikelihood_cholesky!(ll_buffer, S, ỹ, log2π_const, d)
end

function _kalman_loglikelihood_cholesky!(ll_buffer::AbstractVector, S, ỹ, log2π_const::Real, d::Int)
    # 1) Solve S * z = ỹ using the in-place Cholesky factorization of S.
    #    ll_buffer stores z = S \ ỹ.
    ldiv!(ll_buffer, S, ỹ)

    # 2) For S = U' * U (Cholesky), logdet(S) = 2 * sum(log(diag(U))).
    U = S.factors
    logdet_half = zero(eltype(ll_buffer))
    @inbounds for idx in 1:d
        logdet_half += log(U[idx, idx])
    end

    # 3) Multivariate normal log-likelihood:
    #    -0.5 * (n*log(2π) + logdet(S) + ỹ' * (S \ ỹ)).
    return -0.5 * (d * log2π_const + 2 * logdet_half + dot(ỹ, ll_buffer))
end
