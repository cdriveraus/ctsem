"""
    CTSEMReverseScratch(T, n, m, k)

Working storage for the reverse pass, one set per workspace.

The reverse pass is a few dozen small matrix products per row, and written the
obvious way each one allocates its result: on a 24-row, one-latent model that
is 5,454 heap objects per subject sweep, of which `_reverse_update!` and
`_reverse_predict!` alone are 71%. Serially that costs perhaps ten per cent.
Across threads it costs everything, because a Julia process has one heap and
allocation is where the threads queue -- 23 threads on that model run at 0.26x
of one, and 23 *processes* at 22.8x, which is how the heap was identified as
the wall.

So every temporary comes from here instead. The buffers are sized once, at the
full state and manifest dimensions, and a row with missing data takes a view of
the leading block rather than a fresh array. Names match the derivation in
`_reverse_update!` and `_reverse_predict!` one for one, so the buffered code
reads the same as the algebra it replaces.
"""
struct CTSEMReverseScratch{T}
    # --- measurement update, sized (n, m) --------------------------------
    Pr::Matrix{T}
    PHt::Matrix{T}
    S::Matrix{T}
    Sinv::Matrix{T}
    G::Matrix{T}
    M::Matrix{T}
    Sbar::Matrix{T}
    Sbar0::Matrix{T}
    Rbar::Matrix{T}
    Mbar::Matrix{T}
    Pnew::Matrix{T}
    Ps::Matrix{T}
    Gbar::Matrix{T}
    Hbar::Matrix{T}
    Lbar::Matrix{T}
    PHtbar::Matrix{T}
    ytilde::Vector{T}
    ybar::Vector{T}
    alpha::Vector{T}
    alphabar::Vector{T}
    beta::Vector{T}
    xnew::Vector{T}
    # generic products, by shape
    nn1::Matrix{T}
    nn2::Matrix{T}
    nm1::Matrix{T}
    mm1::Matrix{T}
    mm2::Matrix{T}
    # --- prediction substep, sized (n, k) --------------------------------
    Abar::Matrix{T}
    Pbar_new::Matrix{T}
    xbar_new::Vector{T}
    kk1::Matrix{T}
    kk2::Matrix{T}
    kk3::Matrix{T}
    nn3::Matrix{T}
    nn4::Matrix{T}
    nn5::Matrix{T}
    nn6::Matrix{T}
    kk4::Matrix{T}
    kk5::Matrix{T}
    kk6::Matrix{T}
    kv1::Vector{T}
    kv2::Vector{T}
    nv1::Vector{T}
    piv::Vector{Int}
end

function CTSEMReverseScratch(::Type{T}, n::Int, m::Int, k::Int) where {T}
    z(r, c) = zeros(T, r, c)
    v(r) = zeros(T, r)
    return CTSEMReverseScratch{T}(
        z(n, n), z(n, m), z(m, m), z(m, m), z(n, m), z(n, n),
        z(m, m), z(m, m), z(m, m), z(n, n), z(n, n), z(n, n),
        z(n, m), z(m, n), z(m, n), z(n, m),
        v(m), v(m), v(m), v(m), v(m), v(n),
        z(n, n), z(n, n), z(n, m), z(m, m), z(m, m),
        z(n, n), z(n, n), v(n), z(k, k), z(k, k), z(k, k), z(n, n), z(n, n),
        z(n, n), z(n, n), z(k, k), z(k, k), z(k, k),
        v(k), v(k), v(n), zeros(Int, max(n, k)))
end

"""The leading `r x c` block of a scratch buffer, as a view."""
@inline _rs(A::Matrix, r::Int, c::Int) = view(A, 1:r, 1:c)
"""The leading `r` entries of a scratch vector, as a view."""
@inline _rs(a::Vector, r::Int) = view(a, 1:r)

"""`destination .= (A + A') / 2`, in place."""
@inline function _symmetrize_into!(destination, A)
    n = size(A, 1)
    @inbounds for j in 1:n, i in 1:n
        destination[i, j] = (A[i, j] + A[j, i]) / 2
    end
    return destination
end
