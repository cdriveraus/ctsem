# Theta's triangular map (sample_map.jl): its inverse, its pullback, its fit's
# choice of parents and size, and its weights, on closed forms.
#
# A wrong map does not error: the sampler is exact whatever the map is, so a
# mistake in the map shows only as slower sampling. These check the map itself.

using LinearAlgebra, Random, Statistics

const _CT_TRI = ContinuousTimeSEM

# A map with random coefficients on random parents, products included.
function _random_tri(n, rng)
    par = [[j for j in 1:(k - 1) if rand(rng) < 0.6] for k in 1:n]
    coef(sc) = [sc .* randn(rng, length(par[k])) for k in 1:n]
    pairs(sc) = [sc .* randn(rng, length(par[k]) * (length(par[k]) - 1) ÷ 2) for k in 1:n]
    return _CT_TRI._JMapTri(n, par, 0.1 .* randn(rng, n), coef(0.2), coef(0.1), pairs(0.1),
        0.1 .* randn(rng, n), coef(0.2), 0.2 .* randn(rng, n), -0.2 .* rand(rng, n))
end

@testset "the triangular map inverts, and its pullback is its gradient" begin
    rng = Xoshiro(5)
    n = 7
    tri = _random_tri(n, rng)
    w = randn(rng, n)
    u, lj = _CT_TRI._jmap_tri_forward(tri, w)
    @test maximum(abs, _CT_TRI._jmap_tri_inverse(tri, u) .- w) < 1e-12
    # The pullback of gu is the gradient in w of gu'u(w) + log|du/dw|.
    gu = randn(rng, n)
    h(x) = (v = _CT_TRI._jmap_tri_forward(tri, x); dot(gu, v[1]) + v[2])
    fd = [(h(w .+ 1e-6 .* (1:n .== j)) - h(w .- 1e-6 .* (1:n .== j))) / 2e-6 for j in 1:n]
    @test maximum(abs, _CT_TRI._jmap_tri_pullback(tri, w, gu) .- fd) < 1e-6
    # The log Jacobian is the log determinant of du/dw (triangular, so the
    # product of its diagonal), by differences.
    J = reduce(hcat, [(_CT_TRI._jmap_tri_forward(tri, w .+ 1e-6 .* (1:n .== j))[1] .-
        _CT_TRI._jmap_tri_forward(tri, w .- 1e-6 .* (1:n .== j))[1]) ./ 2e-6 for j in 1:n])
    @test abs(logabsdet(J)[1] - lj) < 1e-6
end

@testset "the fit reads the parents a shape needs, as many as the draws support" begin
    rng = Xoshiro(11)
    N = 4000
    # u2 bends with u1; u3 moves with the product u1 u2; u4 is independent.
    u1 = randn(rng, N)
    u2 = 0.8 .* (u1 .^ 2 .- 1) .+ 0.4 .* randn(rng, N)
    u3 = 0.6 .* u1 .* u2 .+ 0.3 .* randn(rng, N)
    u4 = randn(rng, N)
    U = permutedims(hcat(u1, u2, u3, u4))
    U = (U .- mean(U; dims=2)) ./ std(U; dims=2)
    tri = _CT_TRI._jmap_tri_fit(U)
    @test 1 in tri.parents[2]
    @test all(in(tri.parents[3]), [1, 2])
    # The product coefficient of u3's location on (u1, u2) carries the shape.
    P3 = tri.parents[3]
    j = findfirst(==((1, 2)), [(P3[q], P3[r]) for q in eachindex(P3) for r in (q + 1):length(P3)])
    @test j !== nothing && abs(tri.f3[3][j]) > 0.1
    # Held out, the map's latent draws are close to standard normal, and its
    # log density beats the whitening's by a margin.
    v1 = randn(rng, 2000)
    v2 = 0.8 .* (v1 .^ 2 .- 1) .+ 0.4 .* randn(rng, 2000)
    V = permutedims(hcat(v1, v2, 0.6 .* v1 .* v2 .+ 0.3 .* randn(rng, 2000), randn(rng, 2000)))
    V = (V .- mean(permutedims(hcat(u1, u2, u3, u4)); dims=2)) ./
        std(permutedims(hcat(u1, u2, u3, u4)); dims=2)
    W = reduce(hcat, [_CT_TRI._jmap_tri_inverse(tri, view(V, :, t)) for t in axes(V, 2)])
    @test maximum(abs, mean(W; dims=2)) < 0.1
    @test maximum(abs, std(W; dims=2) .- 1) < 0.1
    logq(x) = (w = _CT_TRI._jmap_tri_inverse(tri, x);
        -dot(w, w) / 2 - _CT_TRI._jmap_tri_forward(tri, w)[2])
    gain = mean(logq(view(V, :, t)) + dot(view(V, :, t), view(V, :, t)) / 2 for t in axes(V, 2))
    @test gain > 0.5
    # Too few draws for any coefficient: no map.
    @test _CT_TRI._jmap_tri_fit(U[:, 1:4]) === nothing
end

@testset "a weighted component's likelihood is its draws repeated by their weights" begin
    rng = Xoshiro(2)
    N, m = 40, 2
    x = randn(rng, N)
    S = randn(rng, m + 1, N)
    Phi, Psi = _CT_TRI._jmap_tri_designs(S, [1, 2])
    reps = rand(rng, 1:3, N)
    p = 0.2 .* randn(rng, size(Phi, 2) + size(Psi, 2) + 2)
    idx = reduce(vcat, [fill(i, reps[i]) for i in 1:N])
    rep = _CT_TRI._JMapTriComponent(x[idx], Phi[idx, :], Psi[idx, :], false)
    wtd = _CT_TRI._JMapTriComponent(x, Phi, Psi, false, Float64.(reps))
    @test _CT_TRI._jmap_tri_component_nll(wtd, p) ≈ _CT_TRI._jmap_tri_component_nll(rep, p)
    gr, Hr = _CT_TRI._jmap_tri_component_derivatives(rep, p)
    gw, Hw = _CT_TRI._jmap_tri_component_derivatives(wtd, p)
    @test gw ≈ gr
    @test Hw ≈ Hr
end

@testset "a size's fit laid out for more parents is the same component" begin
    rng = Xoshiro(4)
    N = 60
    S = randn(rng, 6, N)
    x = randn(rng, N)
    P, Q = [2, 5], [1, 2, 4, 5]
    Phi, Psi = _CT_TRI._jmap_tri_designs(S, P)
    PhiQ, PsiQ = _CT_TRI._jmap_tri_designs(S, Q)
    p = 0.3 .* randn(rng, size(Phi, 2) + size(Psi, 2) + 2)
    q = _CT_TRI._jmap_tri_embed(p, P, Q)
    @test length(q) == size(PhiQ, 2) + size(PsiQ, 2) + 2
    @test _CT_TRI._jmap_tri_component_nll(_CT_TRI._JMapTriComponent(x, PhiQ, PsiQ, false), q) ≈
        _CT_TRI._jmap_tri_component_nll(_CT_TRI._JMapTriComponent(x, Phi, Psi, false), p)
    @test _CT_TRI._jmap_tri_embed(p, P, [1, 2, 4]) === nothing
end
