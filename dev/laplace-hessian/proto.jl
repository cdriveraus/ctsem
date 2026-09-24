# Laplace Hessian by central differences of the exact gradient, with each
# differencing point's inner modes started from the base point's modes rather
# than from the origin. Prototype, loaded into Main next to the engine:
#   include("dev/laplace-hessian/proto.jl")
module LaplaceHess

import ContinuousTimeSEM as CT
using LinearAlgebra

# Every unit's mode reset to the base point's before an evaluation, so each
# differencing point starts from the same place whatever was evaluated before
# it: deterministic, and 1e-4 from where the modes were solved from the origin.
function _restore!(L, base)
    for U in eachindex(base)
        L.modes[U] = copy(base[U])
    end
end

function _grad(L, x)
    r = CT.ctsem_evaluate(L, x; gradient=true)
    (g=collect(Float64, r.gradient), its=sum(L.inner_iterations),
     converged=all(L.inner_converged))
end

"""Central-difference Hessian of the exact gradient, warm-started per point."""
function warm_hessian(L, values; step::Real=1e-4)
    x = collect(Float64, values)
    n = length(x)
    CT.ctsem_evaluate(L, x; gradient=false)         # modes at x, from the origin
    base = deepcopy(L.modes)
    prev = CT.ctsem_set_warm_start!(true)
    H = zeros(n, n); its = 0; converged = true
    try
        for j in 1:n
            h = step * max(1.0, abs(x[j]))
            xp = copy(x); xp[j] += h
            xm = copy(x); xm[j] -= h
            _restore!(L, base); p = _grad(L, xp)
            _restore!(L, base); m = _grad(L, xm)
            its += p.its + m.its
            converged &= p.converged && m.converged
            H[:, j] = (p.g .- m.g) ./ (2h)
        end
    finally
        CT.ctsem_set_warm_start!(prev)
        _restore!(L, base)
    end
    (H=(H .+ H') ./ 2, iterations=its, converged=converged)
end

"""
What one gradient costs at a point 1e-4 away from x, cold (from the origin) and
warm (from x's modes), minimum of `reps`; plus the mean inner iterations.
"""
function point_costs(L, values, reps::Integer=3)
    x = collect(Float64, values)
    CT.ctsem_evaluate(L, x; gradient=false)
    base = deepcopy(L.modes)
    y = copy(x); y[1] += 1e-4 * max(1.0, abs(x[1]))
    nunits = length(base)
    tcold = Inf; tval = Inf; twarm = Inf; icold = 0; iwarm = 0
    prev = CT.ctsem_set_warm_start!(false)
    try
        for _ in 1:reps
            t = @elapsed CT.ctsem_evaluate(L, y; gradient=false); tval = min(tval, t)
            t = @elapsed r = _grad(L, y); tcold = min(tcold, t); icold = r.its
        end
        CT.ctsem_set_warm_start!(true)
        for _ in 1:reps
            _restore!(L, base)
            t = @elapsed r = _grad(L, y); twarm = min(twarm, t); iwarm = r.its
        end
    finally
        CT.ctsem_set_warm_start!(prev)
        _restore!(L, base)
    end
    (value=tval, cold=tcold, warm=twarm, cold_iterations=icold / nunits,
     warm_iterations=iwarm / nunits)
end

end # module
