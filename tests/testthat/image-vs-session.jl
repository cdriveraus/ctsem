# For test-julia-precompile.R: the engine's value, gradient and exact Hessian
# at the points the package image's own workload was built from, one line
# each, full precision.
#
#   julia image-vs-session.jl <engine environment> <outfile>
#
# Run once against the engine loaded from its package image and once against a
# workload-free build of the same source, whose code is therefore compiled in
# the session; the test compares the two files. It rebuilds each captured
# model's objective from the recorded calls (precompile_shapes.jl), so it needs
# no R and no fit.
# What the R bridge's Julia server (JuliaConnectoR's main.jl) loads before the
# engine, in its order and from the default environment: the package image is
# valid only in a process that loaded the same packages first, since it was
# built in one (Tables brings OrderedCollections from the default environment,
# and a process without it rejects the image and builds another).
import Pkg
try
    @eval import Tables
catch
end
using InteractiveUtils
import REPL
Pkg.activate(ARGS[1]; io=devnull)
using ContinuousTimeSEM
const CT = ContinuousTimeSEM
const BUILD = (:ctsem_transforms_cached, :ekf_from_columns, :ctsem_objective,
    :ctsem_laplace_objective)

function replay_points(out::AbstractString)
    open(out, "w") do io
        println(io, "workload ", CT._PRECOMPILE_WORKLOAD_ENABLED)
        for name in sort(collect(keys(CT._PRECOMPILE_SHAPES)))
            calls = CT._PRECOMPILE_SHAPES[name]
            results = Vector{Any}(nothing, length(calls))
            seen = 0
            for (k, (fn, args, kw)) in enumerate(calls)
                if fn in BUILD
                    a = Any[CT._precompile_value(x, results) for x in args]
                    kws = Dict{Symbol,Any}()
                    for (key, value) in kw
                        kws[key] = CT._precompile_value(value, results)
                    end
                    results[k] = Base.invokelatest(getfield(CT, fn), a...; kws...)
                elseif fn === :ctsem_evaluate && seen < 2
                    objective = CT._precompile_value(args[1], results)
                    x = Vector{Float64}(args[2])
                    r = Base.invokelatest(CT.ctsem_evaluate, objective, x;
                        gradient=true, contributions=false, gradient_method="adjoint")
                    println(io, name, " value ", k, " ", repr(Float64(r.value)))
                    println(io, name, " gradient ", k, " ", join(repr.(Float64.(r.gradient)), " "))
                    if seen == 0 && objective isa CT.CTSEMObjective
                        H = Base.invokelatest(CT.ctsem_hessian, objective, x)
                        println(io, name, " hessian ", k, " ", join(repr.(vec(H)), " "))
                    end
                    seen += 1
                end
            end
        end
        println(io, "done")
    end
end

replay_points(ARGS[2])
