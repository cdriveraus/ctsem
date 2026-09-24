# Julia engine unit suite against a given tree, in a persistent test environment.
#   julia julia-unit.jl <engine dir> <env dir>
using Pkg
engine, env = ARGS[1], ARGS[2]
Pkg.activate(env)
Pkg.develop(path = engine)
Pkg.add(["DataFrames", "ForwardDiff", "ChainRulesCore", "Optim", "LineSearches",
         "ComponentArrays", "FiniteDiff"])
cd(engine)
println("LOADED"); flush(stdout)
include(joinpath(engine, "test", "runtests.jl"))
println("DONE"); flush(stdout)
