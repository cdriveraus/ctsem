"""
    ContinuousTimeSEM

Continuous-time structural equation modeling utilities in Julia, with helpers
for integration from R.
"""
module ContinuousTimeSEM

# Imports.
#
# Deliberately narrow, and deliberately not re-exported. This module used to
# `@reexport using ComponentArrays, DataFrames, Distributions, ForwardDiff,
# LinearAlgebra, Optim, StaticArrays`, which forced all of them to load and
# dumped their names into every caller's namespace. Distributions, StaticArrays
# and Reexport had no call sites at all; DataFrames had one, as a row container
# in the R interface, and is now a test-only dependency. The closure went from
# 111 packages to 56, a fresh install from 268 MB to 124 MB, and `using
# ContinuousTimeSEM` from ~8.7 s to ~3.8 s -- paid in every R session.
using ComponentArrays, ForwardDiff, LinearAlgebra, Optim

export hello
"""
    hello()

Print a short message and return `nothing`.

This is a lightweight smoke-test helper for checking that the Julia module is
loaded correctly from R.
"""
function hello()
    println("Hello from Julia")
    return nothing
end

export scalar_square
"""
    scalar_square(x)

Return `x^2`.

This small helper is useful as a scalar round-trip test when calling Julia from
R or other host environments.
"""
function scalar_square(x)
    return x^2
end

# Includes
include("small_linalg.jl")
include("parameters.jl")
include("constrain_cor_sqrt.jl")
include("r_interface.jl")
include("helper_functions.jl")
include("ksolve.jl")
include("parameter_transforms.jl")
include("buffered_exponential.jl")
include("buffered_lyapunov.jl")
include("discrete_time_form.jl")
include("log_likelihood.jl")
include("workspace_buffers.jl")
include("kalman_filters.jl")
include("benchmark_setup.jl")
include("ctsem_backend.jl")
include("adjoint_primitives.jl")
include("adjoint_parameters.jl")
include("reverse_scratch.jl")
include("adjoint_ekf.jl")
include("adjoint.jl")
include("summary_matrices.jl")
include("kalman_trace.jl")
include("laplace.jl")
include("quadrature.jl")

# Last, because it exercises everything above it.
include("precompile_workload.jl")


end # module ContinuousTimeSEM
