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
#
# The first line is the exception, and nothing here uses it. It loads what the
# R bridge's Julia server loads before this module (JuliaConnectoR 1.1.6,
# main.jl: Pkg, then REPL and InteractiveUtils), so that the package image is
# compiled against the method table it will be loaded into. Without it the
# precompile workload bought nothing in any R session. LibGit2, under Pkg,
# defines `cconvert(::Type{Ptr{StrArrayStruct}}, ::Vector)`; the image reaches
# `pointer(::Vector)` through the abstract instances inference makes at the
# engine's function barriers (a workspace from an untyped cache, a boxed
# capture); and loaded into a session that already had LibGit2, 4787 of the
# image's compiled methods failed verification through that one edge and were
# compiled again on first use -- 20 s of the workload's own calls, against
# 0.06 s in a bare session. It costs nothing in R, where these are loaded
# already, about 0.3 s in a bare Julia session, and no download: all three are
# standard libraries. test-julia-precompile.R fails if a first fit through
# the bridge starts compiling again.
import Pkg, REPL, InteractiveUtils
using ComponentArrays, ForwardDiff, LinearAlgebra, Optim, SpecialFunctions

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
# First, so that `using Printf` is in scope for every file that reports.
include("progress.jl")
# Checked from progress.jl's per-iteration hooks.
include("interrupt.jl")
include("small_linalg.jl")
include("opcounts.jl")
include("parameters.jl")
include("cov_cache.jl")
include("constrain_cor_sqrt.jl")
include("r_interface.jl")
# Socket options for the R bridge. No engine code depends on it, so it can sit
# anywhere; here it is beside the rest of the R-facing surface.
include("bridge_tuning.jl")
include("helper_functions.jl")
include("ksolve.jl")
include("parameter_transforms.jl")
include("buffered_exponential.jl")
include("frechet_exponential.jl")
include("cov_expm.jl")
include("buffered_lyapunov.jl")
include("series_discretization.jl")
include("discrete_time_form.jl")
include("log_likelihood.jl")
include("workspace_buffers.jl")
include("kalman_filters.jl")
include("ctsem_backend.jl")
include("substep_mesh.jl")
include("adjoint_primitives.jl")
include("adjoint_parameters.jl")
include("reverse_scratch.jl")
# Before adjoint_ekf.jl: the tape holds a vector of binary records, so the
# record type must exist when `CTSEMAdjointTape` is defined.
include("adjoint_binary.jl")
include("adjoint_ekf.jl")
include("adjoint.jl")
include("summary_matrices.jl")
include("kalman_trace.jl")
include("laplace.jl")
# After kalman_trace.jl and laplace.jl: it reads the filter trace on one route
# and the unit curvature on the other.
include("effect_information.jl")
# After laplace.jl: the batch subsets both the marginal and the laplace objective.
include("optimiser.jl")
include("quadrature.jl")
# After quadrature.jl: the continuation places its fixed rules with the leaf
# rule and the Gauss-Hermite grids defined there, and drives them with
# `ctsem_optimize` from ctsem_backend.jl and optimiser.jl.
include("laplace_continuation.jl")
# After quadrature.jl: the binary measurement update integrates the
# observation with the Gauss-Hermite rule defined there.
include("binary_measurement.jl")
# After binary_measurement.jl: the conditional observation model it defines --
# the log likelihood of one observation at a *known* linear predictor -- is
# exactly what a sampled state supplies, so the state path evaluates the same
# kernels the filter integrates.
include("state_sampling.jl")
include("particle_filter.jl")
# After particle_filter.jl, which brings in Random for the effect draws.
include("laplace_loo.jl")
include("sample_density.jl")
include("sample_nuts.jl")
include("sample_adapt.jl")
include("sample_run.jl")

# Last, because it exercises everything above it.
include("precompile_workload.jl")


# The operation counters in opcounts.jl increment during the precompile
# workload, and a `Ref` inside a `const` keeps whatever value it had when the
# package image was written. Zero them so a session starts counting from its
# own first evaluation, not from the workload's -- and the other diagnostic
# counters the workload moves the same way: the covariance- and expm-cache
# counts, the seeded-assembly refusal code and the Laplace gradient fallbacks.
# The OpenBLAS count is read here for the same reason: it belongs to this
# session, not to the process that wrote the image, and LinearAlgebra's own
# `__init__` has set it by now.
function __init__()
    ctsem_reset_opcounts!()
    ctsem_cov_cache_reset_counts()
    ctsem_cov_expm_reset_counts()
    ctsem_laplace_diag_reset!()
    _CTSEM_LAPLACE_FALLBACKS[] = 0
    _CTSEM_BLAS_START[] = LinearAlgebra.BLAS.get_num_threads()
    return nothing
end

end # module ContinuousTimeSEM
