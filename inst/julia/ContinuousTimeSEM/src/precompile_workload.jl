# Compile the whole default fit for a few common model shapes at package build
# time, rather than in the user's first fit.
#
# ## What it replays
#
# `precompile_shapes.jl` holds, for each model listed in
# `tools/generate-precompile-shapes.R`, every engine call a default `ctFit` made,
# captured through the R bridge with the arguments Julia actually received. The
# workload replays them in order. That is the only way to cover what a fit
# compiles: the code is specialised on the model's type and on each call's exact
# argument types, keywords included, and a fit reaches the Newton finish, the
# certification Hessian, the smoother, the effect check and on the Laplace route
# the continuation -- none of which a hand-written "evaluate and take two
# optimiser steps" workload reached. Hand-writing the model was worse still: the
# first version of this file got four details of it wrong and compiled a shape
# no fit asks for.
#
# ## Why a model type can be shared at all
#
# A model's type carries its transform closures, and `_regular_transform_closure`
# lifts the parameter index out of each, so a closure's type depends on its
# template rather than on which parameter it reads. A model is therefore typed by
# its matrix dimensions and the set of templates it uses, not by its layout, and
# one captured model covers every model of the same dimensions and templates.
# Anything else still compiles on first use.
#
# ## What it costs
#
# Build time, once per engine version: about a first fit's compile per model.
# And load time in every session, which is easy to miss: the image holds each
# model's specialisations whether or not the session fits that model, and loading
# them is most of what `using ContinuousTimeSEM` costs. So the list stays short.
#
# ## Failure must not be fatal
#
# A replay that fails degrades to "that model was not precompiled", never to a
# package that will not load -- hence the `try`, and a warning, because the only
# other symptom is a slow first fit. `JULIA_CTSEM_PRECOMPILE=0` skips it all, and
# `CTSEM_PRECOMPILE_WORKLOAD=false` keeps the specs but skips the replay: worth
# setting while editing the engine, since the cache key is a content hash and
# every edit pays the build again. Both are read at precompile time and are not
# part of the cache key: an image built with them set is reused, workload-free,
# by every later session of the same engine source. Rerunning the generator
# with them set, and getting a file identical to the one it loaded, leaves
# exactly that behind; change the source or clear the cache before measuring.

"""
A recorded argument that is the object an earlier recorded call returned: the
`call`-th entry of the same replay.
"""
struct _PrecompileRef
    call::Int
end

include("precompile_shapes.jl")

"""
The `EKFParameters` types the captured models produced.

Recorded because "does precompilation still cover real models?" has no other
symptom: when what R sends drifts, the captured types silently stop matching
and first fits are merely slow again. `ctsem_shape_is_precompiled` reads this.
"""
const _PRECOMPILE_SHAPE_TYPES = Set{Any}()

_precompile_value(x, results) = x isa _PrecompileRef ? results[x.call] : x

"""
    _precompile_replay(calls; build_only=false)

Make the recorded calls in order, each through `Base.invokelatest` so that a
transform closure `eval`ed by an earlier call is callable by a later one, and
each with its keywords in a `Dict` built in the recorded order -- which is how
the bridge passes them, so the keyword `NamedTuple` types match too.

`build_only` makes only the `ekf_from_columns` calls, for the top-level pass
that builds the specs.
"""
function _precompile_replay(calls; build_only::Bool=false)
    results = Vector{Any}(nothing, length(calls))
    for (k, (fn, args, kw)) in enumerate(calls)
        build_only && fn !== :ekf_from_columns && continue
        f = getfield(@__MODULE__, fn)
        a = Any[_precompile_value(x, results) for x in args]
        kws = Dict{Symbol,Any}()
        for (name, value) in kw
            kws[name] = _precompile_value(value, results)
        end
        results[k] = Base.invokelatest(f, a...; kws...)
    end
    return results
end

export ctsem_shape_is_precompiled
"""
    ctsem_shape_is_precompiled(matrix, row, col, parnumber, value, transform,
                               predicttransform, updatetransform, tdtransform)

Does a model built from these columns have one of the types compiled into the
package image?

This is the only observable difference between precompilation working and
precompilation quietly covering nothing: both load fine, and only one of them
makes the first fit fast. The R suite asks this about a real model so that a
drift in what `ctModelWriter` emits fails a test rather than a stopwatch. It
compares the spec type only, so it is necessary for a fast first fit and not
sufficient; the R suite also times one.
"""
function ctsem_shape_is_precompiled(matrix, row, col, parnumber, value, transform,
    predicttransform, updatetransform, tdtransform)
    spec = ekf_from_columns(matrix, row, col, parnumber, value, transform,
        predicttransform, updatetransform, tdtransform)
    return typeof(spec) in _PRECOMPILE_SHAPE_TYPES
end

using PrecompileTools: @compile_workload

# The specs first, in their own top-level statement: `ekf_from_columns` `eval`s
# the transform closures, and building them here fills the transform cache in
# an earlier world than the replay below. (The replay's `invokelatest` would
# cope anyway; this is also where the captured types are recorded.)
const _PRECOMPILE_BUILT = let built = Symbol[]
    if get(ENV, "JULIA_CTSEM_PRECOMPILE", "1") != "0"
        for name in keys(_PRECOMPILE_SHAPES)
            try
                calls = _PRECOMPILE_SHAPES[name]
                results = _precompile_replay(calls; build_only=true)
                for (k, call) in enumerate(calls)
                    call[1] === :ekf_from_columns &&
                        push!(_PRECOMPILE_SHAPE_TYPES, typeof(results[k]))
                end
                push!(built, name)
            catch err
                err isa InterruptException && rethrow()
                @warn "ContinuousTimeSEM: precompile model skipped" name exception = err
            end
        end
    end
    built
end

# A session run through the R bridge has Pkg loaded before this module, and
# the image is only usable there because this module loads it too: see the
# note at `import Pkg` in ContinuousTimeSEM.jl. Until that was found, this
# workload made every session load ~450 MB of specialisations and then
# compiled them again in the first fit.
if get(ENV, "CTSEM_PRECOMPILE_WORKLOAD", "true") != "false"
@compile_workload begin
    for name in _PRECOMPILE_BUILT
        try
            _precompile_replay(_PRECOMPILE_SHAPES[name])
        catch err
            err isa InterruptException && rethrow()
            @warn "ContinuousTimeSEM: precompile replay failed; a first fit of this model will compile" name exception = err
        end
    end
end
end
