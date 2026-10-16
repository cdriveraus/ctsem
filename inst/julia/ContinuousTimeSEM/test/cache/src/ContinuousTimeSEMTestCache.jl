"""
A compile cache for the engine's tests, and nothing else.

Most of the engine suite's time is compiling, not running: test_laplace.jl
runs in seconds once compiled and takes half an hour to compile, each
ForwardDiff oracle compiling the filter and the Laplace unit term for a dual
type of its own. This package runs the test files as its precompile workload,
so its package image keeps the engine code compiled at those dual types, and a
session that loads it before the tests finds that code already compiled.

Only for the tests. No ctsem user loads it: it lives under the engine's
`test/`, which neither the engine copy ctsem installs nor its version hash
covers, and the engine's own image is unchanged by it.

What it can keep is code whose types outlive the session that compiled it:
an engine function, or an engine object, differentiated at a dual type --
`_laplace_nested_gradient`, `ForwardDiff.gradient(objective, x)` on a
`CTSEMObjective`. A closure written in a test file is a new type every time
the file runs, so ForwardDiff over one compiles afresh whatever is cached.

The image is rebuilt whenever the engine's source changes, as every package
that depends on the engine is, so it pays off when the engine is unchanged
between building it and running the tests: running them again, test edits,
R-side work. The test files are read as text, not `include`d, so editing one
does not invalidate the cache; what it then calls that the cache lacks is
compiled at the time, as before.
"""
module ContinuousTimeSEMTestCache

using PrecompileTools
using Test, ContinuousTimeSEM
using LinearAlgebra, Random, Printf, ForwardDiff, DataFrames, ChainRulesCore,
    ComponentArrays

const TESTS = normpath(joinpath(@__DIR__, "..", ".."))

"""The test files whose compiling the cache holds: every one runtests.jl runs,
or those named in `CTSEM_TEST_CACHE_FILES` (comma separated) when the image is
built."""
function _cache_files()
    given = get(ENV, "CTSEM_TEST_CACHE_FILES", "")
    isempty(given) || return String.(strip.(split(given, ",")))
    return sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"),
        readdir(TESTS)))
end

# As text, so a test file is not a dependency of the image; see above.
_run_text(path) = Base.include_string(@__MODULE__, read(path, String), path)

@compile_workload begin
    _run_text(joinpath(TESTS, "table_helpers.jl"))
    for file in _cache_files()
        path = joinpath(TESTS, file)
        isfile(path) || continue
        # A failing test fails the run that reports it, not the cache build.
        try
            @testset "$file" begin
                _run_text(path)
            end
        catch err
            err isa Test.TestSetException || @warn "test cache: $file did not run" err
        end
    end
end

end
