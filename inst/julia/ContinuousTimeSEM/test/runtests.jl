using Test
using ContinuousTimeSEM
# The same packages runtests_parallel.jl's PREAMBLE loads into each worker, so
# both runners give a test file the same names. Without them a file that used
# one without importing it passed locally and errored here, or not, depending
# on whether a file sorting before it happened to import it into Main.
using LinearAlgebra, Random, Printf, ForwardDiff, DataFrames, ChainRulesCore,
    ComponentArrays

# The tests' compile cache (test/cache), when this environment has it: loaded
# first, so the oracles find the engine already compiled at their dual types.
# Building it, when the engine has changed, costs what compiling the suite
# would; with it current, the suite is mostly running.
if Base.identify_package("ContinuousTimeSEMTestCache") !== nothing
    @eval using ContinuousTimeSEMTestCache
end

# Shared, test-only conveniences for writing parameter tables as DataFrames.
# The package itself takes plain column vectors; see src/r_interface.jl.
include("table_helpers.jl")

# Keep the harness maintenance-free: every new test file only needs the
# `test_*.jl` prefix to become part of the package test suite.
const TEST_FILES = sort(filter(
    file -> startswith(basename(file), "test_") && endswith(file, ".jl"),
    readdir(@__DIR__; join=true),
))

# Under GitHub Actions, say which file and which assertion in an annotation
# as well as in the log. A failed run's log is served only to someone who can
# write to the repository; its annotations are served to anyone. Without this
# the entire published account of a red engine run is "Process completed with
# exit code 1", and finding the assertion behind it means reproducing the
# whole suite locally -- which it has cost, at about fifty minutes a time.
const ANNOTATE = get(ENV, "GITHUB_ACTIONS", "") == "true"

# A finished testset holds its children, so the failures are there to be read
# rather than only printed. Recursive because a test file nests its own.
_failures(ts::Test.DefaultTestSet) =
    reduce(vcat, map(_failures, ts.results); init = Any[])
_failures(r::Union{Test.Fail,Test.Error}) = Any[r]
_failures(::Any) = Any[]

# How long each testset took, and how much of that was compiling. Julia keeps a
# testset's wall time but not its compile time, and here compiling is most of
# the cost: test_laplace.jl ran in 3.5 s once compiled, after 1711 s of
# compiling, each ForwardDiff oracle compiling the filter and the Laplace term
# for a dual type of its own. A thin wrapper around the default testset reads
# Julia's compile counter as a testset starts and finishes; every testset in a
# file inherits it, and everything else -- results, printing, failures -- is
# the default testset's. Written to `CTSEM_TEST_TIMINGS` (a CSV path) when that
# is set, and to the job summary under GitHub Actions.
Base.cumulative_compile_timing(true)
_compile_seconds() = Base.cumulative_compile_time_ns()[1] / 1e9

const TIMINGS = NamedTuple{(:file, :depth, :testset, :seconds, :compile),
    Tuple{String,Int,String,Float64,Float64}}[]
const CURRENT_FILE = Ref("")

mutable struct _TimedTestSet <: Test.AbstractTestSet
    inner::Test.DefaultTestSet
    compile0::Float64
end
_TimedTestSet(desc::AbstractString; kwargs...) =
    _TimedTestSet(Test.DefaultTestSet(desc; kwargs...), _compile_seconds())
Test.record(ts::_TimedTestSet, r) = Test.record(ts.inner, r)
if isdefined(Test, :get_rng)
    Test.get_rng(ts::_TimedTestSet) = Test.get_rng(ts.inner)
    Test.set_rng!(ts::_TimedTestSet, rng::Random.AbstractRNG) =
        Test.set_rng!(ts.inner, rng)
end
function Test.finish(ts::_TimedTestSet)
    push!(TIMINGS, (file = CURRENT_FILE[], depth = Test.get_testset_depth(),
        testset = ts.inner.description, seconds = time() - ts.inner.time_start,
        compile = _compile_seconds() - ts.compile0))
    # Records the default testset into the parent, as any nested one does.
    Test.finish(ts.inner)
    return ts
end
_failures(ts::_TimedTestSet) = _failures(ts.inner)

_csv(s::AbstractString) = string('"', replace(s, '"' => "\"\""), '"')
# Capped: the compile counter and the wall clock are read a moment apart.
_share(part, whole) = min(100, round(Int, 100 * part / max(whole, eps())))

function _timing_report(rows)
    path = get(ENV, "CTSEM_TEST_TIMINGS", "")
    if !isempty(path)
        open(path, "w") do io
            println(io, "file,depth,testset,seconds,compile")
            for r in rows
                println(io, join((r.file, r.depth, _csv(r.testset),
                    round(r.seconds; digits = 2), round(r.compile; digits = 2)), ","))
            end
        end
    end
    files = sort(filter(r -> r.depth == 1, rows); by = r -> -r.seconds)
    total = sum(r -> r.seconds, files; init = 0.0)
    slow = first(sort(filter(r -> r.depth == 2, rows); by = r -> -r.seconds), 20)
    @printf("\nEngine tests: %.1f min, %d%% of it compiling. Slowest testsets:\n",
        total / 60, _share(sum(r -> r.compile, files; init = 0.0), total))
    for r in slow
        @printf("  %7.1f s  %3d%% compiling  %s: %s\n", r.seconds,
            _share(r.compile, r.seconds), r.file, r.testset)
    end
    summary = get(ENV, "GITHUB_STEP_SUMMARY", "")
    isempty(summary) && return nothing
    open(summary, "a") do io
        println(io, "### Engine tests: time, and the share spent compiling\n")
        println(io, "| file | seconds | compiling |\n|---|---:|---:|")
        for r in files
            @printf(io, "| %s | %.1f | %d%% |\n", r.file, r.seconds,
                _share(r.compile, r.seconds))
        end
        println(io, "\n| slowest testsets | seconds | compiling |\n|---|---:|---:|")
        for r in slow
            @printf(io, "| %s: %s | %.1f | %d%% |\n", r.file,
                replace(r.testset, '|' => '/'), r.seconds, _share(r.compile, r.seconds))
        end
    end
    return nothing
end

failed_files = String[]
@testset "ContinuousTimeSEM" begin
    for file in TEST_FILES
        # Wrap each file in its own testset so failures point to the subsystem
        # under test while preserving the file-level organization. A nested
        # testset records into its parent rather than throwing, so what it
        # returns can be inspected here and the run still ends red.
        CURRENT_FILE[] = basename(file)
        ts = @testset _TimedTestSet "$(basename(file))" begin
            include(file)
        end
        bad = _failures(ts)
        isempty(bad) && continue
        push!(failed_files, basename(file))
        ANNOTATE || continue
        # GitHub serves ten annotations per run and the summary below has to
        # be one of them, so the detail stops short of that.
        for f in Iterators.take(bad, 3)
            msg = replace(sprint(show, f), r"[\r\n]+" => " ")
            println("::error title=", basename(file), "::", first(msg, 600))
        end
    end
    # Inside the block, because the outermost testset throws as it finishes
    # and nothing after `end` would run.
    _timing_report(TIMINGS)
    if ANNOTATE && !isempty(failed_files)
        println("::error title=Engine suite::", length(failed_files),
            " test files failed: ", join(failed_files, " "))
    end
end
