using Test
using ContinuousTimeSEM

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

failed_files = String[]
@testset "ContinuousTimeSEM" begin
    for file in TEST_FILES
        # Wrap each file in its own testset so failures point to the subsystem
        # under test while preserving the file-level organization. A nested
        # testset records into its parent rather than throwing, so what it
        # returns can be inspected here and the run still ends red.
        ts = @testset "$(basename(file))" begin
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
    if ANNOTATE && !isempty(failed_files)
        println("::error title=Engine suite::", length(failed_files),
            " test files failed: ", join(failed_files, " "))
    end
end
