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

@testset "ContinuousTimeSEM" begin
    for file in TEST_FILES
        # Wrap each file in its own testset so failures point to the subsystem
        # under test while preserving the file-level organization.
        @testset "$(basename(file))" begin
            include(file)
        end
    end
end
