# One process per test file, longest first.
#
# `runtests.jl` includes every file into one process, which is the right default
# -- it is what `Pkg.test()` runs and it keeps the shared model objects shared.
# It is also serial, and on this suite that means 22 minutes of which
# `test_laplace.jl` is twelve, almost all of it Julia compiling the nested dual
# types the exact outer gradient needs.
#
# Nothing in the suite depends on running in one process: each file builds the
# models it uses. So the files can run side by side, and the wall clock drops to
# roughly the longest file plus its compilation. Cost is that each process
# compiles what it touches, so this is a win for a whole-suite run and a loss
# for a single file -- use `Pkg.test()` for one file.
#
#     julia --project=. test/runtests_parallel.jl          # all files
#     julia --project=. test/runtests_parallel.jl laplace  # files matching
#
# `CTSEM_TEST_WORKERS` caps the process count; it defaults to the core count.

using Printf

const HERE = @__DIR__
# `Pkg.test()` builds a temporary environment carrying the test-only
# dependencies; a bare `--project=<package>` has no DataFrames. Point
# `CTSEM_TEST_PROJECT` at an environment that has both, or run this from one.
const PROJECT = get(ENV, "CTSEM_TEST_PROJECT", Base.active_project())

files = sort(filter(readdir(HERE)) do name
    startswith(name, "test_") && endswith(name, ".jl")
end)
if !isempty(ARGS)
    files = filter(f -> any(a -> occursin(a, f), ARGS), files)
end
isempty(files) && (println("no test files matched ", ARGS); exit(1))

# Longest first, so the tail of the run is short files rather than one long one.
# The ordering is a hint from previous runs; an unlisted file sorts as average.
const KNOWN_SECONDS = Dict(
    "test_laplace.jl" => 750, "test_adjoint_gradient_validation.jl" => 195,
    "test_hessian.jl" => 155, "test_quadrature.jl" => 120,
    "test_kalman_trace.jl" => 40, "test_workspace_and_ekf.jl" => 40,
    "test_discrete_time.jl" => 39, "test_threading.jl" => 25,
    "test_ctsem_backend.jl" => 17, "test_adjoint_primitives.jl" => 15,
    "test_subject_gradients.jl" => 11, "test_summary_matrices.jl" => 6)
sort!(files; by = f -> -get(KNOWN_SECONDS, f, 20))

workers = parse(Int, get(ENV, "CTSEM_TEST_WORKERS", string(Sys.CPU_THREADS)))
workers = clamp(workers, 1, length(files))
@printf("%d files over %d processes\n", length(files), workers)

"""
    shards(file, count)

Split one test file's top-level `@testset` blocks into `count` groups.

The wall clock of the whole suite is one file: `test_laplace.jl` runs for ten
minutes while twenty other processes sit idle, so spreading files across cores
stops helping long before the cores run out. Its testsets are independent -- each
takes a fresh Laplace wrapper over a shared model -- so they can be split too.

The scan is lexical rather than a parse: a top-level testset in these files opens
with `@testset` in column one and closes with `end` in column one. Everything
outside those blocks is setup (helper functions, includes) and is kept in *every*
shard, since a shard that dropped it would not run. If the scan finds nothing to
split, the file is returned whole and nothing is lost.

Returns a vector of `(first_line, last_line)` ranges to keep, one per shard.
"""
function shards(file::AbstractString, count::Int)
    lines = readlines(file)
    starts = findall(l -> startswith(l, "@testset"), lines)
    length(starts) >= 2 * count || return [Int[]]
    stops = Int[]
    for s in starts
        stop = findfirst(i -> lines[i] == "end", s:length(lines))
        stop === nothing && return [Int[]]
        push!(stops, s + stop - 1)
    end
    blocks = collect(zip(starts, stops))
    # Round-robin rather than contiguous: neighbouring testsets in these files
    # tend to share a model, so consecutive ones cost about the same and dealing
    # them out spreads the expensive ones evenly.
    return [[i for i in eachindex(blocks) if mod1(i, count) == c] for c in 1:count]
end

"""The source of `file` with only the listed testset blocks kept."""
function shard_source(file::AbstractString, keep::Vector{Int})
    lines = readlines(file)
    isempty(keep) && return join(lines, "
")
    starts = findall(l -> startswith(l, "@testset"), lines)
    stops = [s + findfirst(i -> lines[i] == "end", s:length(lines)) - 1 for s in starts]
    drop = falses(length(lines))
    for (b, (s, e)) in enumerate(zip(starts, stops))
        b in keep && continue
        drop[s:e] .= true
    end
    return join(lines[.!drop], "
")
end

# Each file runs the same preamble `runtests.jl` uses, then just that file.
const PREAMBLE = """
using Test, ContinuousTimeSEM
using LinearAlgebra, Random, Printf, ForwardDiff, DataFrames, ChainRulesCore,
    ComponentArrays
cd(raw"$HERE")
include(joinpath(raw"$HERE", "table_helpers.jl"))
"""

# Split only the file that dominates the wall clock, and only in two.
#
# A shard is a fresh process, so it recompiles the fixtures it touches -- about
# two minutes of the ten this suite takes. That fixed cost is paid per shard, so
# splitting `test_laplace.jl` four ways cost four compilations to save three
# quarters of one file's arithmetic and came out no faster. Two is where the
# arithmetic saved still exceeds the compilation added.
#
# The real lever is not here: it is that a new model shape specialises the whole
# filter, so every process that touches a different model pays again. Removing
# that would shorten this suite and the first fit of a session alike.
# Empty on purpose. Splitting a file across processes was tried at two and four
# ways and neither was faster: a shard is a fresh process, so it recompiles the
# fixtures it touches -- roughly two minutes -- and that fixed cost cancels the
# arithmetic saved. The machinery is kept because it is correct and costs
# nothing while this is empty, but the lever is elsewhere: a new model shape
# specialises the whole filter, so every process touching a different model pays
# again. Fix that and this suite and the first fit of a session both get shorter.
const SHARD_INTO = Dict{String,Int}()

units = Tuple{String,Vector{Int},String}[]
for file in files
    count = get(SHARD_INTO, file, 1)
    groups = count > 1 ? shards(joinpath(HERE, file), count) : [Int[]]
    for (k, keep) in enumerate(groups)
        label = length(groups) > 1 ? "$(file) [$k/$(length(groups))]" : file
        push!(units, (file, keep, label))
    end
end
@printf("%d units over %d processes
", length(units), workers)

results = Dict{String,Tuple{Bool,Float64,String}}()
lock = ReentrantLock()
queue = Channel{Tuple{String,Vector{Int},String}}(length(units))
foreach(u -> put!(queue, u), units)
close(queue)

@sync for _ in 1:workers
    Threads.@spawn for (file, keep, label) in queue
        # Written out and `include`d rather than inlined into the script: a test
        # file declares `const`s at top level, and pasting it inside a `begin`
        # block turns those into local declarations, which is a syntax error.
        # `include` keeps the file's contents at top level, exactly as
        # `runtests.jl` does.
        # The shard goes *beside* its siblings, not in a temp directory: test
        # files reach their fixtures with `include(joinpath(@__DIR__, ...))`,
        # and `@__DIR__` is wherever the file being run happens to live.
        path = joinpath(HERE, file)
        if !isempty(keep)
            path = joinpath(HERE, ".shard-$(getpid())-$(hash(label))-$(file)")
            write(path, shard_source(joinpath(HERE, file), keep))
        end
        script = PREAMBLE * """
        @testset "$(label)" begin
            include(raw"$(path)")
        end
        """
        output = IOBuffer()
        started = time()
        ok = try
            run(pipeline(`$(Base.julia_cmd()) --project=$PROJECT -e $script`;
                stdout=output, stderr=output))
            true
        catch
            false
        end
        elapsed = time() - started
        isempty(keep) || rm(path; force = true)
        Base.lock(lock) do
            results[file] = (ok, elapsed, String(take!(output)))
            @printf("%-42s %s %6.1f s\n", file, ok ? "ok  " : "FAIL", elapsed)
        end
    end
end

failed = [f for f in files if !results[f][1]]
println()
@printf("%d/%d files passed, slowest %.1f s\n", length(files) - length(failed),
    length(files), maximum(r[2] for r in values(results); init=0.0))
for f in failed
    println("\n", "="^70, "\n", f, "\n", "="^70)
    # Only the interesting part: the failures and the summary line.
    for line in split(results[f][3], '\n')
        (occursin("Test Failed", line) || occursin("Error During", line) ||
         occursin("did not pass", line) || occursin("Expression:", line) ||
         occursin("Evaluated:", line) || occursin("Got exception", line) ||
         occursin("ERROR:", line) || occursin("Error:", line) ||
         startswith(strip(line), "nested task error")) && println("  ", line)
    end
end
exit(isempty(failed) ? 0 : 1)
