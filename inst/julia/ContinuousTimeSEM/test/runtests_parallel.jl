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

# Each file runs the same preamble `runtests.jl` uses, then just that file.
const PREAMBLE = """
using Test, ContinuousTimeSEM
using LinearAlgebra, Random, Printf, ForwardDiff, DataFrames, ChainRulesCore,
    ComponentArrays
cd(raw"$HERE")
include(joinpath(raw"$HERE", "table_helpers.jl"))
"""

results = Dict{String,Tuple{Bool,Float64,String}}()
lock = ReentrantLock()
queue = Channel{String}(length(files))
foreach(f -> put!(queue, f), files)
close(queue)

@sync for _ in 1:workers
    Threads.@spawn for file in queue
        script = PREAMBLE * """
        @testset "$(file)" begin
            include(joinpath(raw"$HERE", "$(file)"))
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
