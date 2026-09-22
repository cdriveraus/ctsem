# Find the docstring pairing that parses cleanly and fails at *load*.
#
#     julia dev/check-julia-docstrings.jl
#
# Two docstrings with no expression between them -- which is what you get by
# deleting the thing one of them documented, or by inserting a documented
# function just above an existing docstring -- parse without complaint and then
# fail with
#
#     ERROR: LoadError: cannot document the following expression:
#
# `Meta.parseall` does not catch it, so the only other way to find it is to
# precompile the engine, which is about two hundred seconds. This is instant.
# It cost three load failures in one session, which is why it exists.
#
# A block that is the file's *first* expression is exempt: several engine files
# open with a section overview that documents nothing, and those load. Anywhere
# else, a docstring must be followed by something to document.
#
# Exit status is 1 when anything is flagged, so it can gate a commit.

const SRC = joinpath(@__DIR__, "..", "inst", "julia", "ContinuousTimeSEM", "src")

"""Line numbers on which a docstring block ends, and whether it opened line 1."""
function docstring_ends(lines)
    ends = Tuple{Int,Bool}[]
    open_at = 0
    for (i, l) in pairs(lines)
        # Occurrences of the triple quote on this line. A regex, not byte
        # slicing: these files contain non-ASCII (the engine writes theta and
        # Frechet with their accents), and indexing into the middle of a
        # multi-byte character throws.
        n = length(collect(eachmatch(r"\"\"\"", l)))
        if n >= 2 && open_at == 0
            push!(ends, (i, i == 1))              # one-line docstring
        elseif n == 1
            if open_at == 0
                open_at = i
            else
                push!(ends, (i, open_at == 1))
                open_at = 0
            end
        end
    end
    return ends
end

problems = 0
for file in sort(filter(f -> endswith(f, ".jl"), readdir(SRC)))
    lines = split(read(joinpath(SRC, file), String), r"\r?\n")
    for (e, first_in_file) in docstring_ends(lines)
        first_in_file && continue
        j = e + 1
        while j <= length(lines) && isempty(strip(lines[j]))
            j += 1
        end
        if j <= length(lines) && startswith(lstrip(lines[j]), "\"\"\"")
            println("$file:$e is a docstring documenting the docstring at $j")
            global problems += 1
        end
    end
end

if problems == 0
    println("no back-to-back docstrings")
else
    println("$problems problem(s): each will fail at load, not at parse")
    exit(1)
end
