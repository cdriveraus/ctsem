"""
    CtsemRServer

The R bridge's Julia server -- JuliaConnectoR's `RConnector`, from the files
of the JuliaConnectoR installed with R -- as a package, so that its code is
compiled once into a package image instead of in every session.
JuliaConnectoR's `main.jl` includes the server into `Main` at run time, where
nothing can be cached, and a fresh ctsem session spent about four seconds
compiling it and another two compiling Base for its types
(review/PRECOMPILE-connector-proposal-2026-10-25.md).

ctsem builds this package from a template (R/ctJuliaServer.R): the installed
JuliaConnectoR's `Julia/*.jl` under `src/connector/`, with
`RConnector.jl`'s `import Main.Tables` pointed at this module's `Tables`, and
starts sessions from it (`serve.jl`). A session it serves is the one
`main.jl` would have started, the same server code under another module path.
"""
module CtsemRServer

# What `main.jl` loads before the server: Pkg, Tables and InteractiveUtils.
import Pkg
import Tables
import InteractiveUtils

include(joinpath("connector", "RConnector.jl"))

# The methods a ctsem session compiles in the server, and in Base for the
# server's types, from a session's `--trace-compile`
# (ctsem's dev/rserver-precompile.R). Compiled into this package's image as it
# is built; each read where the server's own module reads it. A statement that
# no longer applies -- a JuliaConnectoR or Julia whose methods differ -- is
# skipped, and its method compiles at run time as before.
const _PRECOMPILE = joinpath(@__DIR__, "precompile_statements.jl")
include_dependency(_PRECOMPILE)
if ccall(:jl_generating_output, Cint, ()) == 1
    for line in eachline(_PRECOMPILE)
        startswith(line, "precompile(") || continue
        try
            Core.eval(RConnector, Meta.parse(line))
        catch
        end
    end
end

end
