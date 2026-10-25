# JuliaConnectoR's main.jl, starting the server from the CtsemRServer package
# (see its module docstring): the names main.jl leaves in `Main` -- `Pkg`,
# `Tables`, `RConnector`, and InteractiveUtils' exports -- bound the same way,
# so the R side, which asks `Main` for `RConnector` and checks `Tables`, finds
# what it would have. Arguments as main.jl's: port, port file, multiclient.
import Pkg
using CtsemRServer
const Tables = CtsemRServer.Tables
const RConnector = CtsemRServer.RConnector
using InteractiveUtils
RConnector.serve(parse(Int, ARGS[1]), portfile = ARGS[2],
   multiclient = (ARGS[3] == "t"))
