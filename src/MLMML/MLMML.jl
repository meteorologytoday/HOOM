module MLMML
using Printf
using Formatting
using SparseArrays

include("../share/constants.jl")
include("extra_constants.jl")

include("Workspace.jl")
include("OceanColumnCollection.jl")

include("trivial_functions.jl")

include("calWeOrMLD.jl")
include("doConvectiveAdjustment.jl")
include("doDiffusion.jl")
include("getIntegratedBuoyancy.jl")
include("stepOceanColumnCollection.jl")


include("setOceanColumn.jl")

end
