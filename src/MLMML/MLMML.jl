module MLMML
using Printf
using Formatting
using SparseArrays
using NCDatasets

include("../share/constants.jl")

include("OceanColumnCollection.jl")

include("trivial_functions.jl")

include("calWeOrMLD.jl")
include("doConvectiveAdjustment.jl")
include("doDiffusion.jl")
include("getIntegratedBuoyancy.jl")
include("stepOceanColumnCollection.jl")


include("setOceanColumn.jl")
include("takeSnapshot.jl")

end
