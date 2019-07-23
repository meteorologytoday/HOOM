using Distributed

@everywhere module SOM

using Printf
using Formatting
using SharedArrays
using Distributed
using NCDatasets

macro hinclude(path)
    return :(include(joinpath(@__DIR__, $path)))
end
    
@hinclude("../share/constants.jl")
@hinclude("OceanColumnCollection.jl")
@hinclude("doHorizontalDiffusion.jl")
@hinclude("stepOceanColumnCollection.jl")

@hinclude("takeSnapshot.jl")

end

