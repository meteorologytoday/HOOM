using Distributed

@everywhere module ESOM

    macro hinclude(path)
        return :(include(joinpath(@__DIR__, $path)))
    end
 
    @hinclude("../share/DisplacedPoleCoordinate.jl")
    @hinclude("../share/MapInfo.jl")

    using .DisplacedPoleCoordinate
    using .ModelMap


    using Printf
    using Formatting
    using SharedArrays
    using Distributed
    using SparseArrays
    using NCDatasets

       
    @hinclude("../share/constants.jl")
    @hinclude("OceanColumnCollection.jl")
    @hinclude("trivial_functions.jl")
    @hinclude("stepOceanColumnCollection.jl")
    @hinclude("takeSnapshot.jl")
    @hinclude("calEkmanTransport.jl")
    @hinclude("doConvectiveAdjustment.jl")

end

