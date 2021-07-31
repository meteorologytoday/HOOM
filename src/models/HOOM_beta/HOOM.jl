using Distributed

macro hinclude(path)
    return :(include(normpath(joinpath(@__DIR__, $path))))
end



@everywhere module HOOM

    using Dates
    using Printf
    using Formatting
    using SharedArrays
    using Distributed
    using SparseArrays
    using NCDatasets

    macro hinclude(path)
        return :(include(normpath(joinpath(@__DIR__, $path))))
    end
 
    @hinclude("../../share/constants.jl")
    @hinclude("../../share/ocean_state_function.jl")

    # classes
    @hinclude("../../share/GridFile.jl")
    @hinclude("../../share/MapInfo.jl")
    @hinclude("../../share/PolelikeCoordinate.jl")
    @hinclude("../../share/BasicMatrixOperators.jl")
    @hinclude("../../share/AdvancedMatrixOperators.jl")

    @hinclude("Env.jl")
#    @hinclude("State.jl")
#    @hinclude("Core.jl")

end



