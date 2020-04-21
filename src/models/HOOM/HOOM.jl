using Distributed

@everywhere module HOOM

    include("include_packages.jl")


    include("../../share/MapInfo.jl")
    include("../../share/PolelikeCoordinate.jl")
    #include("../../share/PolelikeCoordinate.jl")
    include("../../share/constants.jl")
    #include("../../share/ocean_state_function.jl")


    include("OcnEnv.jl")
    include("SharedData.jl")
    include("JobDistributionInfo.jl")
    include("Model.jl")
   
    include("core_func.jl")


    # Sub models
    include("Dyn.jl")

    # Slave
    include("DynSlave.jl")
    #include("DynSlave.jl")

end
