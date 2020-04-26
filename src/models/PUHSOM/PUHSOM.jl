using Distributed

@everywhere include("../../share/PolelikeCoordinate.jl")
@everywhere include("../../share/MapInfo.jl")
@everywhere include("Dyn/Dyn.jl")
@everywhere include("Tmd/Tmd.jl")

@everywhere module PUHSOM

    using Formatting
    using SharedArrays
    using Distributed
    using NCDatasets
    using ..Dyn
    using ..Tmd
    using ..ModelMap
    using ..PolelikeCoordinate

    include("OcnEnv.jl")
    include("SharedData.jl")
    include("JobDistributionInfo.jl")
    include("Model.jl")
   
    include("core_func.jl")

    # Slave
    include("DynSlave.jl")
    include("TmdSlave.jl")

end
