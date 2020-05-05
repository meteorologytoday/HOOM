using Distributed

@everywhere include("../../share/PolelikeCoordinate.jl")
@everywhere include("../../share/GridFiles.jl")
@everywhere include("../../share/RecordTool.jl")
@everywhere include("Dyn/Dyn.jl")
@everywhere include("Tmd/Tmd.jl")

@everywhere module BLOHSOM

    using Formatting
    using SharedArrays
    using Distributed
    using NCDatasets
    using ..Dyn
    using ..Tmd
    using ..GridFiles
    using ..RecordTool
    using ..PolelikeCoordinate

    include("../../share/constants.jl")

    include("Log.jl")
    include("OcnEnv.jl")
    include("DataUnit.jl")
    include("SharedData.jl")
    include("JobDistributionInfo.jl")
    include("DataExchanger.jl")
    include("VerticalAverager.jl")
    

    # Model included last
    include("Model.jl")
    include("core_func.jl")

    # Slave
    include("DynSlave.jl")
    include("DynSlave_func.jl")
    include("TmdSlave.jl")
    include("TmdSlave_func.jl")

    include("varlist.jl")
end
