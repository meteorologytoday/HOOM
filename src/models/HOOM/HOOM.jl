using Distributed

macro hinclude(path)
    return :(include(normpath(joinpath(@__DIR__, $path))))
end



@everywhere module HOOM

    using Printf
    using Formatting
    using SharedArrays
    using Distributed
    using SparseArrays
    using NCDatasets

    macro hinclude(path)
        return :(include(normpath(joinpath(@__DIR__, $path))))
    end
 
    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
       #= 
        return :( for $(esc(idx1))=1:$(esc(ocn)).Nx, $(esc(idx2))=1:$(esc(ocn)).Ny

            $(esc(ocn)).mask[$(esc(idx1)), $(esc(idx2))] == 0.0 && continue
            $(esc(stmts))

        end )
        =#
    end


    @hinclude("../../share/DisplacedPoleCoordinate.jl")
    @hinclude("../../share/MapInfo.jl")
       
    @hinclude("../../share/constants.jl")
    @hinclude("InputFields.jl")
    @hinclude("AccumulativeVariables.jl")
    @hinclude("SpeedUpMtx3D.jl")
    @hinclude("Ocean.jl")

    @hinclude("latent_heat_release_of_freezing.jl")
    @hinclude("columnwise_budget.jl")
    @hinclude("trivial_functions.jl")

    @hinclude("calNewMLD.jl")
    @hinclude("doConvectiveAdjustment.jl")
    @hinclude("doDiffusion.jl")
    @hinclude("mixUnmix.jl")
    @hinclude("calFLDOPartition.jl")
    @hinclude("getIntegratedBuoyancy.jl")
    @hinclude("deep_ocn_correction.jl")
    @hinclude("shortwave_radiation.jl")
    @hinclude("flx_correction.jl")
    
    @hinclude("stepOcean_prepare.jl")
    @hinclude("stepOcean_hz.jl")
    @hinclude("stepOcean_vt.jl")
    @hinclude("advection.jl")


    @hinclude("setOceanColumn.jl")
    @hinclude("takeSnapshot.jl")
    @hinclude("rearrange.jl")

    @hinclude("accumulate.jl")
    @hinclude("driver.jl")

    @hinclude("varlist.jl")



end



