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
 
    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end

    @inline function mul_autoflat!(
        a :: AbstractArray{Float64},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64},
    )
        mul!(view(a, :), b, view(c, :))
    end
 
    @hinclude("../../share/constants.jl")
    @hinclude("../../share/ocean_state_function.jl")

    # classes
    @hinclude("../../share/DisplacedPoleCoordinate.jl")
    @hinclude("../../share/MapInfo.jl")
    @hinclude("Workspace.jl")
    @hinclude("MatrixOperators.jl")
    @hinclude("AdvectionSpeedUpMatrix.jl")
    @hinclude("InputFields.jl")
    @hinclude("AccumulativeVariables.jl")
    @hinclude("Ocean.jl")

    # functions
    @hinclude("latent_heat_release_of_freezing.jl")
    @hinclude("columnwise_budget.jl")
    @hinclude("trivial_functions.jl")

    @hinclude("mld_calculation.jl")
    @hinclude("convective_adjustment.jl")
    @hinclude("diffusion.jl")
    @hinclude("mixUnmix.jl")
    @hinclude("calTsSsMixed.jl")
    @hinclude("calFLDOPartition.jl")
    @hinclude("columnwise_integration.jl")
    @hinclude("deep_ocn_correction.jl")
    @hinclude("shortwave_radiation.jl")
    @hinclude("flx_correction.jl")
    @hinclude("nudge_seaice.jl")
    
    @hinclude("step_ocean_prepare.jl")
    @hinclude("step_ocean_hz.jl")
    @hinclude("step_ocean_vt.jl")
    @hinclude("advection.jl")


    @hinclude("setOceanColumn.jl")
    @hinclude("takeSnapshot.jl")
    @hinclude("rearrange.jl")

    @hinclude("accumulate.jl")
    @hinclude("driver.jl")

    @hinclude("var_desc.jl")
    @hinclude("var_list.jl")




end



