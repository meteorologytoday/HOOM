module Tmd

    @inline function cyc(i::Int64, N::Int64)
        return mod(i-1, N) + 1
    end

    @inline function mul2!(
        a :: AbstractArray{Float64, 2},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64, 2},
    )
        mul!(view(a, :), b, view(c, :))
    end
 
    @inline function mul3!(
        a :: AbstractArray{Float64, 3},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64, 3},
    )
        for k=1:size(a)[3]
            mul!(
                view(view(a, :, :, k), :),
                b,
                view(view(c, :, :, k), :),
            )
        end
    end



    macro fast_extract(model)
        return esc(:( 
            co = $(model).core;
            st = $(model).state;
            dg = $(model).diag;
            fr = $(model).forcing;
            ev = $(model).env;
        ))
    end




    macro loop_hor(model, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(model)).core.valid_idx)[2]

            $(esc(idx1)) = $(esc(model)).core.valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(model)).core.valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end



    using Formatting
    using LinearAlgebra    
    using ..PolelikeCoordinate
    using ..ModelMap
    using Statistics: mean

    include("../../../share/constants.jl")
    include("../../../share/ocean_state_function.jl")

    include("Workspace.jl")
    include("AdvectionSpeedUpMatrix.jl")
    include("TmdEnv.jl")
    include("TmdState.jl")
    include("TmdDiag.jl")
    include("TmdForcing.jl")
    include("TmdCore.jl")
    include("TmdModel.jl")

    # functions
    include("latent_heat_release_of_freezing.jl")
    include("columnwise_budget.jl")
    include("trivial_functions.jl")

    include("mld_calculation.jl")
    include("convective_adjustment.jl")
    include("diffusion.jl")
    include("mixUnmix.jl")
    include("calFLDOPartition.jl")
    include("columnwise_integration.jl")
    #include("deep_ocn_correction.jl")
    include("shortwave_radiation.jl")
    #include("flx_correction.jl")
    
 
    include("varlist.jl")
    include("step_model.jl")
    include("step_tmd_mixed_layer.jl")

    function stepModel!(
        model :: TmdModel,
    )
        reset!(model.core.wksp)
        doMixedLayerDynamics!(model)
        advectTracer!(model)
        
    end


end
