
include("../../../share/PolelikeCoordinate.jl")

module Dyn

    using Formatting
    using LinearAlgebra    
    using ..PolelikeCoordinate
    using Statistics: mean

    include("../../../share/constants.jl")
    include("../../../share/ocean_state_function.jl")




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
 
    const α_AB = 1.0 /  2.0 
    const β_AB = 5.0 / 12.0  
    @inline function ABIII(a, aa, aaa)
        return ( 1 + α_AB + β_AB ) * a - (α_AB + 2*β_AB) * aa + β_AB * aaa
    end
   

    include("../MatrixOperators.jl")
    include("../VerticalAverager.jl")

    include("AdvectionSpeedUpMatrix_dyn.jl")
    include("PhiSolver.jl")
    include("DiffusionSolver.jl")
    include("Workspace.jl")


    include("DynEnv.jl")
    include("DynState.jl")
    include("DynCore.jl")
    include("DynModel.jl")
    include("step_model.jl")


    include("varlist.jl")

    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end



    function stepModel!(
        model :: DynModel,
    )
        reset!(model.core.wksp)
        advectDynamic!(model)
    end


end
