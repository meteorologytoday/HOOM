include("../../share/PolelikeCoordinate.jl")
include("../../share/MapInfo.jl")

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


    using Formatting
    using LinearAlgebra    
    using ..PolelikeCoordinate
    using ..ModelMap
    using Statistics: mean

    include("../../share/constants.jl")
    include("../../share/ocean_state_function.jl")

    include("AdvectionSpeedUpMatrix.jl")
    include("TmdEnv.jl")
    include("TmdState.jl")
    include("TmdCore.jl")
    include("TmdModel.jl")

    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end



    function stepModel!(
        model :: TmdModel,
    )

        #setupFlow!(model.state)
        
        advectTracer!(model, Δt)

    end


end
