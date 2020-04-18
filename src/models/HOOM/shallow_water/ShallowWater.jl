
include("../../../share/DisplacedPoleCoordinate.jl")
include("../../../share/constants.jl")
include("../../../share/ocean_state_function.jl")

module ShallowWater

    using LinearAlgebra    
    using ..DisplacedPoleCoordinate


    function allocate(datakind::Symbol, dtype::DataType, dims... ; func=Main.zeros)
        if datakind == :local
            return func(dtype, dims...)
        elseif datakind == :shared
            return SharedArray{dtype}(dims...)
        else
            ErrorException("Unknown kind: " * string(datakind)) |> throw
        end
    end

    @inline function mul2!(
        a :: AbstractArray{Float64},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64},
    )

        mul!(view(a, :), b, view(c, :))

    end

    include("AdvectionSpeedUpMatrix.jl")
    include("Env.jl")
    include("State.jl")
    include("TracerAdv.jl")
    include("DynamicAdv.jl")
    include("Model.jl")
    include("step_tcr_adv.jl")
    include("step_dyn_adv.jl")

    include("../rearrange.jl")
    include("var_list.jl")

    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end



    function stepModel!(
        model :: Model,
        Δt    :: Float64,
    )

        #setupFlow!(model.state)
        
        advectTracer!(model, Δt)
        #advectDynamic!(model, Δt)

    end


end
