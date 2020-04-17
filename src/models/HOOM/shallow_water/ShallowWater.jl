
include("../../../share/DisplacedPoleCoordinate.jl")

module ShallowWater

    using LinearAlgebra    
    using ..DisplacedPoleCoordinate

    include("AdvectionSpeedUpMatrix.jl")
    include("structures.jl")
    include("advection.jl")

    include("../rearrange.jl")
    include("var_list.jl")

    macro loop_hor(ocn, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(ocn)).valid_idx)[2]

            $(esc(idx1)) = $(esc(ocn)).valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(ocn)).valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end


    function allocate(datakind::Symbol, dtype::DataType, dims... ; func=Main.zeros)
        if datakind == :local
            return func(dtype, dims...)
        elseif datakind == :shared
            return SharedArray{dtype}(dims...)
        else
            ErrorException("Unknown kind: " * string(datakind)) |> throw
        end
    end

    function advectTracer!(
        model   :: Model,
        Δt      :: Float64,
    )
   
        state = model.state
        tcr_adv = model.tcr_adv
        env = model.env

 
        calDIV!(
            gi    = env.gi,
            Nx    = env.Nx,
            Ny    = env.Ny,
            Nz    = env.Nz_av_f,
            u_bnd = state.u_f,
            v_bnd = state.v_f,
            div   = tcr_adv.div,
            mask3 = env.mask3_f,
        )
 
        calVerVelBnd!(
            gi    = env.gi,
            Nx    = env.Nx,
            Ny    = env.Ny,
            Nz    = env.Nz_av_f,
            w_bnd = state.w_f,
            hs    = env.H_f,
            div   = tcr_adv.div,
            mask3 = env.mask3_f,
        )
       

        # Pseudo code
        # 1. calculate tracer flux
        # 2. calculate tracer flux divergence
        calDiffAdv_QUICK_SpeedUp!(model, Δt)
        for x=1:env.NX
            for i = 1:env.Nx, j=1:env.Ny
                for k = 1:env.Nz_av_f[i, j]
                    state.X[k, i, j, x] += Δt * tcr_adv.XFLUX_CONV[k, i, j, x]
                end
            end
        end
    
    end


    function stepModel!(
        model :: Model,
        Δt    :: Float64,
    )

        #setupFlow!(model.state)
        
        advectTracer!(model, Δt)
        #advectDynamic!(model.dyn_adv, model.state, Δt)

    end

    mutable struct ABIIIObj
        gi :: DisplacedPoleCoordinate.GridInfo
        
    
        function ABIIIObj(
            ShallowWater
        )
                        
        end
    end


    function ABIII!(
        o :: ABIIIObj
    )

    end
end
