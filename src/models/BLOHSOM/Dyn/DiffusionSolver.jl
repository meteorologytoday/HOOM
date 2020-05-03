
using LinearAlgebra

mutable struct DiffusionSolver

    UU_length    :: Int64
    VV_length    :: Int64

    M :: DynamicAdvSpeedUpMatrix
   
    D  :: Float64
    Δt :: Float64
    α  :: Float64

    UU_send_U    :: AbstractArray{Float64, 2}
    U_send_UU    :: AbstractArray{Float64, 2}

    VV_send_V    :: AbstractArray{Float64, 2}
    V_send_VV    :: AbstractArray{Float64, 2}

    UU_MoLap_UU    :: AbstractArray{Float64, 2}
    VV_MoLap_VV    :: AbstractArray{Float64, 2}

    tool_mtx

    wksp_UU        :: Array
    wksp_VV        :: Array
 
    function DiffusionSolver(;
        gi             :: PolelikeCoordinate.GridInfo,
        mask2          :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
        D              :: Float64,
        Δt             :: Float64,
        M              :: Union{DynamicAdvSpeedUpMatrix, Nothing} = nothing
    )

        α = 1.0 / (D * Δt)

        if M == nothing
            println("Mask not provided. Use mask provided by GridInfo.")
            M = DynamicAdvSpeedUpMatrix(;
                gi = gi,
                Nz = 1,
                mask2 = mask2,
            )
        end

        Nx = gi.Nx
        Ny = gi.Ny
        Nyp1 = Ny+1

        # need a mask excluding bnd points
        # Create coversion matrix and its inverse
        println(size(M.filter_U))
        println(Nx, ",", Ny)
        U_num        = reshape(M.filter_U * collect(1:(Nx*Ny)), Nx, Ny)
        active_U_num = U_num[ U_num .!= 0.0 ]
        UU_send_U    = M.op.U_I_U[active_U_num, :]
        U_send_UU    = UU_send_U' |> sparse
        UU_I_UU      = UU_send_U * U_send_UU
        UU_Lap_UU    = UU_send_U * M.U_Lap_U * U_send_UU
        UU_MoLap_UU  = UU_Lap_UU - α * UU_I_UU
        MoLap_UU     = lu(UU_MoLap_UU)


        V_num        = reshape(M.filter_V * collect(1:(Nx*Nyp1)), Nx, Nyp1)
        active_V_num = V_num[ V_num .!= 0.0 ]
        VV_send_V    = M.op.V_I_V[active_V_num, :]
        V_send_VV    = VV_send_V' |> sparse
        VV_I_VV      = VV_send_V * V_send_VV
        VV_Lap_VV    = VV_send_V * M.V_Lap_V * V_send_VV
        VV_MoLap_VV  = VV_Lap_VV - α * VV_I_VV
        MoLap_VV     = lu(VV_MoLap_VV)

        tool_mtx = (
            MoLap_VV = MoLap_VV,
            MoLap_UU = MoLap_UU,
        )
        
        UU_length = size(UU_MoLap_UU)[1]
        VV_length = size(VV_MoLap_VV)[1]

        wksp_UU = [ zeros(Float64, UU_length), zeros(Float64, UU_length) ]
        wksp_VV = [ zeros(Float64, VV_length), zeros(Float64, VV_length) ]

        return new(
            UU_length,
            VV_length,
            M,
            D,
            Δt,
            α,
            UU_send_U,
            U_send_UU,
            VV_send_V,
            V_send_VV,
            UU_MoLap_UU,
            VV_MoLap_VV,
            tool_mtx,
            wksp_UU,
            wksp_VV,
        )

    end
end


function solveDiffusion!(
    ds     :: DiffusionSolver,
    grid   :: Symbol,
    input  :: AbstractArray{Float64},
    output :: AbstractArray{Float64},
)

    if grid == :U
        rhs = ds.wksp_UU[1]
        lhs = ds.wksp_UU[2]
        tool = ds.tool_mtx.MoLap_UU
        send    = ds.UU_send_U
        invsend = ds.U_send_UU
    elseif grid == :V
        rhs = ds.wksp_VV[1]
        lhs = ds.wksp_VV[2]
        tool = ds.tool_mtx.MoLap_VV
        send    = ds.VV_send_V
        invsend = ds.V_send_VV
    else
        throw(ErrorException("Unrecognized grid: " * string(grid)))
    end

    mul!(rhs, send, view(input, :))
    rhs .*= -ds.α
    ldiv!(lhs, tool, rhs)
    mul!(view(output, :), invsend, lhs)

end
