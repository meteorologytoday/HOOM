
mutable struct PhiSolver

    TT_length    :: Int64

    M :: DynamicAdvSpeedUpMatrix
    α :: Float64

    TT_send_T    :: AbstractArray{Float64, 2}
    T_send_TT    :: AbstractArray{Float64, 2}
    TT_Lap_TT    :: AbstractArray{Float64, 2}
    TT_MoLap_TT  :: AbstractArray{Float64, 2}
   
    T_Lap_T
    tool_mtx

    function PhiSolver(;
        gi             :: PolelikeCoordinate.GridInfo,
        mask2          :: AbstractArray{Float64, 2},
        α              :: Float64,  # ((Δt)^2 * H)^(-1)
        M              :: Union{DynamicAdvSpeedUpMatrix, Nothing} = nothing
    )
        if M == nothing
            M = DynamicAdvSpeedUpMatrix(;
                gi = gi,
                Nz = 1,
                mask2 = mask2,
            )
        end

        Nx = gi.Nx
        Ny = gi.Ny

        # need a mask excluding bnd points
        #mask2_exclude_bnd = reshape(M.borderfilter_T * view(mask2, :), Nx, Ny)

        # Create coversion matrix and its inverse
        T_num = reshape(collect(1:length(mask2)), size(mask2)...)
        active_num             = T_num[ mask2 .==1 ]
        #active_num_exclude_bnd = T_num[ mask2_exclude_bnd .==1 ]


        TT_send_T  = M.op.T_I_T[active_num, :]
        #cTT_send_T = M.op.T_I_T[active_num_exclude_bnd, :]
        
        dropzeros!(TT_send_T)
        #dropzeros!(cTT_send_T)

        T_send_TT  = sparse(TT_send_T')
        #T_send_cTT = sparse(cTT_send_T')
       
        # identity 
        #cTT_I_cTT = cTT_send_T * T_send_cTT
        TT_I_TT = TT_send_T * T_send_TT

        # Laplacian on T grids with gradient equals zero at boundaries (U, V grid boundaries)
        T_Lap_T   = M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T


        # cTT_Lap_cTT will get incorrect Laplacian since some info is lost during compression
        TT_Lap_TT = TT_send_T * T_Lap_T * T_send_TT

        dropzeros!(TT_Lap_TT)

        #println(Array(TT_Lap_TT))

        # Modified Laplacian
        TT_MoLap_TT = TT_Lap_TT - α * TT_I_TT


#=
        # Testing       
        rr, cc = size(TT_Lap_TT)

#        println(TT_Lap_TT.nzval)

        r_cnt = 0
        for i=1:rr
            if all(TT_Lap_TT[i, :] .<= 1e-13)
                r_cnt += 1
            end
        end

        c_cnt = 0
        for j=1:cc
            if all(TT_Lap_TT[:, j] .<= 1e-13)
                c_cnt += 1
            end
        end

        println("Empty rows: ", r_cnt)
        println("Empty cols: ", c_cnt)
=#

        MoLap = lu(TT_MoLap_TT)
        tool_mtx = (
        #    Lap   = lu(TT_Lap_TT),
            MoLap = lu(TT_MoLap_TT),
        )
        
        return new(
            size(TT_MoLap_TT)[1],
            M,
            α,
        #    cTT_send_T,
            TT_send_T,
        #    T_send_cTT,
            T_send_TT,
            TT_Lap_TT,
            TT_MoLap_TT,
            T_Lap_T,
            tool_mtx,
        )

    end
end
