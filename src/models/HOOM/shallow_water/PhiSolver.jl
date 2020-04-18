
mutable struct ΦSolver
    M :: DynamicAdvSpeedUpMatrix

#    ccT_send_T   :: AbstractArray{Float64, 2}
    cT_send_T   :: AbstractArray{Float64, 2}
#    T_send_ccT   :: AbstractArray{Float64, 2}
    T_send_cT   :: AbstractArray{Float64, 2}
 
    cT_Lap_cT    :: AbstractArray{Float64, 2}
    cT_MoLap_cT  :: AbstractArray{Float64, 2}
   
    T_Lap_T

    function ΦSolver(;
        gi             :: PolelikeCoordinate.GridInfo,
        mask2          :: AbstractArray{Float64, 2},
        α              :: Float64,  # ((Δt)^2 * H)^(-1)
    )
        M = DynamicAdvSpeedUpMatrix(;
            gi = gi,
            Nz = 1,
            mask2 = mask2,
        )

        Nx = gi.Nx
        Ny = gi.Ny

        # need a mask excluding bnd points
        #mask2_exclude_bnd = reshape(M.borderfilter_T * view(mask2, :), Nx, Ny)

        # Create coversion matrix and its inverse
        T_num = reshape(collect(1:length(mask2)), size(mask2)...)
        active_num             = T_num[ mask2 .==1 ]
        #active_num_exclude_bnd = T_num[ mask2_exclude_bnd .==1 ]


        cT_send_T  = M.op.T_I_T[active_num, :]
        #ccT_send_T = M.op.T_I_T[active_num_exclude_bnd, :]
        
        dropzeros!(cT_send_T)
        #dropzeros!(ccT_send_T)

        T_send_cT  = sparse(cT_send_T')
        #T_send_ccT = sparse(ccT_send_T')
       
        # identity 
        #ccT_I_ccT = ccT_send_T * T_send_ccT
        cT_I_cT = cT_send_T * T_send_cT

        # Laplacian on T grids with gradient equals zero at boundaries
        T_Lap_T   = M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T

        # ccT_Lap_ccT will get incorrect Laplacian since some info is lost during compression
        cT_Lap_cT = cT_send_T * T_Lap_T * T_send_cT

        dropzeros!(cT_Lap_cT)

        # Modified Laplacian
        cT_MoLap_cT = cT_Lap_cT - α * cT_I_cT
        
        return new(
            M,
        #    ccT_send_T,
            cT_send_T,
        #    T_send_ccT,
            T_send_cT,
            cT_Lap_cT,
            cT_MoLap_cT,
            T_Lap_T
        )

    end
end
