
using SparseArrays

@inline function speye(dtype, n)
    return spdiagm(0=>ones(dtype, n))
end

mutable struct DynamicAdvSpeedUpMatrix
    
    op       :: MatrixOperators

    # Nomenclature:
    #
    # [new-grid][function][old-grid]
    #
    # U  : U grid
    # V  : V grid
    # T  : T grid
    # S  : S grid
    # RU : Reduced U grid
    # 
    # U_W_T : sending variable westward from T grid to U grid

    U_∂x_T   :: AbstractArray{Float64, 2}   # ∇b 
    V_∂y_T   :: AbstractArray{Float64, 2}   
 
    U_∂x_U   :: AbstractArray{Float64, 2}   # ∇u
    U_∂y_U   :: AbstractArray{Float64, 2}

    V_∂x_V   :: AbstractArray{Float64, 2}   # ∇v
    V_∂y_V   :: AbstractArray{Float64, 2}

    T_DIVx_U :: AbstractArray{Float64, 2}
    T_DIVy_V :: AbstractArray{Float64, 2}
    
    T_Lap_T  :: AbstractArray{Float64, 2}
    U_Lap_U  :: AbstractArray{Float64, 2}
    V_Lap_V  :: AbstractArray{Float64, 2}
 
    U_interp_V :: AbstractArray{Float64, 2}  # interpolation of V grid onto U grid
    V_interp_U :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid

    T_interp_U :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    T_interp_V :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    U_interp_T :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    V_interp_T :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid

    U_f_V      :: AbstractArray{Float64, 2}   # used to get fv on U grid
    V_f_U      :: AbstractArray{Float64, 2}   # used to get fu on V grid
    
    filter_T       :: AbstractArray{Float64, 2}
    filter_U       :: AbstractArray{Float64, 2}
    filter_V       :: AbstractArray{Float64, 2}
    filter_F       :: AbstractArray{Float64, 2}
    borderfilter_T :: AbstractArray{Float64, 2}

    function DynamicAdvSpeedUpMatrix(;
        gi             :: PolelikeCoordinate.GridInfo,
        Nz             :: Int64,
        mask2          :: AbstractArray{Float64, 2},
    )

       #println("Begin")
        Nx = gi.Nx
        Ny = gi.Ny

        #println("Build MatrixOperators")
        @time op = MatrixOperators(Nx=Nx, Ny=Ny, Nz=Nz)

        # making filter
        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt23 = (x,) -> repeat(reshape(x, 1, size(x)...), outer=(Nz, 1, 1))
        cvt23_diagm = (x,) -> spdiagm( 0 => view(cvt23(x), :) )

        #println("Build filters")
        mask3_flat = view(mask2 |> cvt23, :)
        
        onV_if_unblocked_north_onT = op.V_S_T * mask3_flat
        onV_if_unblocked_south_onT = op.V_N_T * mask3_flat
        onU_if_unblocked_east_onT  = op.U_W_T * mask3_flat
        onU_if_unblocked_west_onT  = op.U_E_T * mask3_flat

        V_mask = onV_if_unblocked_north_onT .* onV_if_unblocked_south_onT
        U_mask = onU_if_unblocked_east_onT  .* onU_if_unblocked_west_onT


        onF_if_unblocked_north_onU = op.F_S_U * U_mask
        onF_if_unblocked_south_onU = op.F_N_U * U_mask
        onF_if_unblocked_east_onV  = op.F_W_V * V_mask
        onF_if_unblocked_west_onV  = op.F_E_V * V_mask

        F_mask = onF_if_unblocked_north_onU .* onF_if_unblocked_south_onU .* onF_if_unblocked_east_onV .* onF_if_unblocked_west_onV

        #=
        if_mask_northeast = op.F_SW_T * mask3_flat
        if_mask_northwest = op.F_SE_T * mask3_flat
        if_mask_southeast = op.F_SW_T * mask3_flat
        if_mask_southwest = op.F_SW_T * mask3_flat
        =#

        #println(size(op.V_S_T))
        #println(size(mask3_flat))

        filter_T = spdiagm(0 => mask3_flat)
        filter_V = spdiagm(0 => V_mask)
        filter_U = spdiagm(0 => U_mask)
        filter_F = spdiagm(0 => F_mask)

        borderfilter_T = spdiagm( 0 => ( 
               (op.T_N_T * mask3_flat)
            .* (op.T_S_T * mask3_flat)
            .* (op.T_E_T * mask3_flat)
            .* (op.T_W_T * mask3_flat)
        ))

        #println("if_mask_north_V", size(if_mask_north_V))
        #println("filter_V", size(filter_V))

        # ===== [ BEGIN face area and lengths on U V ] =====
        
        Δx_U = gi.dx_w
        Δy_U = gi.DY
 
        Δx_V = gi.DX
        Δy_V = hcat(gi.dy_s, gi.dy_n[:, end])

        V_Δx_V    = (Δx_V         |> cvt23_diagm)
        U_Δy_U    = (Δy_U         |> cvt23_diagm)

        U_invΔx_U = (Δx_U.^(-1)   |> cvt23_diagm)
        U_invΔy_U = (Δy_U.^(-1)   |> cvt23_diagm)
        
        V_invΔx_V = (Δx_V.^(-1)   |> cvt23_diagm)
        V_invΔy_V = (Δy_V.^(-1)   |> cvt23_diagm)

        Δσ_U = Δx_U .* Δy_U
        Δσ_V = Δx_V .* Δy_V

        U_invΔσ_U = Δσ_U.^(-1) |> cvt23_diagm
        V_invΔσ_V = Δσ_V.^(-1) |> cvt23_diagm

        # ===== [ END face area and lengths on U V ] =====
        
        # ===== [ BEGIN face area and lengths on T ] =====

        Δx_T = gi.dx_c
        Δy_T = gi.dy_c
        Δσ_T = Δx_T .* Δy_T

        T_Δx_T = (Δx_T |> cvt23_diagm)
        T_Δy_T = (Δy_T |> cvt23_diagm)
        
        T_invΔx_T = (Δx_T.^(-1) |> cvt23_diagm)
        T_invΔy_T = (Δy_T.^(-1) |> cvt23_diagm)
        T_invΔσ_T = (Δσ_T.^(-1) |> cvt23_diagm)

        # ===== [ END face area and lengths on T ] =====

        # ===== [ BEGIN length on F_grid ] =====
        # F grid, tricky. Need to worry about north and south boundary since F grid size is ill-defined.
        Δx_F = (circshift(gi.DX, (1, 0)) + gi.DX) / 2.0 

        Δy_F = zeros(Float64, size(Δx_F)...)
        Δy_F[:, 2:Ny] = (gi.DY[:, 1:end-1] + gi.DY[:, 2:end]) / 2.0
        Δy_F[:, 1   ] = gi.DY[:,   1] / 2.0
        Δy_F[:, end ] = gi.DY[:, end] / 2.0
       
        invΔx_F = Δx_F.^(-1)
        invΔy_F = Δy_F.^(-1)

        if any(isnan.(invΔx_F)) || any(isnan.(invΔy_F))
            throw(ErrorException("Contains NaN. Maybe Δx_F or Δy_F contains zero length."))
        end

        F_Δx_F = Δx_F |> cvt23_diagm
        F_Δy_F = Δy_F |> cvt23_diagm

        F_invΔx_F = invΔx_F |> cvt23_diagm
        F_invΔy_F = invΔy_F |> cvt23_diagm

        # ===== [ END length on F_grid ] =====

        #println("Making various operators")

        # MAGIC!!
        U_∂x_T = filter_U * U_invΔx_U * (op.U_W_T - op.U_E_T) ; dropzeros!(U_∂x_T);
        V_∂y_T = filter_V * V_invΔy_V * (op.V_S_T - op.V_N_T) ; dropzeros!(V_∂y_T);

        # notice here no filter is needed. No extrapolated information.
        T_∂x_U = T_invΔx_T * (op.T_W_U - op.T_E_U) ; dropzeros!(T_∂x_U);
        T_∂y_V = T_invΔy_T * (op.T_S_V - op.T_N_V) ; dropzeros!(T_∂y_V);

        U_∂x_U = filter_U * U_invΔx_U * (op.U_W_U - op.U_E_U) / 2.0 ; dropzeros!(U_∂x_U);
        U_∂y_U = filter_U * U_invΔy_U * (op.U_S_U - op.U_N_U) / 2.0 ; dropzeros!(U_∂y_U);

        V_∂x_V = filter_V * V_invΔx_V * (op.V_W_V - op.V_E_V) / 2.0 ; dropzeros!(V_∂x_V);
        V_∂y_V = filter_V * V_invΔy_V * (op.V_S_V - op.V_N_V) / 2.0 ; dropzeros!(V_∂y_V);

        F_∂y_U = filter_F * F_invΔy_F * (op.F_S_U - op.F_N_U) ; dropzeros!(F_∂y_U); 
        F_∂x_V = filter_F * F_invΔx_F * (op.F_W_V - op.F_E_V) ; dropzeros!(F_∂x_V); 

        T_DIVx_U = filter_T * T_invΔσ_T * ( op.T_W_U - op.T_E_U  ) * U_Δy_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V = filter_T * T_invΔσ_T * ( op.T_S_V - op.T_N_V  ) * V_Δx_V  ; dropzeros!(T_DIVy_V);
        
        U_DIVx_T = filter_U * U_invΔσ_U * ( op.U_W_T - op.U_E_T  ) * T_Δy_T  ; dropzeros!(U_DIVx_T);
        V_DIVy_T = filter_V * V_invΔσ_V * ( op.V_S_T - op.V_N_T  ) * T_Δx_T  ; dropzeros!(V_DIVy_T);

        V_DIVx_F = filter_V * V_invΔσ_V * ( op.V_W_F - op.V_E_F  ) * F_Δy_F  ; dropzeros!(V_DIVx_F);
        U_DIVy_F = filter_U * U_invΔσ_U * ( op.U_S_F - op.U_N_F  ) * F_Δx_F  ; dropzeros!(U_DIVy_F);

        T_Lap_T   = T_DIVx_U * U_∂x_T + T_DIVy_V * V_∂y_T ; dropzeros!(T_Lap_T);
        U_Lap_U   = filter_U * ( U_DIVx_T * T_∂x_U + U_DIVy_F * F_∂y_U ) ; dropzeros!(U_Lap_U);
        V_Lap_V   = filter_V * ( V_DIVx_F * F_∂x_V + V_DIVy_T * T_∂y_V ) ; dropzeros!(V_Lap_V);

        # Just for test. This is not real divergence and P grid is needed
#        U_Lap_U = filter_U * ( U_∂x_T * T_∂x_U + ( U_invΔy_U * (op.U_S_U - op.U_I_U) - U_invΔy_U * (op.U_I_U - op.U_N_U) ) ) ; dropzeros!(U_Lap_U);
#        V_Lap_V = filter_V * ( V_∂y_T * T_∂y_V + ( V_invΔx_V * (op.V_W_V - op.V_I_V) - V_invΔx_V * (op.V_I_V - op.V_E_V) ) ) ; dropzeros!(V_Lap_V);

        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse

            #=
            @time for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    m_t[:, i] ./= wgt
                end
            end
            =#
            #for (i, wgt) in enumerate(wgts)
            #    if wgt != 0
            #        m[i, :] ./= wgt
            #    end
            #end

            # Hack into sparse matrix to speed up. Speed up like 10 times at least
            for (i, wgt) in enumerate(wgts)
                if wgt == 1
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] .= 0
                elseif wgt > 1
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] ./= wgt
                end
            end
          
            return transpose(m_t) |> sparse
     
        end

        ones_U = ones(Float64, op.U_pts)
        ones_V = ones(Float64, op.V_pts)
        ones_T = ones(Float64, op.T_pts)

        U_interp_V = (op.U_SW_V + op.U_SE_V + op.U_NW_V + op.U_NE_V) * filter_V     
        U_interp_V = selfDivision(U_interp_V, ones_V)

        V_interp_U = (op.V_SW_U + op.V_SE_U + op.V_NW_U + op.V_NE_U) * filter_U
        V_interp_U = selfDivision(V_interp_U, ones_U)
       
        V_interp_T = (op.V_S_T + op.V_N_T) * filter_T
        V_interp_T = selfDivision(V_interp_T, ones_T)

        U_interp_T = (op.U_W_T + op.U_E_T) * filter_T
        U_interp_T = selfDivision(U_interp_T, ones_T)
 
        T_interp_U = (op.T_W_U + op.T_E_U) * filter_U
        T_interp_U = selfDivision(T_interp_U, ones_U)
 
        T_interp_V = (op.T_S_V + op.T_N_V) * filter_V
        T_interp_V = selfDivision(T_interp_V, ones_V)
 
        #println("Making coriolis operators")
        f = gi.c_f |> cvt23

        # Dangerous. Coriolis force can result in energy leak.
        #Need to deal with it carefully
#        U_f_U = filter_U * spdiagm( 0 => view(U_interp_T * view(f, :), :) )
#        V_f_V = filter_V * spdiagm( 0 => view(V_interp_T * view(f, :), :) )

        T_f_T = filter_T * spdiagm( 0 => view(f, :) )
        
        # imagine term fv act on U grid
        # filter_U and filter_V have already been applied in U_f_U and V_f_V
#        U_f_V = filter_U * U_f_U * U_interp_V
#        V_f_U = filter_V * V_f_V * V_interp_U
 
#        U_f_V = filter_U * U_interp_V * V_f_V
#        V_f_U = filter_V * V_interp_U * U_f_U
 

        U_f_V = filter_U * U_interp_T * T_f_T * T_interp_V
        V_f_U = filter_V * V_interp_T * T_f_T * T_interp_U

       


        


        dropzeros!(U_∂x_U);

       # dropzeros!(U_f_U)
       # dropzeros!(V_f_V)
        dropzeros!(U_f_V)
        dropzeros!(V_f_U)

        return new(
            op,
            U_∂x_T, V_∂y_T,
            U_∂x_U, U_∂y_U,
            V_∂x_V, V_∂y_V,
            T_DIVx_U, T_DIVy_V,
            T_Lap_T, U_Lap_U, V_Lap_V,
            U_interp_V, V_interp_U,
            T_interp_U, T_interp_V,
            U_interp_T, V_interp_T,
            U_f_V, V_f_U,
            filter_T, filter_U, filter_V, filter_F,
            borderfilter_T,
        )
    end
end
