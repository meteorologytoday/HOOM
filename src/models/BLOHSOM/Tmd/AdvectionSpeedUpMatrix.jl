


mutable struct AdvectionSpeedUpMatrix

    op       :: MatrixOperators

    T_DIVx_U    :: AbstractArray{Float64, 2}
    T_DIVy_V    :: AbstractArray{Float64, 2}
    T_DIVz_W    :: AbstractArray{Float64, 2}
    
    U_∂x_T      :: AbstractArray{Float64, 2}
    V_∂y_T      :: AbstractArray{Float64, 2}
    W_∂z_T      :: AbstractArray{Float64, 2}

    T_∂x_U      :: AbstractArray{Float64, 2}
    T_∂y_V      :: AbstractArray{Float64, 2}
    T_∂z_W      :: AbstractArray{Float64, 2}


    U_interp_T :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    V_interp_T :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    W_interp_T :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid

    filter_T       :: AbstractArray{Float64, 2}
    filter_U       :: AbstractArray{Float64, 2}
    filter_V       :: AbstractArray{Float64, 2}
    filter_W       :: AbstractArray{Float64, 2}

    borderfilter_T :: AbstractArray{Float64, 2}
    
    T_Δvol_T :: AbstractArray{Float64, 2}
    
    Δx_U :: AbstractArray{Float64, 3}
    Δy_V :: AbstractArray{Float64, 3}
    Δz_W :: AbstractArray{Float64, 3}


    function AdvectionSpeedUpMatrix(;
        gi             :: PolelikeCoordinate.GridInfo,
        Nz             :: Int64,
        Nz_av          :: AbstractArray{Int64, 2},
        mask3          :: AbstractArray{Float64, 3},
        noflux_x_mask3 :: AbstractArray{Float64, 3},
        noflux_y_mask3 :: AbstractArray{Float64, 3},
        Δz_T           :: AbstractArray{Float64, 3},      # thickness of T grid
        Δz_W           :: AbstractArray{Float64, 3},      # thickness of W grid
    )


        println("TODO: there should be no horizontal flux in the last z grid (bottom)")
        println("so as to avoid instability caused by diffusion in a thin artial gridbox.")
        println("Possible solution is to do vertical diffusion implicitly.")

        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt23 = (x,) -> repeat(reshape(x, 1, size(x)...), outer=(Nz, 1, 1))
        cvt23_diagm = (x,) -> spdiagm( 0 => view(cvt23(x), :) )
        cvt3_diagm = (x,) -> spdiagm( 0 => view(x, :) )

        Nx = gi.Nx
        Ny = gi.Ny
        
        @time op = MatrixOperators(Nx=Nx, Ny=Ny, Nz=Nz)
         
        mask3_flat = view(mask3,  :)

        onV_if_unblocked_north_onT = op.V_S_T  * mask3_flat
        onV_if_unblocked_south_onT = op.V_N_T  * mask3_flat
        onU_if_unblocked_east_onT  = op.U_W_T  * mask3_flat
        onU_if_unblocked_west_onT  = op.U_E_T  * mask3_flat
        onW_if_unblocked_up_onT    = op.W_DN_T * mask3_flat
        onW_if_unblocked_dn_onT    = op.W_UP_T * mask3_flat

        V_mask = onV_if_unblocked_north_onT .* onV_if_unblocked_south_onT
        U_mask = onU_if_unblocked_east_onT  .* onU_if_unblocked_west_onT
        W_mask = onW_if_unblocked_up_onT    .* onW_if_unblocked_dn_onT

        filter_T = spdiagm(0 => mask3_flat)
        filter_V = spdiagm(0 => V_mask)
        filter_U = spdiagm(0 => U_mask)
        filter_W = spdiagm(0 => W_mask)

        borderfilter_T = spdiagm( 0 => ( 
               (op.T_N_T  * mask3_flat)
            .* (op.T_S_T  * mask3_flat)
            .* (op.T_E_T  * mask3_flat)
            .* (op.T_W_T  * mask3_flat)
            .* (op.T_UP_T * mask3_flat)
            .* (op.T_DN_T * mask3_flat)
        ))

        # ===== [ BEGIN face area and lengths on U V ] =====
        
        Δx_U = gi.dx_w |> cvt23
        Δy_U = gi.DY   |> cvt23
 
        Δx_V = gi.DX   |> cvt23
        Δy_V = (hcat(gi.dy_s, gi.dy_n[:, end])) |> cvt23

        V_Δx_V    = (Δx_V         |>  cvt3_diagm)
        U_Δy_U    = (Δy_U         |>  cvt3_diagm)
        W_Δz_W    = (Δz_W         |>  cvt3_diagm)

        U_invΔx_U = (Δx_U.^(-1)   |>  cvt3_diagm)
        U_invΔy_U = (Δy_U.^(-1)   |>  cvt3_diagm)
        
        V_invΔx_V = (Δx_V.^(-1)   |>  cvt3_diagm)
        V_invΔy_V = (Δy_V.^(-1)   |>  cvt3_diagm)

        W_invΔz_W = (Δz_W.^(-1)   |>  cvt3_diagm)

        Δσ_U = Δx_U .* Δy_U
        Δσ_V = Δx_V .* Δy_V

        U_invΔσ_U = Δσ_U.^(-1) |> cvt3_diagm
        V_invΔσ_V = Δσ_V.^(-1) |> cvt3_diagm

        # ===== [ END face area and lengths on U V ] =====
        
        # ===== [ BEGIN face area and lengths on T ] =====

        Δx_T = gi.dx_c |> cvt23
        Δy_T = gi.dy_c |> cvt23

        Δσ_T = Δx_T .* Δy_T 

        T_Δx_T = (Δx_T |>  cvt3_diagm)
        T_Δy_T = (Δy_T |>  cvt3_diagm)
        T_Δz_T = (Δz_T |>  cvt3_diagm)
        
        T_invΔx_T = (Δx_T.^(-1) |> cvt3_diagm)
        T_invΔy_T = (Δy_T.^(-1) |> cvt3_diagm)
        T_invΔσ_T = (Δσ_T.^(-1) |> cvt3_diagm)

        T_invΔz_T = (Δz_T.^(-1) |>  cvt3_diagm)

        T_Δvol_T = T_Δx_T * T_Δy_T * T_Δz_T

        # Δz is special. Need to clear NaN
        function clearNaN!(m)
            for i = 1:length(m.nzval)
                if isnan(m.nzval[i])
                    m.nzval[i] = 0
                end
            end
            dropzeros!(m)
        end

        clearNaN!(T_Δz_T)
        clearNaN!(W_Δz_W)
        clearNaN!(T_invΔz_T)
        clearNaN!(W_invΔz_W)

        clearNaN!(T_Δvol_T)

#        any(Δz_T .== 0) && throw(ErrorException("0!"))
#        any(Δz_W .== 0) && throw(ErrorException("00!"))
        
#        isnan(sum(T_invΔz_T)) && throw(ErrorException("000!"))
#        isnan(sum(W_invΔz_W)) && throw(ErrorException("0000!"))

        # ===== [ END face area and lengths on T ] =====

        # ===== [ BEG making matrix ] =====
        # MAGIC!!

        T_DIVx_U = filter_T * T_invΔσ_T * ( op.T_W_U - op.T_E_U    ) * U_Δy_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V = filter_T * T_invΔσ_T * ( op.T_S_V - op.T_N_V    ) * V_Δx_V  ; dropzeros!(T_DIVy_V);
        T_DIVz_W = filter_T * T_invΔz_T * ( op.T_DN_W - op.T_UP_W  )           ; dropzeros!(T_DIVz_W);

        U_∂x_T = filter_U * U_invΔx_U * (op.U_W_T  - op.U_E_T)                 ; dropzeros!(U_∂x_T);
        V_∂y_T = filter_V * V_invΔy_V * (op.V_S_T  - op.V_N_T)                 ; dropzeros!(V_∂y_T);
        W_∂z_T = filter_W * W_invΔz_W * (op.W_DN_T - op.W_UP_T)                ; dropzeros!(W_∂z_T);

        T_∂x_U  = filter_T * T_invΔx_T * ( op.T_W_U - op.T_E_U )               ; dropzeros!(T_∂x_U);
        T_∂y_V  = filter_T * T_invΔy_T * ( op.T_S_V - op.T_N_V )               ; dropzeros!(T_∂x_U);
        T_∂z_W  = filter_T * T_invΔz_T * ( op.T_DN_W - op.T_UP_W )               ; dropzeros!(T_∂x_U);



        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse
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

        ones_T = ones(Float64, op.T_pts)

        U_interp_T = (op.U_W_T + op.U_E_T) * filter_T
        U_interp_T = selfDivision(U_interp_T, ones_T)

        V_interp_T = (op.V_S_T + op.V_N_T) * filter_T
        V_interp_T = selfDivision(V_interp_T, ones_T)

        W_interp_T = (op.W_DN_T + op.W_UP_T) * filter_T
        W_interp_T = selfDivision(W_interp_T, ones_T)

        return new(
            op,

            T_DIVx_U,
            T_DIVy_V,
            T_DIVz_W,

            U_∂x_T,
            V_∂y_T,
            W_∂z_T,

            T_∂x_U,
            T_∂y_V,
            T_∂z_W,

            U_interp_T,
            V_interp_T,
            W_interp_T,

            filter_T,
            filter_U,
            filter_V,
            filter_W,

            borderfilter_T,

            T_Δvol_T,

            Δx_U,
            Δy_V,
            Δz_W,


        )
    end
end
