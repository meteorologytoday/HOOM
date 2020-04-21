
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
    borderfilter_T :: AbstractArray{Float64, 2}

    function DynamicAdvSpeedUpMatrix(;
        gi             :: PolelikeCoordinate.GridInfo,
        Nz             :: Int64,
        mask2          :: AbstractArray{Float64, 2},
    )

       #println("Begin")
        Nx = gi.Nx
        Ny = gi.Ny

        println("Build MatrixOperators")
        @time op = MatrixOperators(Nx=Nx, Ny=Ny, Nz=Nz)

        # making filter
        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt23 = (x,) -> repeat(reshape(x, 1, size(x)...), outer=(Nz, 1, 1))
        cvt23_diagm = (x,) -> spdiagm( 0 => view(cvt23(x), :) )

        println("Build filters")
        mask3_flat = view(mask2 |> cvt23, :)
        
        if_mask_north_V = op.V_S_T * mask3_flat
        if_mask_south_V = op.V_N_T * mask3_flat
        if_mask_east_U  = op.U_W_T * mask3_flat
        if_mask_west_U  = op.U_E_T * mask3_flat


        #println(size(op.V_S_T))
        #println(size(mask3_flat))

        filter_T = spdiagm(0 => mask3_flat)
        filter_V = spdiagm(0 => if_mask_north_V .* if_mask_south_V)
        filter_U = spdiagm(0 => if_mask_east_U  .* if_mask_west_U)

        borderfilter_T = spdiagm( 0 => ( 
               (op.T_N_T * mask3_flat)
            .* (op.T_S_T * mask3_flat)
            .* (op.T_E_T * mask3_flat)
            .* (op.T_W_T * mask3_flat)
        ))

        #println("if_mask_north_V", size(if_mask_north_V))
        #println("filter_V", size(filter_V))

        # Some face area and lengths

        U_invΔx_U = (gi.dx_w.^(-1) |> cvt23_diagm)
        U_invΔy_U = (gi.DY.^(-1)   |> cvt23_diagm)
        
        V_invΔx_V = (gi.DX.^(-1)   |> cvt23_diagm)
        V_invΔy_V = (hcat(gi.dy_s, gi.dy_n[:, end]).^(-1) |> cvt23_diagm)

        V_Δx_V    = (gi.DX       |> cvt23_diagm)
        U_Δy_U    = (gi.DY       |> cvt23_diagm)
        T_invΔσ_T = (gi.dσ.^(-1) |> cvt23_diagm)

        println("Making various operators")
        # MAGIC!!
        U_∂x_T = filter_U * U_invΔx_U * (op.U_W_T - op.U_E_T) ; dropzeros!(U_∂x_T);
        V_∂y_T = filter_V * V_invΔy_V * (op.V_S_T - op.V_N_T) ; dropzeros!(V_∂y_T);

        U_∂x_U = filter_U * U_invΔx_U * (op.U_W_U - op.U_E_U) / 2.0 ; dropzeros!(U_∂x_U);
        U_∂y_U = filter_U * U_invΔy_U * (op.U_S_U - op.U_N_U) / 2.0 ; dropzeros!(U_∂y_U);

        V_∂x_V = filter_V * V_invΔx_V * (op.V_W_V - op.V_E_V) / 2.0 ; dropzeros!(V_∂x_V);
        V_∂y_V = filter_V * V_invΔy_V * (op.V_S_V - op.V_N_V) / 2.0 ; dropzeros!(V_∂y_V);
 
        T_DIVx_U = filter_T * T_invΔσ_T * ( op.T_W_U - op.T_E_U  ) * U_Δy_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V = filter_T * T_invΔσ_T * ( op.T_S_V - op.T_N_V  ) * V_Δx_V  ; dropzeros!(T_DIVy_V);


        function selfDivision!(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse

            #=
            @time for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    m_t[:, i] ./= wgt
                end
            end
            =#

            # Hack into sparse matrix to speed up. Speed up like 10 times at least
            @time for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg: _end] ./= wgt
                end
            end
           
     
        end

        ones_U = ones(Float64, op.U_pts)
        ones_V = ones(Float64, op.V_pts)
        ones_T = ones(Float64, op.T_pts)

        U_interp_V = (op.U_SW_V + op.U_SE_V + op.U_NW_V + op.U_NE_V) * filter_V     
        selfDivision!(U_interp_V, ones_V)
        dropzeros!(U_interp_V)

        V_interp_U = (op.V_SW_U + op.V_SE_U + op.V_NW_U + op.V_NE_U) * filter_U
        selfDivision!(V_interp_U, ones_U)
        dropzeros!(V_interp_U)
       
        V_interp_T = (op.V_S_T + op.V_N_T) * filter_T
        selfDivision!(V_interp_T, ones_T)
        dropzeros!(V_interp_T)

        U_interp_T = (op.U_W_T + op.U_E_T) * filter_T
        selfDivision!(U_interp_T, ones_T)
        dropzeros!(U_interp_T)
 
        T_interp_U = (op.T_W_U + op.T_E_U) * filter_U
        selfDivision!(T_interp_U, ones_U)
        dropzeros!(T_interp_U)
 
        T_interp_V = (op.T_S_V + op.T_N_V) * filter_V
        selfDivision!(T_interp_V, ones_V)
        dropzeros!(T_interp_V)
 
        #println("Making coriolis operators")
        f = gi.c_f |> cvt23
        U_f_U = filter_U * spdiagm( 0 => view(U_interp_T * view(f, :), :) )
        V_f_V = filter_V * spdiagm( 0 => view(V_interp_T * view(f, :), :) )
        
        # imagine term fv act on U grid
        # filter_U and filter_V have already been applied in U_f_U and V_f_V
        U_f_V = filter_U * U_f_U * U_interp_V
        V_f_U = filter_V * V_f_V * V_interp_U
        
        dropzeros!(U_f_U)
        dropzeros!(V_f_V)
        dropzeros!(U_f_V)
        dropzeros!(V_f_U)

        return new(
            op,
            U_∂x_T, V_∂y_T,
            U_∂x_U, U_∂y_U,
            V_∂x_V, V_∂y_V,
            T_DIVx_U, T_DIVy_V,
            U_interp_V, V_interp_U,
            T_interp_U, T_interp_V,
            U_interp_T, V_interp_T,
            U_f_V, V_f_U,
            filter_T, filter_U, filter_V,
            borderfilter_T,
        )
    end
end
