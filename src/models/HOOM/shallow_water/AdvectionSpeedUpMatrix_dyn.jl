
using SparseArrays

@inline function flat_i(
    k :: Int64,
    i :: Int64,
    j :: Int64,
    Nz :: Int64,
    Nx :: Int64,
    Ny :: Int64,
)
    return k + (i-1) * Nz + (j-1) * Nz * Nx
end

@inline function cyc(i::Int64, N::Int64)
    return mod(i-1, N) + 1
end


mutable struct DynamicAdvSpeedUpMatrix

    DIV_X    :: AbstractArray{Float64, 2}
    DIV_Y    :: AbstractArray{Float64, 2}

    # notation : T grid to U grid written as UT  (matrix arrangement)

    dx_UT   :: AbstractArray{Float64, 2}   # ∇b 
    dy_VT   :: AbstractArray{Float64, 2}   
 
    dx_UU   :: AbstractArray{Float64, 2}   # ∇u
    dy_UU   :: AbstractArray{Float64, 2}

    dx_VV   :: AbstractArray{Float64, 2}   # ∇v
    dy_VV   :: AbstractArray{Float64, 2}

    interp_VU :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid
    interp_UV :: AbstractArray{Float64, 2}  # interpolation of V grid onto U grid

    f_U2V      :: AbstractArray{Float64, 2}   # used to get fu on V grid
    f_V2U      :: AbstractArray{Float64, 2}   # used to get fv on U grid

    function AdvectionSpeedUpMatrix(;
        gi             :: DisplacedPoleCoordinate.GridInfo,
        Nx             :: Int64,
        Ny             :: Int64,
        Nz             :: Int64,
        mask2          :: AbstractArray{Float64, 2},
        noflux_x_mask2 :: AbstractArray{Float64, 2},
        noflux_y_mask2 :: AbstractArray{Float64, 2},
    )

        elm_max = (Nz+1)*(Nx+1)*(Ny+1) * 2 
        I = zeros(Int64, elm_max)
        J = zeros(Int64, elm_max)
        V = zeros(Float64, elm_max)
        idx = 0

        function add!(i::Int64, j::Int64, v::Float64)
            idx += 1
            I[idx] = i
            J[idx] = j
            V[idx] = v
        end

        function getSparse!(m::Int64, n::Int64)
            s = sparse(view(I, 1:idx), view(J, 1:idx), view(V, 1:idx), m, n)
            idx = 0
            return s 
        end


        speye = (dtype, n) -> spdiagm(0=>ones(dtype, n))

        elm_max = (Nz+1)*(Nx+1)*(Ny+1) * 2 

        # Making operator
        U_pts = (Nx+1) * Ny * Nz
        V_pts = Nx * (Ny+1) * Nz
        T_pts = Nx * Ny * Nz

        # I_U , Ie_U
        I_UU = speye(Float64, U_pts)
        I_VV = speye(Float64, V_pts)
        I_TT = speye(Float64, T_pts)

        num_U = zeros(Int64, Nz, Nx+1, Ny)
        num_V = zeros(Int64, Nz, Nx, Ny+1)
        num_T = zeros(Int64, Nz, Nx,   Ny)

        num_U[:] = 1:length(num_U)
        num_V[:] = 1:length(num_V)
        num_T[:] = 1:length(num_T)
        
        U = num_U * 0
        V = num_V * 0
        T = num_T * 0

        U_flat = view(U, :)
        V_flat = view(V, :)
        T_flat = view(T, :)

        # East West Identity Matrix
        U[:, 1:end-1, :] = num_T; U[:, end, :] = view(U, :, 1, :);    E_UT = I_TT[U_flat, :]; dropzeros!(E_UT); U .= 0
        U[:, 2:end, :]   = num_T; U[:, 1,   :] = view(U, :, Nx+1, :); W_UT = I_TT[U_flat, :]; dropzeros!(W_UT); U .= 0

        U[:] = circshift(num_U, (0, -1, 0)); E_UU = I_UU[U_flat, :]; dropzeros!(E_UU); U .= 0
        U[:] = circshift(num_U, (0,  1, 0)); W_UU = I_UU[U_flat, :]; dropzeros!(W_UU); U .= 0

        V[:] = circshift(num_V, (0, -1, 0)); E_VV = I_VV[V_flat, :]; dropzeros!(E_VV); V .= 0
        V[:] = circshift(num_V, (0,  1, 0)); W_VV = I_VV[V_flat, :]; dropzeros!(W_VV); V .= 0
 
        # North South Identity Matrix
        V[:, :, 1:Ny]   = num_T; V[:, :, Ny+1] = view(V, :, :, Ny);  N_VT = I_TT[V_flat, :]; dropzeros!(N_VT); V .= 0
        V[:, :, 2:Ny+1] = num_T; V[:, :,    1] = view(V, :, :,  2);  S_VT = I_TT[V_flat, :]; dropzeros!(S_VT); V .= 0

        U[:, :, 1:Ny-1] = view(num_U, :, :, 2:Ny); U[:, :, Ny] = view(U, :, :, Ny-1);  N_UU = I_UU[U_flat, :]; dropzeros!(N_UU); U .= 0
        U[:, :, 2:Ny] = view(num_U, :, :, 1:Ny-1); U[:, :, 1] = view(U, :, :, 2);      S_UU = I_UU[U_flat, :]; dropzeros!(S_UU); U .= 0
        
        V[:, :, 1:Ny]   = view(num_V, :, :, 2:Ny+1); V[:, :, Ny+1] = view(V, :, :, Ny);  N_VV = I_VV[V_flat, :]; dropzeros!(N_VV); V .= 0
        V[:, :, 2:Ny+1] = view(num_V, :, :, 1:Ny  ); V[:, :,    1] = view(V, :, :,  2);  S_VV = I_VV[V_flat, :]; dropzeros!(S_VV); V .= 0


        inv_dx_UU = spdiam( 0 => gi.)

        # MAGIC!!
        dx_UT = inv_dx_UU * (E_UT - W_UT)
        dy_VT = inv_dy_VV * (N_VT - S_VT)

        dx_UU = inv_dx_UU * (E_UU - W_UU) / 2.0
        dy_UU = inv_dy_UU * (N_UU - S_UU) / 2.0

        dx_VV = inv_dx_VV * (E_VV - W_VV) / 2.0 
        dy_VV = inv_dy_VV * (N_VV - S_VV) / 2.0 

        interp_UV = E_UV
        interp_VU = E_VU
        
    f_U2V      :: AbstractArray{Float64, 2}   # used to get fu on V grid
    f_V2U      :: AbstractArray{Float64, 2}   # used to get fv on U grid





        # ===== [BEGIN] Making interp matrix =====
        # x
        for i=1:Nx+1, j=1:Ny  # iterate through bounds
            for k=1:Nz[cyc(i, Nx), j]  # Bounds Nx+1 is the same as the bound 1
                if noflux_x_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i           , j, Nz, Nx+1, Ny)
                    ic_e = flat_i(k, cyc(i  ,Nx) , j, Nz, Nx  , Ny)
                    ic_w = flat_i(k, cyc(i-1,Nx) , j, Nz, Nx  , Ny)

                    #u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
                    add!(ib, ic_w, 1.0 - gi.weight_e[i, j])
                    add!(ib, ic_e, gi.weight_e[i, j])
                end
            end
        end
        mtx_interp_U = getSparse!(Nz * (Nx+1) * Ny, Nz * Nx * Ny)

        # y
        for i=1:Nx, j=2:Ny   # iterate through bounds
            for k=1:Nz[i, j]
               if noflux_y_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i, j  , Nz, Nx, Ny+1)
                    ic_n = flat_i(k, i, j  , Nz, Nx, Ny  )
                    ic_s = flat_i(k, i, j-1, Nz, Nx, Ny  )

                    #v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
                    add!(ib, ic_s, 1.0 - gi.weight_n[i, j])
                    add!(ib, ic_n, gi.weight_n[i, j])

                end

            end
        end
        mtx_interp_V = getSparse!(Nz * Nx * (Ny+1), Nz * Nx * Ny)

#=
        for i=1:Nx+1, j=1:Ny  # iterate through bounds
            for k=1:Nz[cyc(i, Nx), j]  # Bounds Nx+1 is the same as the bound 1
                if noflux_x_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i           , j, Nz, Nx+1, Ny)
                    ic_e = flat_i(k, cyc(i  ,Nx) , j, Nz, Nx  , Ny)
                    ic_w = flat_i(k, cyc(i-1,Nx) , j, Nz, Nx  , Ny)

                    #u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
                    mtx_interp_U[ic_w, ib] = 1.0 - gi.weight_e[i, j] 
                    mtx_interp_U[ic_e, ib] = gi.weight_e[i, j]
                end
            end
        end

        # y
        for i=1:Nx, j=2:Ny   # iterate through bounds
            for k=1:Nz[i, j]
               if noflux_y_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i, j  , Nz, Nx, Ny+1)
                    ic_n = flat_i(k, i, j  , Nz, Nx, Ny  )
                    ic_s = flat_i(k, i, j-1, Nz, Nx, Ny  )

                    #v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
                    mtx_interp_V[ic_s, ib] = 1.0 - gi.weight_n[i, j] 
                    mtx_interp_V[ic_n, ib] = gi.weight_n[i, j]
                end

            end
        end
=#
        # ===== [END] Making interp matrix =====

        println("Making Divergence Matrix")
        # ===== [BEGIN] Making divergent matrix =====

        # x
        for i=1:Nx, j=1:Ny  # iterate through face centers
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # X direction
                ib_e   = flat_i(k, i+1, j, Nz, Nx+1, Ny)
                ib_w   = flat_i(k, i  , j, Nz, Nx+1, Ny)

                add!(ic, ib_e,   gi.DY[i+1, j] / gi.dσ[i, j])
                add!(ic, ib_w, - gi.DY[i  , j] / gi.dσ[i, j])

            end
        end
        mtx_DIV_X = getSparse!(Nz * Nx * Ny, Nz * (Nx+1) * Ny)

        # y
        for i=1:Nx, j=1:Ny  # iterate through face centers
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # Y direction
                ib_n   = flat_i(k, i, j+1, Nz, Nx, Ny+1)
                ib_s   = flat_i(k, i, j  , Nz, Nx, Ny+1)


                #add!(ic, ib_n,   gi.DX[i, j] / gi.dσ[i, j])
                add!(ic, ib_n,   gi.DX[i, j+1] / gi.dσ[i, j])
                add!(ic, ib_s, - gi.DX[i, j  ] / gi.dσ[i, j])

                #add!(ic, ib_n, 0.0)#   gi.DX[i, j+1] / gi.dσ[i, j])
                #add!(ic, ib_s, 0.0)#- gi.DX[i, j  ] / gi.dσ[i, j])


            end
        end
        mtx_DIV_Y = getSparse!(Nz * Nx * Ny, Nz * Nx * (Ny+1))

        # z
        for i=1:Nx, j=1:Ny  # iterate through face centers
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # Z direction
                ib_t   = flat_i(k  , i, j, Nz+1, Nx, Ny)
                ib_b   = flat_i(k+1, i, j, Nz+1, Nx, Ny)

                add!(ic, ib_t,   1.0 / hs[k, i, j])
                add!(ic, ib_b, - 1.0 / hs[k, i, j])

            end
        end
        mtx_DIV_Z = getSparse!(Nz * Nx * Ny, (Nz+1) * Nx * Ny)
        # ===== [END] Making divergent matrix =====
        
        println("Making GRAD Matrix")
        # ===== [BEGIN] Making GRAD matrix =====
        # x
        for i=1:Nx+1, j=1:Ny  # iterate through bounds
            for k=1:Nz[cyc(i, Nx), j]  # Bounds Nx+1 is the same as the bound 1
                if noflux_x_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i           , j, Nz, Nx+1, Ny)
                    ic_e = flat_i(k, cyc(i  ,Nx) , j, Nz, Nx  , Ny)
                    ic_w = flat_i(k, cyc(i-1,Nx) , j, Nz, Nx  , Ny)
                    
                    # ( qs[k, i, j] - qs[k, i-1, j] ) / gi.dx_w[i, j] 
                    add!(ib, ic_e,   1.0 / gi.dx_w[cyc(i, Nx), j])
                    add!(ib, ic_w, - 1.0 / gi.dx_w[cyc(i, Nx), j])
                end
            end
        end
        mtx_GRAD_X = getSparse!(Nz * (Nx+1) * Ny, Nz * Nx * Ny)

        # y
        for i=1:Nx, j=2:Ny   # iterate through bounds
            for k=1:Nz[i, j]
               if noflux_y_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i, j  , Nz, Nx, Ny+1)
                    ic_n = flat_i(k, i, j  , Nz, Nx, Ny  )
                    ic_s = flat_i(k, i, j-1, Nz, Nx, Ny  )

                    # ( qs[k, i, j] - qs[k, i, j-1] ) / gi.dy_s[i, j]
                    add!(ib, ic_n,   1.0 / gi.dy_s[i, j])
                    add!(ib, ic_s, - 1.0 / gi.dy_s[i, j])
                end
            end
        end
        mtx_GRAD_Y = getSparse!(Nz * Nx * (Ny+1), Nz * Nx * Ny)

#        println("GRADXY")
#        println(mtx_GRAD_Y)
#        println(mtx_GRAD_X)


        # z
        for i=1:Nx, j=1:Ny   # iterate through bounds

            if mask3[1, i, j] == 0.0
                continue
            end
        
            _Nz = Nz[i, j]

            # The frist layer of w is zero -- means do not assign any value to this row
 
            # Assign from the second row
            for k=2:_Nz
                ib   = flat_i(k  , i, j, Nz+1, Nx, Ny)
                ic_t = flat_i(k-1, i, j, Nz  , Nx, Ny)
                ic_b = flat_i(k  , i, j, Nz  , Nx, Ny)

                # ( qs[k-1, i, j] - qs[k, i, j] ) / Δzs[k-1, i, j]
                add!(ib, ic_t,   1.0 / Δzs[k-1, i, j])
                add!(ib, ic_b, - 1.0 / Δzs[k-1, i, j])


                # Bottom is the same as the layer right above it.
                # This is necessary because the thickness of last layer might be
                # very thin due to topography to create instability during
                # doing adveciton.
                if k == _Nz
                    ib_b = flat_i(k+1, i, j, Nz+1, Nx, Ny)
                    add!(ib_b, ic_t,   1.0 / Δzs[k-1, i, j])
                    add!(ib_b, ic_b, - 1.0 / Δzs[k-1, i, j])
                end
            end

        end
        mtx_GRAD_Z = getSparse!((Nz+1) * Nx * Ny, Nz * Nx * Ny)
        # ===== [END] Making GRAD matrix =====
        
        println("Making CURV Matrix")
        # ===== [BEGIN] Making CURV matrix =====

        # x
        for i=1:Nx, j=1:Ny
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # X direction
                # CURV_x[k, i, j] = ( GRAD_bnd_x[k, i+1, j  ] - GRAD_bnd_x[k  , i, j] ) / gi.dx_c[i, j]
                ib_e   = flat_i(k, i+1, j, Nz, Nx+1, Ny)
                ib_w   = flat_i(k, i  , j, Nz, Nx+1, Ny)

                add!(ic, ib_e,   1.0 / gi.dx_c[i, j])
                add!(ic, ib_w, - 1.0 / gi.dx_c[i, j])

            end
        end
        mtx_CURV_X = getSparse!(Nz * Nx * Ny, Nz * (Nx+1) * Ny)

        for i=1:Nx, j=1:Ny
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # Y direction
                # CURV_y[k, i, j] = ( GRAD_bnd_y[k, i  , j+1] - GRAD_bnd_y[k  , i, j] ) / gi.dy_c[i, j]
                ib_n   = flat_i(k, i, j+1, Nz, Nx, Ny+1)
                ib_s   = flat_i(k, i, j  , Nz, Nx, Ny+1)

                add!(ic, ib_n,   1.0 / gi.dy_c[i, j])
                add!(ic, ib_s, - 1.0 / gi.dy_c[i, j])

            end
        end
        mtx_CURV_Y = getSparse!(Nz * Nx * Ny, Nz * Nx * (Ny+1))


        for i=1:Nx, j=1:Ny
            for k=1:Nz[i, j]
                if mask3[k, i, j] == 0.0
                    break
                end

                ic = flat_i(k, i, j, Nz, Nx  , Ny)

                # Z direction
                # CURV_z[k, i, j] = ( GRAD_bnd_z[k, i  , j  ] - GRAD_bnd_z[k+1, i, j] ) / hs[k, i, j]
                ib_t   = flat_i(k  , i, j, Nz+1, Nx, Ny)
                ib_b   = flat_i(k+1, i, j, Nz+1, Nx, Ny)

                add!(ic, ib_t,   1.0 / hs[k, i, j])
                add!(ic, ib_b, - 1.0 / hs[k, i, j])
            end
        end
        mtx_CURV_Z = getSparse!(Nz * Nx * Ny, (Nz+1) * Nx * Ny)
        # ===== [END] Making CURV matrix =====
        =#
        return new(
            mtx_interp_U,
            mtx_interp_V,
            mtx_DIV_X, 
            mtx_DIV_Y,  
            mtx_DIV_Z,  
            mtx_GRAD_X,  
            mtx_GRAD_Y,  
            mtx_GRAD_Z,  
            mtx_CURV_X,  
            mtx_CURV_Y,  
            mtx_CURV_Z,  
        )
    end
end
