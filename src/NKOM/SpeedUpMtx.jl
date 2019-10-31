
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


mutable struct AdvectionSpeedUpMatrix
#    mtx_GRAD :: AbstractArray{Float64, 2}
#    mtx_CURV :: AbstractArray{Float64, 2}
    mtx_interp_U :: AbstractArray{Float64, 2}
    mtx_interp_V :: AbstractArray{Float64, 2}
    mtx_DIV2_X   :: AbstractArray{Float64, 2}
    mtx_DIV2_Y   :: AbstractArray{Float64, 2}

#    mtx_DIV3 :: AbstractArray{Float64, 2}

    function AdvectionSpeedUpMatrix(;
        gi             :: DisplacedPoleCoordinate.GridInfo,
        Nx             :: Int64,
        Ny             :: Int64,
        Nz_bone        :: Int64,
        Nz             :: AbstractArray{Int64, 2},
        mask3          :: AbstractArray{Float64, 3},
        noflux_x_mask3 :: AbstractArray{Float64, 3},
        noflux_y_mask3 :: AbstractArray{Float64, 3},
        Δzs            :: AbstractArray{Float64, 3},
        hs             :: AbstractArray{Float64, 3},
    )
   # mtx_GRAD_x = spzeros(Float64, Nx * Ny * Nz_bone, (Nx+1) * Ny * Nz_bone)
   # mtx_GRAD_y = spzeros(Float64, Nx * Ny * Nz_bone, Nx * (Ny+1) * Nz_bone)
   # mtx_GRAD_z = spzeros(Float64, Nx * Ny * Nz_bone, Nx * Ny * (Nz_bone+1))
    
   # mtx_CURV :: AbstractArray{Float64, 2}
     
        mtx_interp_U = spzeros(Float64, Nz_bone * (Nx+1) * Ny    , Nz_bone * Nx * Ny)
        mtx_interp_V = spzeros(Float64, Nz_bone * Nx     * (Ny+1), Nz_bone * Nx * Ny)
        mtx_DIV_U    = spzeros(Float64, Nz_bone * Nx * Ny        , Nz_bone * (Nx+1) * Ny)
        mtx_DIV_V    = spzeros(Float64, Nz_bone * Nx * Ny        , Nz_bone * Nx * (Ny+1))

        println("Size mtx_interp_U: ", size(mtx_interp_U))
        println("Size mtx_interp_V: ", size(mtx_interp_V))

        # ===== [BEGIN] Making interp matrix =====
        # x
        for i=1:Nx+1, j=1:Ny
            for k=1:Nz[cyc(i, Nx), j]
                if noflux_x_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i           , j, Nz_bone, Nx+1, Ny)
                    ic_e = flat_i(k, cyc(i  ,Nx) , j, Nz_bone, Nx  , Ny)
                    ic_w = flat_i(k, cyc(i-1,Nx) , j, Nz_bone, Nx  , Ny)

                    #u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
                    mtx_interp_U[ib, ic_w] = 1.0 - gi.weight_e[i, j] 
                    mtx_interp_U[ib, ic_e] = gi.weight_e[i, j]
                end
            end
        end

        # y
        for i=1:Nx, j=2:Ny
            for k=1:Nz[i, j]
               if noflux_y_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i, j  , Nz_bone, Nx, Ny+1)
                    ic_n = flat_i(k, i, j  , Nz_bone, Nx, Ny  )
                    ic_s = flat_i(k, i, j-1, Nz_bone, Nx, Ny  )

                    #v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
                    mtx_interp_V[ib, ic_s] = 1.0 - gi.weight_n[i, j] 
                    mtx_interp_V[ib, ic_n] = gi.weight_n[i, j]
                end

            end
        end
        # ===== [END] Making interp matrix =====

        # ===== [BEGIN] Making divergent matrix =====
        # x
        for i=1:Nx, j=1:Ny
            for k=1:Nz[cyc(i, Nx), j]
                if mask3[k, i, j] == 0.0
                    break
                end
                ic = flat_i(k, i, j, Nz_bone, Nx  , Ny)
                ib_w   = flat_i(k, i  , j, Nz_bone, Nx+1, Ny)
                ib_e   = flat_i(k, i+1, j, Nz_bone, Nx+1, Ny)

                    #u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
                    mtx_interp_U[ib, ic_w] = 1.0 - gi.weight_e[i, j] 
                    mtx_interp_U[ib, ic_e] = gi.weight_e[i, j]
                end
            end
        end

        # y
        for i=1:Nx, j=2:Ny
            for k=1:Nz[i, j]
               if noflux_y_mask3[k, i, j] != 0.0
                    ib   = flat_i(k, i, j  , Nz_bone, Nx, Ny+1)
                    ic_n = flat_i(k, i, j  , Nz_bone, Nx, Ny  )
                    ic_s = flat_i(k, i, j-1, Nz_bone, Nx, Ny  )

                    #v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
                    mtx_interp_V[ib, ic_s] = 1.0 - gi.weight_n[i, j] 
                    mtx_interp_V[ib, ic_n] = gi.weight_n[i, j]
                end

            end
        end

        # 
#    local tmp = tmp_σ = 0.0
    for i=1:Nx, j=1:Ny

        w_bnd[1, i, j] = 0.0

        for k=1:Nz[i, j]

            if mask3[k, i, j] == 0.0
                break
            end
            
            div[k, i, j] =  (  
                u_bnd[k, i+1, j  ]  * gi.DY[i+1, j  ]
              - u_bnd[k, i,   j  ]  * gi.DY[i  , j  ]
              + v_bnd[k, i,   j+1]  * gi.DX[i  , j+1]
              - v_bnd[k, i,   j  ]  * gi.DX[i  , j  ]
            ) / gi.dσ[i, j]

            w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
        end

#        tmp   += w_bnd[Nz[i, j]+1, i, j] * gi.dσ[i, j]
#        tmp_σ += gi.dσ[i, j]
    end

#    println("tmp: ", tmp, "; tmp_σ: ", tmp_σ, "; Average w: ", tmp/tmp_σ)


        # ===== [END] Making divergent matrix =====
        

        return new(
            mtx_interp_U,
            mtx_interp_V,
        )
    end
end
