using LinearAlgebra: mul!

function calWeightedQuantity(;
    top     :: Float64,
    bot     :: Float64,
    split_z :: Float64,
    zs      :: AbstractArray{Float64, 1},
    hs      :: AbstractArray{Float64, 1},
    layer   :: Integer,
)

    Δh     = ocn.hs[layer, i, j]
    Δh_top = ocn.zs[layer, i, j] - split_z
    Δh_bot = Δh - Δh_top

    return ( Δh_top * top + Δh_bot * bot ) / Δh

end


function stepOcean_prepare!(ocn::Ocean; cfgs...)

    adv_scheme = cfgs[:adv_scheme]

    if adv_scheme == :static
        return
    end

    # Transform input wind stress vector first
    DisplacedPoleCoordinate.project!(ocn.gi, ocn.in_flds.taux, ocn.in_flds.tauy, ocn.τx, ocn.τy, direction=:Forward)

    if adv_scheme == :ekman_all_in_ML
        
        @loop_hor ocn i j let
            
            ϵ = ocn.ϵs[i, j]
            f = ocn.fs[i, j]
            
            τx = ocn.τx[i, j]
            τy = ocn.τy[i, j]
            
            h_ML = ocn.h_ML[i, j]
            Nz   = ocn.Nz[i, j] 
            s2ρh = ρ * h_ML * (ϵ^2.0 + f^2.0)
            
            ek_u = (ϵ * τx + f * τy) / s2ρh
            ek_v = (ϵ * τy - f * τx) / s2ρh

            FLDO = ocn.FLDO[i, j]

            if FLDO == -1
                ocn.u[:, i, j] .= ek_u
                ocn.v[:, i, j] .= ek_v
            else
                ocn.u[1:FLDO-1, i, j] .= ek_u
                ocn.v[1:FLDO-1, i, j] .= ek_v

                ocn.u[FLDO:Nz, i, j] .= 0.0
                ocn.v[FLDO:Nz, i, j] .= 0.0
            end
            
        end
    elseif adv_scheme == :testu
 
        @loop_hor ocn i j let
            for k = 1:ocn.Nz[i, j]
                ocn.u[k, i, j] = 2π/365/86400.0 * cos(ocn.gi.c_lat[i, j]) * Re
                ocn.v[k, i, j] = 0.0
            end
        end

        #ocn.u .= 0.0
        #ocn.v[1, :, :] .= 0.1
        
    elseif adv_scheme == :ekman_codron2012_partition

        H_ek =  50.0
        H_rf = 250.0   # Codron (2012) suggests 150 - 350 meters. Here I take the average.
        H_total = H_ek + H_rf

        bot_lay_ek = getLayerFromDepth(
            z  = - H_ek,
            zs = ocn.zs_bone,  
            Nz = ocn.Nz_bone,
        )

        bot_lay_rf = getLayerFromDepth(
            z  = - H_total,
            zs = ocn.zs_bone,  
            Nz = ocn.Nz_bone,
        )

        @loop_hor ocn i j let

            M̃ = (ocn.τx[i, j] + ocn.τy[i, j] * im) / (ρ * (ocn.ϵs[i, j] + ocn.fs[i, j] * im) )

            ṽ_ek =   M̃ / H_ek
            ṽ_rf = - M̃ / H_rf

            u_ek, v_ek = real(ṽ_ek), imag(ṽ_ek)
            u_rf, v_rf = real(ṽ_rf), imag(ṽ_rf)

            if bot_lay_ek == -1
            
                ocn.u[:, i, j] .= u_ek
                ocn.v[:, i, j] .= v_ek

            else

                ocn.u[1:bot_lay_ek, i, j] .= u_ek
                ocn.v[1:bot_lay_ek, i, j] .= v_ek

                # Mix the top of RF layer
                Δh     = ocn.hs[bot_lay_ek, i, j]
                Δh_top = H_ek + ocn.zs[bot_lay_ek, i, j]
                Δh_bot = Δh - Δh_top

                ocn.u[bot_lay_ek, i, j] = (Δh_top * u_ek + Δh_bot * u_rf) / Δh
                ocn.v[bot_lay_ek, i, j] = (Δh_top * v_ek + Δh_bot * v_rf) / Δh

                if bot_lay_ek < ocn.Nz[i, j] # Bottom layers exists
                    if bot_lay_rf == -1
                       ocn.u[bot_lay_ek+1:end, i, j] .= u_rf
                       ocn.v[bot_lay_ek+1:end, i, j] .= v_rf
                    else
                       ocn.u[bot_lay_ek+1:bot_lay_rf, i, j] .= u_rf
                       ocn.v[bot_lay_ek+1:bot_lay_rf, i, j] .= v_rf

                        # Mix the bottom of RF layer
                        Δh     = ocn.hs[bot_lay_rf, i, j]
                        Δh_top = H_total + ocn.zs[bot_lay_rf, i, j]
                        Δh_bot = Δh - Δh_top

                        ocn.u[bot_lay_rf, i, j] = Δh_top * u_rf / Δh
                        ocn.v[bot_lay_rf, i, j] = Δh_top * v_rf / Δh

                    end
                end

            end
        end

    elseif adv_scheme == :ekman_simple_partition
        @loop_hor ocn i j let

            h_ML = ocn.h_ML[i, j]
            s̃ = ocn.ϵs[i, j] + ocn.fs[i, j] * im
            H̃ = √(ocn.K_v / s̃)
            H = abs(H̃)
            
            M̃ = (ocn.τx[i, j] + ocn.τy[i, j] * im) / (ρ * s̃)

            depth_rf = 500.0

            H_ek = max(h_ML, 2*H)
            H_rf = max(depth_rf - H_ek, 100.0 )

            ṽ_ek =   M̃ / H_ek
            ṽ_rf = - M̃ / H_rf


            u_ek, v_ek = real(ṽ_ek), imag(ṽ_ek)
            u_rf, v_rf = real(ṽ_rf), imag(ṽ_rf)


            bot_lay_ek = getLayerFromDepth(
                z  = - H_ek,
                zs = ocn.cols.zs[i, j],  
                Nz = ocn.Nz[i, j],
            )

            bot_lay_rf = getLayerFromDepth(
                z  = - H_ek - H_rf,
                zs = ocn.cols.zs[i, j],  
                Nz = ocn.Nz[i, j],
            )

#            if i == 1
#                println(format("[{:03d}] H_ek={:f}, bot_layer_ek={:d}, bot_layer_rf={:d}", j, H_ek, bot_lay_ek, bot_lay_rf))
#                println("zs: ", ocn.cols.zs[i, j][1:10])
#            end
#            if (i, j) == (48, 89)
#                println("u_ek = ", u_ek, "; v_ek=", v_ek)
#                println("u_rf = ", u_rf, "; v_rf=", v_rf)
#            end
            if bot_lay_ek == -1
            
                ocn.u[:, i, j] .= u_ek
                ocn.v[:, i, j] .= v_ek

            else

                ocn.u[1:bot_lay_ek, i, j] .= u_ek
                ocn.v[1:bot_lay_ek, i, j] .= v_ek

                # Mix the top of RF layer
                Δh     = ocn.hs[bot_lay_ek, i, j]
                Δh_top = H_ek + ocn.zs[bot_lay_ek, i, j]
                Δh_bot = Δh - Δh_top

                ocn.u[bot_lay_ek, i, j] = (Δh_top * u_ek + Δh_bot * u_rf) / Δh
                ocn.v[bot_lay_ek, i, j] = (Δh_top * v_ek + Δh_bot * v_rf) / Δh

                if bot_lay_ek < ocn.Nz[i, j] && bot_lay_ek < bot_lay_rf

                    if bot_lay_rf == -1
                       ocn.u[bot_lay_ek+1:end, i, j] .= u_rf
                       ocn.v[bot_lay_ek+1:end, i, j] .= v_rf
                    else
                       ocn.u[bot_lay_ek+1:bot_lay_rf, i, j] .= u_rf
                       ocn.v[bot_lay_ek+1:bot_lay_rf, i, j] .= v_rf

                        # Mix the bottom of RF layer
                        Δh     = ocn.hs[bot_lay_rf, i, j]
                        Δh_top = H_ek + H_rf + ocn.zs[bot_lay_rf, i, j]
                        Δh_bot = Δh - Δh_top

                        ocn.u[bot_lay_rf, i, j] = Δh_top * u_rf / Δh
                        ocn.v[bot_lay_rf, i, j] = Δh_top * v_rf / Δh

                    end
                end

            end
        end

    #=
    elseif adv_scheme == :test
        ocn.u .= 0.0
        ocn.v .= 0.0
    elseif adv_scheme == :testusin
        @loop_hor ocn i j let
            for k = 1:ocn.Nz[i, j]
                ocn.u[k, i, j] = .1 * exp(ocn.zs[k, i, j]/50.0) * sin(ocn.gi.c_lon[i, j])
                ocn.v[k, i, j] = .1 * exp(ocn.zs[k, i, j]/50.0) * cos(ocn.gi.c_lat[i, j])
            end
        end

        #ocn.u .= 0.0
        #ocn.v[1, :, :] .= 0.1
    =#
    else
        throw(ErrorException("Unknown advection scheme: " * string(adv_scheme)))
    end


    #println("calHorVelBnd!")    
    #= 
    calHorVelBnd!(
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        weight_e = ocn.gi.weight_e,
        weight_n = ocn.gi.weight_n,
        u     = ocn.u,
        v     = ocn.v,
        u_bnd = ocn.u_bnd,
        v_bnd = ocn.v_bnd,
        mask3 = ocn.mask3,
        noflux_x_mask3 = ocn.noflux_x_mask3,
        noflux_y_mask3 = ocn.noflux_y_mask3,
    )


    calDIV!(
        gi    = ocn.gi,
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        u_bnd = ocn.u_bnd,
        v_bnd = ocn.v_bnd,
        div   = ocn.div,
        mask3 = ocn.mask3,
    )
    =#

    #println("calHorVelBnd with spmtx")
    mul!(view(ocn.u_bnd, :), ocn.ASUM.mtx_interp_U, view(ocn.u, :))
    mul!(view(ocn.v_bnd, :), ocn.ASUM.mtx_interp_V, view(ocn.v, :))
    mul!(view(ocn.div, :), ocn.ASUM.mtx_DIV_X, view(ocn.u_bnd, :))
    mul!(view(ocn.workspace2, :), ocn.ASUM.mtx_DIV_Y, view(ocn.v_bnd, :))

    ocn.div .+= ocn.workspace2

    calVerVelBnd!(
        gi    = ocn.gi,
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        w_bnd = ocn.w_bnd,
        hs    = ocn.hs,
        div   = ocn.div,
        mask3 = ocn.mask3,
    )
#=
    println("Max nswflx: ", maximum(ocn.in_flds.nswflx[ocn.mask_idx]))
    println("Min nswflx: ", minimum(ocn.in_flds.nswflx[ocn.mask_idx]))
    println("Max swflx: ", maximum(abs.(ocn.in_flds.swflx[ocn.mask_idx])))
    println("Max tau : ", sqrt(maximum(abs.(ocn.in_flds.taux[ocn.mask_idx].^2.0 + ocn.in_flds.tauy[ocn.mask_idx].^2.0))))
    println("Max taux: ", maximum(abs.(ocn.in_flds.taux[ocn.mask_idx])))
    println("Max tauy: ", maximum(abs.(ocn.in_flds.tauy[ocn.mask_idx])))
 
    println("Max τx: ", maximum(abs.(ocn.τx[ocn.mask_idx])))
    println("Max τy: ", maximum(abs.(ocn.τy[ocn.mask_idx])))
    println("Max u: ", maximum(abs.(ocn.u)))
    println("Max v: ", maximum(abs.(ocn.v)))
    println("Max w_bnd: ", maximum(abs.(ocn.w_bnd)))
=#
    #=    
    # Calculate ∇⋅v
    for k=1:ocn.Nz_bone
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.u[k],  ocn.lays.v[k],  ocn.lays.div[k], ocn.lays.mask3[k])
    end

    # Calculate w
    @loop_hor ocn i j let

        Nz = ocn.Nz[i, j]

#=
        ocn.w[1, i, j] = 0.0

#        for k = 2:Nz+1
#            ocn.w[k, i, j] = ocn.w[k-1, i, j] + ocn.div[k-1, i, j]
#        end

        for k = 2:Nz
            ocn.w[k, i, j] = ocn.w[k-1, i, j] + (ocn.hs[k-1, i, j] * ocn.div[k-1, i, j] + ocn.hs[k, i, j] * ocn.div[k, i, j]) / 2.0
        end
=#

        ocn.w_bnd[1, i, j] = 0.0

        for k = 2:Nz+1
            Δw = ocn.hs[k-1, i, j] * ocn.div[k-1, i, j]
            ocn.w_bnd[k, i, j] = ocn.w_bnd[k-1, i, j] + Δw
            ocn.w[k-1, i, j]   = ocn.w_bnd[k-1, i, j] + Δw / 2.0
        end
        
    end

    #ocn.w .= -1e-4
    =#
end

