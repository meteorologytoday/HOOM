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

#=
            if (i, j) == (67, 57)
                println("M̃ = ", M̃)
                println("H_ek = ", H_ek)
                println("H_rf = ", H_rf)
                println("h_ML = ", h_ML)
                println("ṽ_ek = ", ṽ_ek)
                println("ṽ_rf = ", ṽ_rf)
            end
  =#      
        end

    elseif adv_scheme == :test
        println("HERE")
        ocn.u .= 0.1
        ocn.v .= 0.0
    elseif adv_scheme == :testusin
        @loop_hor ocn i j let
            for k = 1:ocn.Nz[i, j]
                ocn.u[k, i, j] = .1 * exp(ocn.zs[k, i, j]/50.0) * sin(ocn.mi.xc[i, j] * π/180.0)
            end
        end

        ocn.u .= 0.0
        ocn.v[1, :, :] .= 0.1

#    else
#        throw(ErrorException("Unknown advection scheme: " * string(adv_scheme)))
    end


     
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
    )

    calVerVelBnd!(
        gi    = ocn.gi,
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        u_bnd = ocn.u_bnd,
        v_bnd = ocn.v_bnd,
        w_bnd = ocn.w_bnd,
        hs    = ocn.hs,
        div   = ocn.div,
        mask3 = ocn.mask3,
    )
   
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

