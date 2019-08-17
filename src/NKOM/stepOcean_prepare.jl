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
            p̃ = exp(- ocn.h_ML[i, j] * H̃)
            
            M̃ = (ocn.τx[i, j] + ocn.τy[i, j] * im) / (ρ * s̃)
            M̃_DO = M̃ * p̃

            ṽ_ML = (M̃ - M̃_DO) / h_ML

            u_ML, v_ML = real(ṽ_ML), imag(ṽ_ML)

            FLDO = ocn.FLDO[i, j]

            if FLDO == -1
            
                ocn.u[:, i, j] .= u_ML
                ocn.v[:, i, j] .= v_ML

            else

                ocn.u[:, i, j] .= 0.0
                ocn.v[:, i, j] .= 0.0

                if FLDO > 1
                    ocn.u[1:FLDO-1, i, j] .= u_ML
                    ocn.v[1:FLDO-1, i, j] .= v_ML
                end

                H *= 3
                
                eklayer = getLayerFromDepth(
                    z  = - h_ML - H,
                    zs = ocn.cols.zs[i, j],  
                    Nz = ocn.Nz[i, j],
                )
            
                ṽ_DO = (M̃ - M̃_DO) / H

                if abs(ṽ_DO) > abs(ṽ_ML)
                    # this means the bottome boundary condition is
                    # not realistic. So just set ṽ_DO = 0.0
                    ṽ_DO = 0.0
                    #println("Ekman b.c. break at bottom of (",i,",",j,"). Set flow to zero.")
                end

                u_DO, v_DO = real(ṽ_DO), imag(ṽ_DO)

                Δh     = ocn.hs[FLDO, i, j]
                Δh_top = h_ML + ocn.zs[FLDO, i, j]

                if eklayer == FLDO
                    # H is too thin that it lies entirely in the same layer
                    Δh_bot = H
                else
                    Δh_bot = Δh - Δh_top
                end

                ocn.u[FLDO, i, j] = (Δh_top * u_ML + Δh_bot * u_DO) / Δh
                ocn.v[FLDO, i, j] = (Δh_top * v_ML + Δh_bot * v_DO) / Δh

                if FLDO != ocn.Nz[i, j]

                    if eklayer == -1
                        ocn.u[FLDO+1:end, i, j] .= u_DO
                        ocn.v[FLDO+1:end, i, j] .= v_DO
                    else
                        ocn.u[FLDO+1:eklayer, i, j] .= u_DO
                        ocn.v[FLDO+1:eklayer, i, j] .= v_DO
                    end

                end
                


            end
            
        end

    end
        
    # Calculate ∇⋅v
    for k=1:ocn.Nz_bone
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.u[k],  ocn.lays.v[k],  ocn.lays.div[k], ocn.lays.mask3[k])
    end

    # Calculate w
    @loop_hor ocn i j let

        Nz = ocn.Nz[i, j]
        ocn.w[1, i, j] = 0.0
        for k = 2:Nz+1
            ocn.w[k, i, j] = ocn.w[k-1, i, j] + ocn.div[k-1, i, j]
        end
    end

    
end

