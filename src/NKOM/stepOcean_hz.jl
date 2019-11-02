function stepOcean_Flow!(
    ocn  :: Ocean;
    cfgs...
)
    
    adv_scheme = cfgs[:adv_scheme]
    do_convadjust = cfgs[:do_convadjust]
    Δt         = cfgs[:Δt]

    if adv_scheme == :static
        return
    end

    # Determine the temperature / salinity of FLDO layer
    println("mixFLDO")
    @time @loop_hor ocn i j let

        ocn.ΔT[i, j] = mixFLDO!(
            qs   = ocn.cols.Ts[i, j],
            zs   = ocn.cols.zs[i, j],
            hs   = ocn.cols.hs[i, j],
            q_ML = ocn.T_ML[i, j],
            FLDO = ocn.FLDO[i, j],
            FLDO_ratio_top = ocn.FLDO_ratio_top[i, j],
            FLDO_ratio_bot = ocn.FLDO_ratio_bot[i, j],
        )

        ocn.ΔS[i, j] = mixFLDO!(
            qs   = ocn.cols.Ss[i, j],
            zs   = ocn.cols.zs[i, j],
            hs   = ocn.cols.hs[i, j],
            q_ML = ocn.S_ML[i, j],
            FLDO = ocn.FLDO[i, j],
            FLDO_ratio_top = ocn.FLDO_ratio_top[i, j],
            FLDO_ratio_bot = ocn.FLDO_ratio_bot[i, j],
        )

    end

    # Pseudo code
    # 1. assign velocity field
    # 2. calculate temperature & salinity flux
    # 3. calculate temperature & salinity flux divergence
    # Gov eqn adv + diff: ∂T/∂t = - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )

#=
    if any(ocn.Ts .>= 100) 
        throw(ErrorException("Ts > 100"))
    end

    if any(ocn.Ss .>= 100) 
        throw(ErrorException("Ss > 100"))
    end
=#

    println("##### QUICK NO SPEEDUP")
    @time calDiffAdv_QUICK!(
        ocn,
        qs          = ocn.Ts,
        wq_bnd      = ocn.wT_bot,
        dΔqdt       = ocn.dΔTdt,
        FLUX_CONV   = ocn.TFLUX_CONV,
        FLUX_CONV_h = ocn.TFLUX_CONV_h,
        FLUX_DEN_x  = ocn.TFLUX_DEN_x,
        FLUX_DEN_y  = ocn.TFLUX_DEN_y,
        FLUX_DEN_z  = ocn.TFLUX_DEN_z,
        Dh = ocn.Dh_T,
        Dv = ocn.Dv_T,
   )
#= 
    println("##### QUICK SPEEDUP")
    @time calDiffAdv_QUICK_SpeedUp!(
        ocn,
        qs          = ocn.Ts,
        wq_bnd      = ocn.wT_bot,
        dΔqdt       = ocn.dΔTdt,
        FLUX_CONV   = ocn.TFLUX_CONV,
        FLUX_CONV_h = ocn.TFLUX_CONV_h,
        FLUX_DEN_x  = ocn.TFLUX_DEN_x,
        FLUX_DEN_y  = ocn.TFLUX_DEN_y,
        FLUX_DEN_z  = ocn.TFLUX_DEN_z,
        Dh = ocn.Dh_T,
        Dv = ocn.Dv_T,
    )
  =#  
    calDiffAdv_QUICK!(
        ocn,
        qs          = ocn.Ss,
        wq_bnd      = ocn.wS_bot,
        dΔqdt       = ocn.dΔSdt,
        FLUX_CONV   = ocn.SFLUX_CONV,
        FLUX_CONV_h = ocn.SFLUX_CONV_h,
        FLUX_DEN_x  = ocn.SFLUX_DEN_x,
        FLUX_DEN_y  = ocn.SFLUX_DEN_y,
        FLUX_DEN_z  = ocn.SFLUX_DEN_z,
        Dh = ocn.Dh_S,
        Dv = ocn.Dv_S,
    )

#=
        if any(ocn.TFLUX_CONV .>= 0.1) 
            throw(ErrorException("TFLUX_CONV > .1"))
        end



    println("Before Ts: ", ocn.Ts[1:5, 48, 89])
    println("Before Ss: ", ocn.Ss[1:5, 48, 89])
=#

    println("ADD")
    @time @loop_hor ocn i j let
 
        Nz = ocn.Nz[i, j]
        zs   = ocn.cols.zs[i, j]
        hs   = ocn.cols.hs[i, j]
        h_ML = ocn.h_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        
        for k = 1:ocn.Nz[i, j]
#=
            if any(ocn.Ts[k, i, j] >= 60.0) 
                println("Weird Temp: ", ocn.Ts[k, i, j], "; k: ", k , ", i:  ", i, ", j:", j)
                throw(ErrorException("TFLUX_CONV > .1"))
            end

            if any(ocn.SFLUX_CONV[k, i, j] >= 60.0/1e5) 
                println("Weird Sali: ", ocn.Ss[k, i, j], "; k: ", k , ", i:  ", i, ", j:", j)
                println("SFLUX_CONV = ", ocn.SFLUX_CONV[k, i, j], "; k: ", k , ", i:  ", i, ", j:", j)
                println("SFLUX_x = ", ocn.SFLUX_DEN_x[k, i:i+1, j])
                println("SFLUX_y = ", ocn.SFLUX_DEN_y[k, i, j:j+1])
                println("SFLUX_z = ", ocn.SFLUX_DEN_z[k:k+1, i, j])
                println("SFLUX_z = ", ocn.SFLUX_DEN_z[:, i, j])
                println("w_bnd   = ", ocn.w_bnd[:, i, j])
                println("Ss_vert = ", ocn.Ss[:, i, j])
                throw(ErrorException("SFLUX_CONV > .1"))
            end
=#
            ocn.Ts[k, i, j] += Δt * ocn.TFLUX_CONV[k, i, j]
            ocn.Ss[k, i, j] += Δt * ocn.SFLUX_CONV[k, i, j]




        end


        # Adjust ΔT, ΔS
        ocn.ΔT[i, j] += ocn.dΔTdt[i, j] * Δt 
        ocn.ΔS[i, j] += ocn.dΔSdt[i, j] * Δt 
#=
        if (i, j) == (48, 89)
            println("dΔTdt = ", ocn.dΔTdt[i, j])
            println("FLDO  = ", ocn.FLDO[i, j])
            println("TFLUX_CONV: ", ocn.TFLUX_CONV[1:6, i, j])
            println("TFLUX_CONV_h: ", ocn.TFLUX_CONV_h[1:6, i, j])
            println("SFLUX_CONV: ", ocn.SFLUX_CONV[1:6, i, j])
            println("SFLUX_CONV_h: ", ocn.SFLUX_CONV_h[1:6, i, j])
            println("Before unmix Ts: ", ocn.Ts[1:5, i, j])
            println("Before unmix After Ss: ", ocn.Ss[1:5, i, j])
            println("Before unmix After bs: ", ocn.bs[1:5, i, j])
            println("u: ", ocn.u_bnd[3, i:i+1, j])
            println("v: ", ocn.v_bnd[3, i, j:j+1])
        end

=#
        ocn.T_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ts[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔT[i, j],
        )

        ocn.S_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ss[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔS[i, j],
        )

#=
        if (i, j) == (48, 89)
            println("dΔTdt = ", ocn.dΔTdt[i, j])
            println("FLDO  = ", ocn.FLDO[i, j])
            println("After unmix Ts: ", ocn.Ts[1:5, i, j])
            println("After unmix Ss: ", ocn.Ss[1:5, i, j])
            println("After unmix bs: ", ocn.bs[1:5, i, j])
        end
=#

        OC_updateB!(ocn, i, j)

        if do_convadjust
            OC_doConvectiveAdjustment!(ocn, i, j)
        end


    end
    println("DONE")

end
