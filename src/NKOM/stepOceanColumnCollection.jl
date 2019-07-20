using Statistics

function stepOceanColumnCollection!(
    occ           :: OceanColumnCollection;
    use_qflx      :: Bool,
    use_h_ML      :: Bool,
    rad_scheme    :: Symbol,    # whether to absorb radiation totally at the surface layer or not
    Δt            :: Float64,
    diffusion_Δt  :: Float64,
    relaxation_Δt :: Float64,
    do_diffusion  :: Bool = true, 
    do_relaxation :: Bool = true, 
    do_convadjust :: Bool = true,

)

    fric_u  = occ.in_flds.fric_u
    ifrac   = occ.in_flds.ifrac
    weighted_fric_u = occ.in_flds.weighted_fric_u
    taux    = occ.in_flds.taux
    tauy    = occ.in_flds.tauy
    swflx   = occ.in_flds.swflx
    nswflx  = occ.in_flds.nswflx
    frwflx  = occ.in_flds.frwflx
    qflx    = occ.in_flds.qflx
    h_ML    = occ.in_flds.h_ML

    #x = nswflx + swflx
    #println("Mean of energy flx: ", mean(x[isfinite.(x)]))
    #println("use_h_ML: ", use_h_ML)

    # It is assumed here that buoyancy has already been updated.
    for grid_idx in 1:size(occ.valid_idx)[2]


        i = occ.valid_idx[1, grid_idx]
        j = occ.valid_idx[2, grid_idx]

        zs = occ.zs_vw[i, j]
        Nz = occ.Nz[i, j]

        #fric_u[i, j] = √( √(taux[i, j]^2.0 + tauy[i, j]^2.0) / NKOM.ρ) * 0.0
        fric_u[i, j] = sqrt( sqrt(taux[i, j]^2.0 + tauy[i, j]^2.0) / NKOM.ρ )
        weighted_fric_u[i, j] = fric_u[i, j] * (1.0 - ifrac[i, j])


        # Pseudo code
        # Current using only Euler forward scheme:
        # 1. Determine h at t+Δt
        # 2. Determine how many layers are going to be
        #    taken away by ML.
        # 3. Cal b at t+Δt for both ML and DO
        # 4. Detect if it is buoyantly stable.
        #    Correct it (i.e. convection) if it is not.
        # 5. If convection happens, redetermine h.

        # p.s.: Need to examine carefully about the
        #       conservation of buoyancy in water column


        surf_Tnswflx = ( nswflx[i, j] + ( ( use_qflx ) ? qflx[i, j] : 0.0 )) / (ρ*c_p) 
        surf_Tswflx  = swflx[i, j] / (ρ*c_p)
        surf_Jflx    = g*α*surf_Tswflx
        surf_Sflx    = - frwflx[i, j] * S_surf_avg
        surf_bflx    = g * ( α * surf_Tnswflx - β * surf_Sflx )
        
        old_FLDO = occ.FLDO[i, j]
        old_h_ML = occ.h_ML[i, j]
        Δb = (old_FLDO != -1) ? occ.b_ML[i, j] - occ.bs[old_FLDO, i, j] : 0.0


        # After convective adjustment, there still might
        # be some numerical error making Δb slightly negative
        # (the one I got is like -1e-15 ~ -1e-8). So I set a
        # tolarence δb = 3e-6 ( 0.001 K => 3e-6 m/s^2 ).
        if Δb < 0.0 && -Δb <= 3e-6
            Δb = 0.0
        end

        new_h_ML = old_h_ML

        if use_h_ML # h_ML is datastream

            new_h_ML = h_ML[i, j]

        else        # h_ML is prognostic
            
            new_h_ML = calNewMLD(;
                h_ML   = old_h_ML,
                Bf     = surf_bflx + surf_Jflx * occ.R,
                J0     = surf_Jflx * (1.0 - occ.R),
                fric_u = weighted_fric_u[i, j],
                Δb     = Δb,
                f      = occ.fs[i, j],
                Δt     = Δt,
                ζ      = occ.ζ,
            )

        end

        new_h_ML = boundMLD(new_h_ML; h_ML_max=occ.h_ML_max[i, j], h_ML_min=occ.h_ML_min[i, j])

        # ML
        #      i: Calculate integrated buoyancy that should
        #         be conserved purely through entrainment
        #     ii: Add to total buoyancy

        # If new_h_ML < old_h_ML, then the FLDO layer should get extra T or S due to mixing

        if new_h_ML < old_h_ML

            new_FLDO = getFLDO(zs=zs, h_ML=new_h_ML, Nz=Nz)

            if old_FLDO == -1

                occ.Ts[new_FLDO:Nz, i, j] .= occ.T_ML[i, j]
                occ.Ss[new_FLDO:Nz, i, j] .= occ.S_ML[i, j]

            else
                FLDO_Δz =  -zs[old_FLDO+1] - old_h_ML
                retreat_Δz =  old_h_ML - ( (new_FLDO == old_FLDO) ? new_h_ML : (-zs[old_FLDO]) )

                occ.Ts[new_FLDO, i, j] = (
                    occ.Ts[old_FLDO, i, j] * FLDO_Δz + occ.T_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)

                occ.Ss[new_FLDO, i, j] = (
                    occ.Ss[old_FLDO, i, j] * FLDO_Δz + occ.S_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)
            end
        end

        # Shortwave radiation is not included yet
        new_S_ML = (OC_getIntegratedSalinity(   occ, i, j; target_z = -new_h_ML) - surf_Sflx * Δt) / new_h_ML
        new_T_ML = (OC_getIntegratedTemperature(occ, i, j; target_z = -new_h_ML) - surf_Tnswflx * Δt) / new_h_ML

        OC_setMixedLayer!(
            occ, i, j;
            T_ML=new_T_ML,
            S_ML=new_S_ML,
            h_ML=new_h_ML,
        )

        # Shortwave radiation
        if rad_scheme == :exponential
            occ.T_ML[i, j] += - occ.R * surf_Tswflx * Δt / new_h_ML
            OC_doShortwaveRadiation!(occ, i, j; Tswflx=surf_Tswflx, Δt=Δt) 
        elseif rad_scheme == :step
            occ.T_ML[i, j] += - surf_Tswflx * Δt / new_h_ML
        end

 
       
        # Climatology relaxation
        if do_relaxation
            OC_doNewtonianRelaxation_T!(occ, i, j; Δt=relaxation_Δt)
            OC_doNewtonianRelaxation_S!(occ, i, j; Δt=relaxation_Δt)
        end

        if do_diffusion
            OC_doDiffusion_EulerBackward!(occ, i, j; Δt=diffusion_Δt)
        end

        OC_updateB!(occ, i, j)

        # TODO: convective adjustment cannot break h_ML_max
        if do_convadjust
            OC_doConvectiveAdjustment!(occ, i, j;)
        end


        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - occ.T_ML[i, j]) * ρ * c_p * occ.h_ML[i, j] / Δt
        occ.T_ML[i, j] = max(T_sw_frz, occ.T_ML[i, j])
        
    end

    updateB!(occ)

    return 0
end


