function doMixedLayerDynamics!(
    m :: TmdModel,
)

    @fast_extract m

    ifrac   = fr.ifrac

    taux    = fr.τx
    tauy    = fr.τy

    swflx   = fr.swflx
    nswflx  = fr.nswflx
    vsflx   = fr.vsflx
    qflx_T  = fr.qflx_T
    qflx_S  = fr.qflx_S

    Δt      = ev.Δt_substep        

    # It is assumed here that buoyancy has already been updated.
    @loop_hor m i j let

        OC_updateB!(m, i, j)
        
        if ev.convective_adjustment
            OC_doConvectiveAdjustment!(m, i, j;)
        end


        zs = co.cols.z_bnd_av[i, j]
        Nz = ev.Nz_av[i, j]

        fric_u          = √( √(taux[i, j]^2.0 + tauy[i, j]^2.0) / ρ_sw)
        weighted_fric_u = fric_u * (1.0 - ifrac[i, j])

#=
        if (i, j) == (35, 20)
            println("Just get in ")
            println(st.T[:, i, j]) 
            println(st.S[:, i, j]) 
            println(st.b[:, i, j]) 
        end
=#

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

        old_FLDO = st.FLDO[i, j]
        old_h_ML = st.h_ML[i, j]
        old_T_ML = st.T_ML[i, j]
        old_S_ML = st.S_ML[i, j]

        α = TS2α(old_T_ML, old_S_ML) 
        β = TS2β(old_T_ML, old_S_ML) 

        surf_Tnswflx = nswflx[i, j] / ρc_sw 
        surf_Tswflx  = swflx[i, j] / ρc_sw
        surf_Jflx    = g * α * surf_Tswflx
        #surf_Sflx    = - frwflx[i, j] * ocn.S_ML[i, j] / ρ_fw 
        surf_Sflx    = vsflx[i, j]
        surf_bflx    = g * ( α * surf_Tnswflx - β * surf_Sflx )
        
        co.XFLUX_top[i, j, 2] = surf_Sflx

        new_h_ML = old_h_ML

        if ev.MLT_scheme == :prescribe # h_ML is datastream

            new_h_ML = fr.h_ML[i, j]

        else        # h_ML is prognostic
 

            Δb = (old_FLDO == -1 ) ? 0.0 : st.b_ML[i, j] - st.b[old_FLDO, i, j]

            # After convective adjustment, there still might
            # be some numerical error making Δb slightly negative
            # (the one I got is like -1e-15 ~ -1e-8). So I set a
            # tolarence δb = 3e-6 ( 0.001 K => 3e-6 m/s^2 ).
            if -3e-6 < Δb < 0.0
                Δb = 0.0
            end

            #=
            if Δb < -3e-6
                println(format("[MLDynamics] At {:d}, {:d}: Δb = {:f}. FLDO = {:d}", i, j, Δb, old_FLDO))
                println("See if T_ML and Ts are inconsistent: ")
                println("T_ML: ", ocn.T_ML[i, j])
                println("Ts  : ", ocn.Ts[:, i, j])
                println(format("ΔT  : {:e}", ocn.T_ML[i, j] - ocn.Ts[1, i, j]))
                println("See if S_ML and Ss are inconsistent: ")
                println("S_ML: ", ocn.S_ML[i, j])
                println("Ss  : ", ocn.Ss[:, i, j])
                println(format("ΔS  : {:e}", ocn.S_ML[i, j] - ocn.Ss[1, i, j]))
            end
=#
#            if Δb < 0.0
#                FLDO = ocn.FLDO[i, j]

#=
                println(format("({:d},{:d}) Averge sense Δb={:f}", i, j, Δb))
                println(format("({:d},{:d}) Jump sense Δb={:f}", i, j, (FLDO != -1) ? ocn.b_ML[i, j] - ocn.bs[FLDO, i, j] : 999 ))
                println("##### Ts: ", ocn.Ts[:, i, j])
                println("##### Ss: ", ocn.Ss[:, i, j])
                println("##### bs: ", ocn.bs[:, i, j])
=#
#            end

            new_h_ML, st.h_MO[i, j] = calNewMLD(;
                h_ML   = old_h_ML,
                Bf     = surf_bflx + surf_Jflx * ev.R,
                J0     = surf_Jflx * (1.0 - ev.R),
                fric_u = weighted_fric_u,
                Δb     = Δb,
                f      = ev.gi.c_f[i, j],
                Δt     = Δt,
                ζ      = ev.ζ,
                h_max  = ev.h_ML_max[i, j],
                we_max = ev.we_max,
            )
            
        end

        new_h_ML = boundMLD(new_h_ML; h_ML_max=ev.h_ML_max[i, j], h_ML_min=ev.h_ML_min[i, j])

        if new_h_ML == 0.0
            throw(ErrorException("zero h_ML"))
        end

        # ML
        #      i: Calculate integrated buoyancy that should
        #         be conserved purely through entrainment
        #     ii: Add to total buoyancy

        # If new_h_ML < old_h_ML, then the FLDO layer should get extra T or S due to mixing


        if new_h_ML < old_h_ML

            new_FLDO = getFLDO(zs=zs, h_ML=new_h_ML, Nz=Nz)

            if old_FLDO == -1

                # Mixing does not happen because FLDO does not exist in this case
                st.T[new_FLDO:Nz, i, j] .= st.T_ML[i, j]
                st.S[new_FLDO:Nz, i, j] .= st.S_ML[i, j]

            else
                FLDO_Δz =  -old_h_ML - zs[old_FLDO+1]
                retreat_Δz =  old_h_ML - ( (new_FLDO == old_FLDO) ? new_h_ML : (-zs[old_FLDO]) )

                st.T[old_FLDO, i, j] = (
                    st.T[old_FLDO, i, j] * FLDO_Δz + st.T_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)

                st.S[old_FLDO, i, j] = (
                    st.S[old_FLDO, i, j] * FLDO_Δz + st.S_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)
            end
        end

        if_entrainment = new_h_ML > old_h_ML

        # Calculate the effect of entrainment on SSS
        new_int_S_ML = OC_getIntegratedSalinity(   m, i, j; target_z = -new_h_ML)
        new_S_ML = new_int_S_ML / new_h_ML
        #ocn.dSdt_ent[i, j] = (if_entrainment) ? (new_S_ML - old_S_ML) / Δt : 0.0

        # Add in external surface flux effect on SSS
        new_S_ML = (new_int_S_ML - surf_Sflx * Δt) / new_h_ML

        # Calculate the effect of entrainment on SST
        new_int_T_ML = OC_getIntegratedTemperature(  m, i, j; target_z = -new_h_ML)
        new_T_ML = new_int_T_ML / new_h_ML
        #ocn.dTdt_ent[i, j] = (if_entrainment) ? (new_T_ML - old_T_ML) / Δt : 0.0

        # Add in external surface flux effect on SST. Shortwave radiation is not included yet
        new_T_ML = (new_int_T_ML - surf_Tnswflx * Δt) / new_h_ML

        # Q-flux 
        if ev.use_Qflux

            new_T_ML += qflx_T[i, j] * Δt / (ρc_sw * new_h_ML)
            new_S_ML += qflx_S[i, j] * Δt / new_h_ML

        end

#=
        if (i, j) == (35, 20)
            println("before setMixedLayer")
            println(st.T[:, i, j]) 
            println(st.S[:, i, j]) 
            println(st.b[:, i, j]) 
        end
=#

        # Update mixed-layer
        OC_setMixedLayer!(
            m, i, j;
            T_ML=new_T_ML,
            S_ML=new_S_ML,
            h_ML=new_h_ML,
        )

        # Shortwave radiation
        if ev.radiation_scheme == :exponential_decay
            FLDO = st.FLDO[i, j]
            st.T_ML[i, j] += - ev.R * surf_Tswflx * Δt / new_h_ML
            st.T[1:((FLDO == -1) ? Nz : FLDO-1 ), i, j] .= st.T_ML[i, j]
            OC_doShortwaveRadiation!(m, i, j; Tswflx=(1.0 - ev.R) * surf_Tswflx, Δt=Δt)

        elseif ev.radiation_scheme == :step
            FLDO = st.FLDO[i, j]
            st.T_ML[i, j] += - surf_Tswflx * Δt / new_h_ML
            st.T[1:((FLDO == -1) ? Nz : FLDO-1 ), i, j] .= st.T_ML[i, j]

        end

#=
        if (i, j) == (35, 20)
            println("before updateB")
            println(st.T[:, i, j]) 
            println(st.S[:, i, j]) 
            println(st.b[:, i, j]) 
        end
=#
        OC_updateB!(m, i, j)
        
        if ev.convective_adjustment
            if OC_doConvectiveAdjustment!(m, i, j;) && (i,j)==(35,20)
                println("!!!!! ", i, "; ", j)
            end
        end

    end

end

