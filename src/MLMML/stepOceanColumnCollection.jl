"""

    stepOceanColumn!(;

        occ    :: OceanColumnCollection;
        fric_u :: AbstractArray{Float64, 2}, # Currently assumed to be u10
        swflx  :: AbstractArray{Float64, 2}, # Shortwave     energy flux at the surface (     J  / s m^2)
        nswflx :: AbstractArray{Float64, 2}, # Non-shortwave energy flux at the surface (     J  / s m^2)
        frwflx :: AbstractArray{Float64, 2}, # Freshwater           flux at the surface (     m  / s m^2)
        Δt     :: Float64, 

    )

# Description
This function update the OceanColumnCollection forward in time.

"""
function stepOceanColumnCollection!(
    occ    :: OceanColumnCollection;
    fric_u :: AbstractArray{Float64, 2}, # Currently assumed to be u10
    swflx  :: AbstractArray{Float64, 2}, # Shortwave     energy flux at the surface (     J  / s m^2)
    nswflx :: AbstractArray{Float64, 2}, # Non-shortwave energy flux at the surface (     J  / s m^2)
    frwflx :: AbstractArray{Float64, 2}, # Freshwater           flux at the surface (     m  / s m^2)
    Δt     :: Float64, 
)

    # It is assumed here that buoyancy has already been updated.
    #@time @sync @distributed  for idx in CartesianIndices((1:occ.Nx, 1:occ.Ny))
    for idx in CartesianIndices((1:occ.Nx, 1:occ.Ny))

        i = idx[1]
        j = idx[2]

        if occ.mask[i, j] == 0.0
            continue
        end

        zs = occ.zs_vw[i, j]
        Nz = occ.Nz[i, j]

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

        #println("### h: ", oc.h)
        #println("FLDO:", oc.FLDO)

        total_Tflx = ( swflx[i, j] + nswflx[i, j] ) / (ρ*c_p) 
        total_Sflx = - frwflx[i, j] * S_surf_avg
        total_bflx = g * ( α * total_Tflx - β * total_Sflx )
        
        old_FLDO = occ.FLDO[i, j]
        old_h_ML = occ.h_ML[i, j]
        Δb = (old_FLDO != -1) ? occ.b_ML[i, j] - occ.bs[i, j, old_FLDO] : 0.0

        # After convective adjustment, there still might
        # be some numerical error making Δb slightly negative
        # (the one I got is like -1e-15). So I set a tolarence
        # ( 0.001 K ≈ 3e-6 m/s^2 ).
        if Δb < 0.0 && -Δb <= 3e-6
            Δb = 0.0
        end



        #fric_u = getFricU(ua=ua)
        flag, val = calWeOrMLD(;
            h_ML   = old_h_ML,
            B      = total_bflx,
            fric_u = fric_u[i, j],
            Δb     = Δb
        )
        #println("Before:" , oc.bs[10], "; oc.FLDO = ", oc.FLDO, "; Δb = ", Δb)
        #println("B: ", total_flx , ";Δb: ", Δb , "; fric_u: ", fric_u[i, j])

        # 1
        if flag == :MLD
            we = 0.0
            new_h_ML  = val
        elseif flag == :we
            we = val 
            new_h_ML = old_h_ML + Δt * we
        end

        #println("h_ML_min: ", occ.h_ML_min, "; h_ML_max: ", occ.h_ML_max)    
        new_h_ML = boundMLD(new_h_ML; h_ML_max=occ.h_ML_max[i, j], h_ML_min=occ.h_ML_min[i, j])

        #println("flag: ", String(flag), "; val: ", val, "; new_h_ML: ", new_h_ML)
        # 2
        # 3

        # ML
        #      i: Calculate integrated buoyancy that should
        #         be conserved purely through entrainment
        #     ii: Add to total buoyancy


        # If new_h_ML < old_h_ML, then the FLDO layer should get extra T or S due to mixing

        if new_h_ML < old_h_ML

            new_FLDO = getFLDO(zs=zs, h_ML=new_h_ML, Nz=Nz)

            if old_FLDO == -1

                occ.Ts[i, j, new_FLDO:Nz] .= occ.T_ML[i, j]
                occ.Ss[i, j, new_FLDO:Nz] .= occ.S_ML[i, j]

            else
                FLDO_Δz =  -zs[old_FLDO+1] - old_h_ML
                retreat_Δz =  old_h_ML - ( (new_FLDO == old_FLDO) ? new_h_ML : (-zs[old_FLDO]) )

                occ.Ts[i, j, new_FLDO] = (
                    occ.Ts[i, j, old_FLDO] * FLDO_Δz + occ.T_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)

                occ.Ss[i, j, new_FLDO] = (
                    occ.Ss[i, j, old_FLDO] * FLDO_Δz + occ.S_ML[i, j] * retreat_Δz
                ) / (FLDO_Δz + retreat_Δz)
            end
        end

        #println("target_z: ", -new_h_ML)
        
        new_T_ML = (OC_getIntegratedTemperature(occ, i, j; target_z = -new_h_ML) - total_Tflx * Δt) / new_h_ML
        new_S_ML = (OC_getIntegratedSalinity(   occ, i, j; target_z = -new_h_ML) - total_Sflx * Δt) / new_h_ML


        OC_setMixedLayer!(
            occ, i, j;
            T_ML=new_T_ML,
            S_ML=new_S_ML,
            h_ML=new_h_ML,
        )
 
       
        # Climatology relaxation
        if occ.Ts_clim != nothing
            OC_doNewtonianRelaxation!(occ, i, j; Δt=Δt)
        end

        if occ.Ss_clim != nothing
            OC_doNewtonianRelaxation!(occ, i, j; Δt=Δt)
        end


        OC_doDiffusion_EulerBackward!(occ, i, j; Δt=Δt)

        OC_updateB!(occ, i, j)
        OC_doConvectiveAdjustment!(occ, i, j;)


        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - occ.T_ML[i, j]) * ρ * c_p * occ.h_ML[i, j] / Δt
        occ.T_ML[i, j] = max(T_sw_frz, occ.T_ML[i, j])
        
    end

    updateB!(occ)
end


