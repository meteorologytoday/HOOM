function stepOceanColumnCollection!(
    occ    :: OceanColumnCollection;
    τ_x :: AbstractArray{Float64, 2},
    τ_y :: AbstractArray{Float64, 2},
    swflx  :: AbstractArray{Float64, 2}, # Shortwave     energy flux at the surface (     J  / s m^2)
    nswflx :: AbstractArray{Float64, 2}, # Non-shortwave energy flux at the surface (     J  / s m^2)
    frwflx :: AbstractArray{Float64, 2}, # Freshwater           flux at the surface (     m  / s m^2)
    Δt     :: Float64, 
)

    # It is assumed here that buoyancy has already been updated.

    # Pseudo code
    # 1. Do heating
    # 2. Do Ekman transport
    # 3. Do horizontal diffusion
    # 4. Do convective adjustment


    # Gov Eqn: ∂T1/∂t = - Fsurf / (ρ0 cp H1) - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )
    # Gov Eqn: ∂T2/∂t =                      - 1 / (ρ H2) ( ∇⋅(M2 T2) + (∇⋅M1) Tmid )

    wksp = occ.wksp
    
    DisplacedPoleCoordinate.project!(occ.gi, τx, τy, wksp.τx, wksp.τy, direction=:Forward)
    calEkmanTransport!(occ, wksp.τx, wksp.τy, wksp.Mx, wksp.My, occ.fs, occ.ϵs)
    
    # Calculate ∇⋅M1
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.Mx, wksp.My, wksp.div_M1)
    
    # Calculate ( M1 T1 ), ( M2 T2 ) 
    for i = 1:occ.Nx, j = 1:occ.Ny

        if occ.mask[i, j] == 0
            continue
        end

        Mx = wksp.Mx[i, j]
        My = wksp.My[i, j]

        wksp.Mx_T1[i, j] =   Mx * occ.Ts[i, j, 1]
        wksp.My_T1[i, j] =   My * occ.Ts[i, j, 1]
        wksp.Mx_T2[i, j] = - Mx * occ.Ts[i, j, 2]
        wksp.My_T2[i, j] = - My * occ.Ts[i, j, 2]

        wksp.Mx_S1[i, j] =   Mx * occ.Ss[i, j, 1]
        wksp.My_S1[i, j] =   My * occ.Ss[i, j, 1]
        wksp.Mx_S2[i, j] = - Mx * occ.Ss[i, j, 2]
        wksp.My_S2[i, j] = - My * occ.Ss[i, j, 2]

    end 

    # Calculate ∇⋅(M1 T1), ∇⋅(M2 T2), ∇⋅(M1 S1), ∇⋅(M2 S2)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.Mx_T1, wksp.My_T1, wksp.div_MT1)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.Mx_T2, wksp.My_T2, wksp.div_MT2)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.Mx_S1, wksp.My_S1, wksp.div_MS1)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.Mx_S2, wksp.My_S2, wksp.div_MS2)

    # Step forward temperature
    for i = 1:occ.Nx, j = 1:occ.Ny

        if occ.mask[i, j] == 0
            continue
        end

        total_Tflx = ( swflx[i, j] + nswflx[i, j] ) / (ρ*c_p) 

        Tmid = occ.Ts[i, j, ( wksp.div_M1 < 0 ) ? 1 : 2 ]
        occ.Ts[i, j, 1] += ( - (swflx[i, j] + nswflx[i, j]) / c_p - (wksp.div_MT1[i, j] - wksp.div_M1[i, j] * Tmid)) / (ρ * occ.hs[1]) * Δt
        occ.Ts[i, j, 2] +=                                        - (wksp.div_MT2[i, j] + wksp.div_M1[i, j] * Tmid)  / (ρ * occ.hs[2]) * Δt

        Smid = occ.Ss[i, j, ( wksp.div_M1 < 0 ) ? 1 : 2 ]
        occ.Ss[i, j, 1] += ( frwflx[i, j] * S_surf_avg * ρ - (wksp.div_MS1[i, j] - wksp.div_M1[i, j] * Smid)) / (ρ * occ.hs[1]) * Δt
        occ.Ss[i, j, 2] +=                                 - (wksp.div_MS2[i, j] + wksp.div_M1[i, j] * Smid)  / (ρ * occ.hs[2]) * Δt


        OC_updateB!(occ, i, j)
        OC_doConectiveAdjustment!(occ, i, j)
        
        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - occ.Ts[i, j, 1]) * ρ * c_p * occ.hs[1] / Δt
        occ.Ts[i, j, 1] = max(T_sw_frz, occ.Ts[i, j, 1])

    end
 
end


