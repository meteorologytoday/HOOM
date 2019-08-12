function stepOcean_hz!(
    ocn  :: Ocean;
    cfgs...
)
    
    # Transform input wind stress vector first
    DisplacedPoleCoordinate.project!(ocn.gi, ocn.in_flds.τx, ocn.in_flds.τy, ocn.τx, ocn.τy, direction=:Forward)

    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        ϵ = ocn.ϵs[i, j]
        f = ocn.fs[i, j]

        τx = ocn.τx[i, j]
        τy = ocn.τy[i, j]

        h_ML = ocn.h_ML[i, j]
    
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

            ocn.u[FLDO, i, j] = ek_u * (h_ML - )
            ocn.u[FLDO+
            


    end
        
        # Pseudo code
        # 1. assign velocity field
        # 2. calculate temperature & salinity flux
        # 3. calculate temperature & salinity flux divergence
        # Gov eqn adv + diff: ∂T/∂t = - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )

    # Calculate ∇⋅v
    DisplacedPoleCoordinate.divergence2!(ocn.gi, ocn.u, ocn.v, wksp.div; mask=ocn.mask)
    
    # Calculate ∇⋅(vT)
    DisplacedPoleCoordinate.divergence2!(ocn.gi, ocn.uT, ocn.vT, wksp.divTF; mask=ocn.mask)
    # Calculate ∇⋅(vS)
    DisplacedPoleCoordinate.divergence2!(ocn.gi, ocn.uS, ocn.vS, wksp.divSF; mask=ocn.mask)
    
    # Calculate ∇²T, ∇²S
    DisplacedPoleCoordinate.cal∇²!(occ.gi, ocn.lays.Ts[k], ocn.∇²T ; mask=ocn.mask)

# DONE. Let vertical to the rest. 
    
    

    # Calculate ( M1 T1 ), ( M2 T2 ) 
    for i = 1:occ.Nx, j = 1:occ.Ny

        if occ.mask[i, j] == 0
            continue
        end

        M1x = wksp.M1x[i, j]
        M1y = wksp.M1y[i, j]

        wksp.M1x_T1[i, j] =   M1x * occ.Ts[i, j, 1]
        wksp.M1y_T1[i, j] =   M1y * occ.Ts[i, j, 1]
        wksp.M1x_T2[i, j] = - M1x * occ.Ts[i, j, 2]
        wksp.M1y_T2[i, j] = - M1y * occ.Ts[i, j, 2]

        wksp.M1x_S1[i, j] =   M1x * occ.Ss[i, j, 1]
        wksp.M1y_S1[i, j] =   M1y * occ.Ss[i, j, 1]
        wksp.M1x_S2[i, j] = - M1x * occ.Ss[i, j, 2]
        wksp.M1y_S2[i, j] = - M1y * occ.Ss[i, j, 2]

    end 

    # Calculate ∇⋅(M1 T1), ∇⋅(M2 T2), ∇⋅(M1 S1), ∇⋅(M2 S2)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.M1x_T1, wksp.M1y_T1, wksp.div_M1T1; mask=occ.mask)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.M1x_T2, wksp.M1y_T2, wksp.div_M1T2; mask=occ.mask)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.M1x_S1, wksp.M1y_S1, wksp.div_M1S1; mask=occ.mask)
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.M1x_S2, wksp.M1y_S2, wksp.div_M1S2; mask=occ.mask)
    
    # Calculate ∇²T1, ∇²T2
    DisplacedPoleCoordinate.cal∇²!(occ.gi, view(occ.Ts, :, :, 1), wksp.∇²T1; mask=occ.mask)
    DisplacedPoleCoordinate.cal∇²!(occ.gi, view(occ.Ts, :, :, 2), wksp.∇²T2; mask=occ.mask)


    # Step forward temperature
    for i = 1:occ.Nx, j = 1:occ.Ny

        if occ.mask[i, j] == 0
            continue
        end

        active_layer = ( wksp.div_M1[i, j] < 0 ) ? 1 : 2

        Tmid = occ.Ts[i, j, active_layer]
        occ.Ts[i, j, 1] += ( ( - (swflx[i, j] + nswflx[i, j]) / c_p - (wksp.div_M1T1[i, j] - wksp.div_M1[i, j] * Tmid)) / (ρ * occ.hs[1]) + occ.Kh_T * wksp.∇²T1[i, j] ) * Δt
        occ.Ts[i, j, 2] += ( - (wksp.div_M1T2[i, j] + wksp.div_M1[i, j] * Tmid)  / (ρ * occ.hs[2]) + occ.Kh_T * wksp.∇²T2[i, j] ) * Δt


end
