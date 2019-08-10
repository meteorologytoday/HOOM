function stepOcean_hz!(
    ocn  :: Ocean;
    cfgs...
)

    # Pseudo code
    # 1. assign velocity field
    # 2. calculate temperature & salinity flux
    # 3. calculate temperature & salinity flux divergence
    # 4. 


    #=
    for grid_idx in 1:size(ocn.valid_idx)[2]


        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]


        ocn.T_ML[i, j] += 0.5


    end
    =#
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
    calEkmanTransport!(occ, wksp.τx, wksp.τy, wksp.M1x, wksp.M1y, occ.fs, occ.ϵs)
    
    # Calculate ∇⋅M1
    DisplacedPoleCoordinate.divergence2!(occ.gi, wksp.M1x, wksp.M1y, wksp.div_M1; mask=occ.mask)
    
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
