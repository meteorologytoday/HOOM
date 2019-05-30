function stepOceanColumnCollection!(
    occ    :: OceanColumnCollection;
    fric_x :: AbstractArray{Float64, 2},
    fric_y :: AbstractArray{Float64, 2},
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


    for



    #@time @sync @distributed  for idx in CartesianIndices((1:occ.Nx, 1:occ.Ny))
    for idx in CartesianIndices((1:occ.Nx, 1:occ.Ny))

        i = idx[1]
        j = idx[2]

        if occ.mask[i, j] == 0
            continue
        end


        # Surface fluxes

        total_Tflx = ( swflx[i, j] + nswflx[i, j] ) / (ρ*c_p) 
        total_Sflx = - frwflx[i, j] * S_surf_avg
        total_bflx = g * ( α * total_Tflx - β * total_Sflx )
       
        occ.Ts[i, j, 1] -= total_Tflx * Δt / occ.hs[1]
        occ.Ss[i, j, 1] -= total_Sflx * Δt / occ.hs[1]


        #  


        OC_setMixedLayer!(
            occ, i, j;
            T_ML=new_T_ML,
            S_ML=new_S_ML,
            h_ML=new_h_ML,
        )
 
       
        OC_doDiffusion_EulerBackward!(occ, i, j; Δt=Δt)

        OC_updateB!(occ, i, j)

        # TODO: convective adjustment cannot break h_ML_max
        OC_doConvectiveAdjustment!(occ, i, j;)


        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - occ.T_ML[i, j]) * ρ * c_p * occ.h_ML[i, j] / Δt
        occ.T_ML[i, j] = max(T_sw_frz, occ.T_ML[i, j])
        
    end

    updateB!(occ)
end


