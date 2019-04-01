function stepOceanColumnCollection!(
    occ  :: OceanColumnCollection;
    eflx :: AbstractArray{Float64, 2}, # Total downward energy flux ( J  / s m^2) downward => (-)
    Δt   :: Float64, 
)

    @sync @distributed  for idx in CartesianIndices((1:occ.Nx, 1:occ.Ny))

        i = idx[1]
        j = idx[2]

        if occ.mask[i, j] == 0.0
            continue
        end

        occ.T_ML[i, j] += -eflx[i, j] / ( occ.h_ML * ρ * c_p ) * Δt


        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - occ.T_ML[i, j]) * ρ * c_p * occ.h_ML / Δt
        occ.T_ML[i, j] = max(T_sw_frz, occ.T_ML[i, j])
        
    end
end


