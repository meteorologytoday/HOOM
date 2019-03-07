"""
    stepOceanColumnCollection!(;
        occ  :: OceanColumnCollection,
        F    :: Array{Float64},
        idx  :: Integer,
        Δt   :: Float64 
    )

# Description
This function update the OceanColumnCollection forward in time.

"""
function stepOceanColumnCollection!(;
    occ :: OceanColumnCollection,
    F   :: Array{Float64},
    idx :: Integer,
    Δt  :: Float64, 
)

    for i = 1 : occ.N_ocs

        # need to derive we, h
        occ.b_ML[i] += Δt * (
                - occ.we[i, idx] * (occ.b_ML[i] - occ.b_DO[i])
                + αgρc * (F[i] + occ.Q_ML[i, idx])
            ) / occ.h_ML[i, idx]

        T = occ.b_ML[i] / (α * g) + T_ref
        occ.qflx2atm[i] = (T_sw_frz - T) * ρ * c_p * occ.h_ML[i, idx] / Δt
        occ.b_ML[i] = max(b_sw_frz, occ.b_ML[i])
 
        #println(occ.we[i, idx], "; ", F[i] + occ.Q_ML[i, idx])
    end

end


