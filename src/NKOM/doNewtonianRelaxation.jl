function OC_doNewtonianRelaxation_T!(
    occ :: OceanColumnCollection,
    i   :: Integer,
    j   :: Integer;
    Δt  :: Float64,
)

    doNewtonianRelaxation!(
        qs      = occ.Ts_vw[i, j],
        qs_clim = occ.Ts_clim_vw[i, j],
        FLDO    = occ.FLDO[i, j],
        Nz      = occ.Nz[i, j],
        τ       = occ.Ts_clim_relax_time,
        Δt      = Δt,
    )

end

function OC_doNewtonianRelaxation_S!(
    occ :: OceanColumnCollection,
    i   :: Integer,
    j   :: Integer;
    Δt  :: Float64,
)

    doNewtonianRelaxation!(
        qs      = occ.Ss_vw[i, j],
        qs_clim = occ.Ss_clim_vw[i, j],
        FLDO    = occ.FLDO[i, j],
        Nz      = occ.Nz[i, j],
        τ       = occ.Ss_clim_relax_time,
        Δt      = Δt,
    )

end

"""

    This function newtonian-relaxes `qs` to `qs_clim` with e-folding time `τ`.
    It has to be noted that mixed-layer is not relaxed.

    Also, Euler backward integration scheme is used.

"""
function doNewtonianRelaxation!(;
    qs         :: AbstractArray{Float64, 1},
    qs_clim    :: AbstractArray{Float64, 1},
    FLDO       :: Integer,
    Nz         :: Integer,
    τ          :: Float64,
    Δt         :: Float64,
)

    if τ > 0.0

        r = Δt / τ
        if FLDO != -1
            for i = FLDO:Nz
                qs[i] = (qs[i] + r * qs_clim[i]) / (1+r)
            end
        end

    elseif τ = 0.0
 
        if FLDO != -1
            for i = FLDO:Nz
                qs[i] = qs_clim[i]
            end
        end
   
    end
end

