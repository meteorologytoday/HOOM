function OC_doShortwaveRadiation!(
    occ    :: OceanColumnCollection,
    i      :: Integer,
    j      :: Integer;
    Tswflx :: Float64,
    Δt     :: Float64,
)
    doShortwaveRadiation!(
        Tswflx = Tswflx,
        Ts   = occ.Ts_vw[i, j],
        zs   = occ.zs_vw[i, j],
        hs   = occ.hs_vw[i, j],
        rad_decay_coes  = occ.rad_decay_coes_vw[i, j],
        rad_absorp_coes = occ.rad_absorp_coes_vw[i, j],
        T_ML = occ.T_ML[i, j],
        h_ML = occ.h_ML[i, j],
        Nz   = occ.Nz[i, j],
        FLDO = occ.FLDO[i, j],
        ζ    = occ.ζ,
        Δt   = Δt,

    )
end


function doShortwaveRadiation!(;
    Tswflx          :: Float64,
    Ts              :: AbstractArray{Float64, 1}, 
    zs              :: AbstractArray{Float64, 1}, 
    hs              :: AbstractArray{Float64, 1}, 
    rad_decay_coes  :: AbstractArray{Float64, 1}, 
    rad_absorp_coes :: AbstractArray{Float64, 1}, 
    T_ML            :: Float64,
    h_ML            :: Float64,
    Nz              :: Integer,
    FLDO            :: Integer,
    ζ               :: Float64,
    Δt              :: Float64,
)

    
    # ===== [BEGIN] Mixed layer =====

    if FLDO == -1      # Entire ocean column is mixed-layer
        T_ML += - Tswflx * Δt / h_ML
        return
    end


    rad_decay_coes_ML = exp(-γ*h_ML)
    T_ML += - Tswflx * (1.0 - rad_decay_coes_ML) * Δt / h_ML
    
    # ===== [END] Mixed layer =====

    # ===== [BEGIN] FLDO layer =====

    h_FLDO = -h_ML - zs[FLDO+1]

    if FLDO == Nz  # FLDO is last layer
        Ts[FLDO] += - Tswflx * rad_decay_coes_ML * Δt / h_FLDO
        return
    else
        Ts[FLDO] += - Tswflx * (rad_decay_coes_ML - rad_decay_coes[FLDO+1]) * Δt / h_FLDO
    end
    # ===== [END] FLDO layer =====

    # ===== [BEGIN] Rest layers =====

    for k=FLDO+1:Nz-1
        Ts[k] += - Tswflx * rad_decay_coes[k] * rad_absorp_coes[k] * Δt / hs[k]
    end

    # ===== [END] Rest layers =====

end

