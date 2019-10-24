function calNetEnergyBudget!(
    ocn::Ocean;
    cfgs...
)

    swflx    = ocn.in_flds.swflx
    nswflx   = ocn.in_flds.nswflx
    frwflx   = ocn.in_flds.frwflx
    qflx     = ocn.in_flds.qflx

    TFLUX_DIV_implied  = ocn.TFLUX_DIV_implied
    qflx2atm           = ocn.qflx2atm
    wT                 = ocn.wT
    Q_clim             = ocn.Q_clim
    dHdt               = ocn.dHdt

    @loop_hor ocn i j let
        TFLUX_DIV_implied[i, j] = - Q_clim[i, j] - ( nswflx[i, j] + swflx[i, j] ) + wT[i, j] * ρc + max(qflx2atm[i, j], 0.0) - dHdt[i, j]
    end

    if cfgs[:qflx_scheme] == :energy_flux
        @loop_hor ocn i j let
            TFLUX_DIV_implied[i, j] -= qflx[i, j]
        end
    end

end
