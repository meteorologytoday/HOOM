function calNetTEMPBudget!(
    ocn::Ocean;
    cfgs...
)

    swflx    = ocn.in_flds.swflx
    nswflx   = ocn.in_flds.nswflx
    frwflx   = ocn.in_flds.frwflx
    qflx     = ocn.in_flds.qflx

    TFLUX_DIV_implied  = ocn.TFLUX_DIV_implied
    qflx2atm           = ocn.qflx2atm
    wT_bot             = ocn.wT_bot
    TSAS_clim          = ocn.TSAS_clim
    dTEMPdt            = ocn.dTEMPdt

    @loop_hor ocn i j let
        #TFLUX_DIV_implied[i, j] = - TSAS_clim[i, j] - ( nswflx[i, j] + swflx[i, j] ) + wT_bot[i, j] * ρc + max(qflx2atm[i, j], 0.0) - dTEMPdt[i, j]
        TFLUX_DIV_implied[i, j] =  ( - ( nswflx[i, j] + swflx[i, j] ) + max(qflx2atm[i, j], 0.0) ) / ρc + TSAS_clim[i, j] + wT_bot[i, j] - dTEMPdt[i, j]
    end

    if cfgs[:qflx_scheme] == :energy_flux
        @loop_hor ocn i j let
            TFLUX_DIV_implied[i, j] +=  - qflx[i, j] / ρc
        end
    end

end

function calNetSALTBudget!(
    ocn::Ocean;
    cfgs...
)

    SFLUX_DIV_implied  = ocn.SFLUX_DIV_implied
    wS_top             = ocn.wS_top
    wS_bot             = ocn.wS_bot
    SSAS_clim          = ocn.SSAS_clim
    dSALTdt            = ocn.dSALTdt

    @loop_hor ocn i j let
        SFLUX_DIV_implied[i, j] = SSAS_clim[i, j] - wS_top[i, j] + wS_bot[i, j] - dSALTdt[i, j]
    end

end
