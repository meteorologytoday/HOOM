function calNetTEMPBudget!(
    ocn::Ocean;
    cfgs...
)

    swflx    = ocn.in_flds.swflx
    nswflx   = ocn.in_flds.nswflx
    frwflx   = ocn.in_flds.frwflx
    qflx_T   = ocn.in_flds.qflx_T

    TFLUX_DIV_implied  = ocn.TFLUX_DIV_implied
    qflx2atm           = ocn.qflx2atm
    TFLUX_bot          = ocn.TFLUX_bot
    TSAS_clim          = ocn.TSAS_clim
    dTEMPdt            = ocn.dTEMPdt

    @loop_hor ocn i j let
        #TFLUX_DIV_implied[i, j] = - TSAS_clim[i, j] - ( nswflx[i, j] + swflx[i, j] ) + TFLUX_bot[i, j] * ρc + max(qflx2atm[i, j], 0.0) - dTEMPdt[i, j]
        TFLUX_DIV_implied[i, j] =  ( - ( nswflx[i, j] + swflx[i, j] ) + max(qflx2atm[i, j], 0.0) ) / ρc + TSAS_clim[i, j] + TFLUX_bot[i, j] - dTEMPdt[i, j]
    end

    if cfgs[:do_qflx]
        @loop_hor ocn i j let
            TFLUX_DIV_implied[i, j] +=  - qflx_T[i, j] / ρc
        end
    end

end

function calNetSALTBudget!(
    ocn::Ocean;
    cfgs...
)

    SFLUX_DIV_implied  = ocn.SFLUX_DIV_implied
    SFLUX_top             = ocn.SFLUX_top
    SFLUX_bot             = ocn.SFLUX_bot
    SSAS_clim          = ocn.SSAS_clim
    dSALTdt            = ocn.dSALTdt
    
    qflx_S   = ocn.in_flds.qflx_S

    @loop_hor ocn i j let
        SFLUX_DIV_implied[i, j] = SSAS_clim[i, j] - SFLUX_top[i, j] + SFLUX_bot[i, j] - dSALTdt[i, j]
    end

    if cfgs[:do_qflx]
        @loop_hor ocn i j let
            SFLUX_DIV_implied[i, j] +=  - qflx_S[i, j] / ρ
        end
    end


end
