
function calImplied∂TEMP∂t!(
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

function calImplied∂SALT∂t!(
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


function calDirect∂TEMP∂t!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let
        
        tmp_TEMP = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1])
        ocn.TEMP[i, j], ocn.dTEMPdt[i, j] = tmp_TEMP, (tmp_TEMP - ocn.TEMP[i, j]) / Δt

    end


    #=
    tmp_budget = 0.0
    tmp_σ = 0.0
    @loop_hor ocn i j let
        #tmp_budget += (ocn.dTEMPdt[i, j] - ocn.wT[i, j] * ρc ) * ocn.mi.area[i, j]
        #tmp_σ      += ocn.mi.area[i, j]

        tmp_budget += (ocn.dTEMPdt[i, j] - ocn.wT[i, j] * ρc ) * ocn.gi.dσ[i, j]
        tmp_σ      += ocn.gi.dσ[i, j]

    end

    println("total balance check inside NKOM (W / m^2): ", tmp_budget / tmp_σ)
    
    =#

end

function calTEMP!(ocn::Ocean)

    @loop_hor ocn i j let
        ocn.TEMP[i, j] = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1])
    end
end

function calDirect∂SALT∂t!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let
        
        tmp_SALT = OC_getIntegratedSalinity(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1])
        ocn.SALT[i, j], ocn.dSALTdt[i, j] = tmp_SALT, (tmp_SALT - ocn.SALT[i, j]) / Δt

    end
    #=
    tmp_budget = 0.0
    tmp_σ = 0.0
    @loop_hor ocn i j let

        tmp_budget += (ocn.dSALTdt[i, j] - ocn.wS[i, j] ) * ocn.gi.dσ[i, j]
        tmp_σ      += ocn.gi.dσ[i, j]

    end

    println("total SALT balance check inside NKOM (W / m^2): ", tmp_budget / tmp_σ)
    
    =#

end

function calSALT!(ocn::Ocean)

    @loop_hor ocn i j let
        
        ocn.SALT[i, j] = OC_getIntegratedSalinity(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1])

    end

end
