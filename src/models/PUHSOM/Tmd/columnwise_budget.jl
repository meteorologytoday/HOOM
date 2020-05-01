
function calImplied_dintTdt!(
    m :: TmdModel;
)

    @fast_extract m

    swflx    = fo.swflx
    nswflx   = fo.nswflx
    qflx_T   = fo.qflx_T

    TFLUX_DIV_implied  = dg.TFLUX_DIV_implied
    qflx2atm           = fo.qflx2atm
    TSAS_wr            = fo.TSAS_wr
    TFLUX_bot          = view(c.XFLUX_bot, :, :, 1)

    dintTdt            = dg.dintTdt

    @loop_hor m i j let
        TFLUX_DIV_implied[i, j] =  ( - ( nswflx[i, j] + swflx[i, j] ) + max(qflx2atm[i, j], 0.0) ) / ρc_sw + TSAS_wr[i, j] + TFLUX_bot[i, j] - dintTdt[i, j]
    end

    if ev.use_Qflux
        @loop_hor m i j let
            TFLUX_DIV_implied[i, j] +=  qflx_T[i, j] / ρc_sw
        end
    end

end

function calImplied_dintSdt!(
    m :: TmdModel;
)
    @fast_extract m

    SFLUX_DIV_implied  = dg.SFLUX_DIV_implied
    SFLUX_top          = view(c.XFLUX_top, :, :, 2)
    SFLUX_bot          = view(c.XFLUX_bot, :, :, 2)
    SSAS_wr            = fo.SSAS_wr
    dintSdt            = dg.dintSdt
    
    qflx_S   = fo.qflx_S

    @loop_hor m i j let
        SFLUX_DIV_implied[i, j] = SSAS_wr[i, j] - SFLUX_top[i, j] + SFLUX_bot[i, j] - dintSdt[i, j]
    end

    if ev.use_Qflx
        @loop_hor m i j let
            SFLUX_DIV_implied[i, j] += qflx_S[i, j]
        end
    end


end


function calDirect_dintTdt!(
    m :: TmdModel;
    Δt::Float64
)

    @fast_extract m

    zs = co.cols.z_bnd_av    
    @loop_hor m i j let
        
        tmp_TEMP = OC_getIntegratedTemperature(m, i, j; target_z = zs[i, j][ev.Nz_av[i, j]+1])
        dg.intT[i, j], dg.dintTdt[i, j] = tmp_TEMP, (tmp_TEMP - dg.intT[i, j]) / Δt

    end

end

function cal_intT!(
    m :: TmdModel
)
    @fast_extract m

    zs = co.cols.z_bnd_av    
    @loop_hor m i j let
        dg.intT[i, j] = OC_getIntegratedTemperature(ocn, i, j; target_z = zs[i, j][ev.Nz_av[i, j]+1])
    end
end

function calDirect_dintSdt!(
    m::TmdModel;
    Δt::Float64
)

    @fast_extract m

    zs = co.cols.z_bnd_av
    @loop_hor m i j let
        
        tmp_SALT = OC_getIntegratedSalinity(m, i, j; target_z = zs[i, j][ev.Nz_av[i, j]+1])
        dg.intS[i, j], dg.dintSdt[i, j] = tmp_SALT, (tmp_SALT - dg.intS[i, j]) / Δt

    end
end

function cal_intS!(
    m :: TmdModel
)

    @fast_extract m
    zs = co.cols.z_bnd_av
    @loop_hor m i j let
        dg.intS[i, j] = OC_getIntegratedSalinity(m, i, j; target_z = zs[i, j][ev.Nz_av[i, j]+1])
    end

end
