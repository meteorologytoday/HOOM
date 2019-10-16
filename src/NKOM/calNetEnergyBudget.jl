function calNetEnergyBudget!(
    ocn::Ocean;
    cfgs...
)

    swflx   = ocn.in_flds.swflx
    nswflx  = ocn.in_flds.nswflx
    frwflx  = ocn.in_flds.frwflx
    qflx    = ocn.in_flds.qflx

    neb = ocn.neb
    qflx2atm = ocn.qflx2atm
    wT = ocn.wT
    Q_clim = ocn.Q_clim

    @loop_hor ocn i j let
        neb[i, j] = Q_clim[i, j] - ( nswflx[i, j] + swflx[i, j] ) - wT[i, j] * ρc + max(qflx2atm[i, j], 0.0)
    end

    if cfgs[:qflx_scheme] == :energy_flux
        @loop_hor ocn i j let
            neb[i, j] -= qflx[i, j]
        end
    end

end
