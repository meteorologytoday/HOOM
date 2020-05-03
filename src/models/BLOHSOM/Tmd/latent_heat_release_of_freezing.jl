function cleanQflx2atm!(m :: TmdModel)

    f = m.forcing

    @loop_hor m.core i j let

        f.qflx2atm[i, j] = 0.0 

    end

end

function calLatentHeatReleaseOfFreezing!(
    m :: TmdModel;
    Δt::Float64,
)
    f = m.forcing
    s = m.state
    c = m.core
    @loop_hor c i j let

        #OC_updateB!(ocn, i, j)

        #=
        old_FLDO=ocn.FLDO[i,j]
        Δb = (old_FLDO == -1 ) ? 0.0 : ocn.b_ML[i, j] - ocn.bs[old_FLDO, i, j]

        if Δb < -3e-6
            println(format("[Qflx2ATM] [BEFORE] At {:d}, {:d}: Δb = {:f}. Even after adjustment.", i, j, Δb))
        end
        =#

        ocn.qflx2atm[i, j] = (T_sw_frz - ocn.T_ML[i, j]) * ρc_sw * ocn.h_ML[i, j] / Δt
        ocn.qflx2atm_pos[i, j] = max(ocn.qflx2atm[i,j], 0.0)
        ocn.qflx2atm_neg[i, j] = min(ocn.qflx2atm[i,j], 0.0)

        if ocn.T_ML[i, j] < T_sw_frz
            ocn.T_ML[i, j] = T_sw_frz
            ocn.Ts[1:((ocn.FLDO[i, j] == -1) ? ocn.Nz[i, j] : ocn.FLDO[i, j]-1 ), i, j] .= T_sw_frz
            OC_updateB!(ocn, i, j)
        end

    end

    if m.env.convective_adjustment
        @loop_hor c i j let
            OC_doConvectiveAdjustment!(m, i, j)
        end
    end

end
