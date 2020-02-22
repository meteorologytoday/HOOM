function cleanQflx2atm!(ocn::Ocean)

    @loop_hor ocn i j let

        ocn.qflx2atm[i, j] = 0.0 

    end

end

function calLatentHeatReleaseOfFreezing!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let

        #OC_updateB!(ocn, i, j)

        #=
        old_FLDO=ocn.FLDO[i,j]
        Δb = (old_FLDO == -1 ) ? 0.0 : ocn.b_ML[i, j] - ocn.bs[old_FLDO, i, j]

        if Δb < -3e-6
            println(format("[Qflx2ATM] [BEFORE] At {:d}, {:d}: Δb = {:f}. Even after adjustment.", i, j, Δb))
        end
        =#

        ocn.qflx2atm[i, j] = (T_sw_frz - ocn.T_ML[i, j]) * ρc * ocn.h_ML[i, j] / Δt

        if ocn.T_ML[i, j] < T_sw_frz
            ocn.T_ML[i, j] = T_sw_frz
            ocn.Ts[1:((ocn.FLDO[i, j] == -1) ? ocn.Nz[i, j] : ocn.FLDO[i, j]-1 ), i, j] .= T_sw_frz
            OC_updateB!(ocn, i, j)
        end

        #=
        old_FLDO=ocn.FLDO[i,j]
        Δb = (old_FLDO == -1 ) ? 0.0 : ocn.b_ML[i, j] - ocn.bs[old_FLDO, i, j]

        if Δb < -3e-6
            println(format("[Qflx2ATM] [AFTER] At {:d}, {:d}: Δb = {:f}. After adjustment.", i, j, Δb))
        end
        =#
    end
end
