

function cleanQflx2atm!(ocn::Ocean)

    @loop_hor ocn i j let

        ocn.qflx2atm[i, j] = 0.0 

    end

end

function calQflx2atm!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let

        ocn.qflx2atm[i, j] = (T_sw_frz - ocn.T_ML[i, j]) * ρ * c_p * ocn.h_ML[i, j] / Δt 
        ocn.T_ML[i, j] = max(T_sw_frz, ocn.T_ML[i, j])
 
        OC_updateB!(ocn, i, j)
    end



end
