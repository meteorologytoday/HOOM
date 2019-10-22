function calH_dHdT!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let
        
        tmp_H = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1]) * ρc
        ocn.H[i, j], ocn.dHdt[i, j] = tmp_H, (tmp_H - ocn.H[i, j]) / Δt

    end



    tmp_budget = 0.0
    tmp_σ = 0.0
    @loop_hor ocn i j let
        #tmp_budget += (ocn.dHdt[i, j] - ocn.wT[i, j] * ρc ) * ocn.mi.area[i, j]
        #tmp_σ      += ocn.mi.area[i, j]

        tmp_budget += (ocn.dHdt[i, j] - ocn.wT[i, j] * ρc ) * ocn.gi.dσ[i, j]
        tmp_σ      += ocn.gi.dσ[i, j]

    end

    println("total balance check inside NKOM (W / m^2): ", tmp_budget / tmp_σ)
    


end

function calH!(ocn::Ocean)

    @loop_hor ocn i j let
        ocn.H[i, j] = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1]) * ρc
    end
end
