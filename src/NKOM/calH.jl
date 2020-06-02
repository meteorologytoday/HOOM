function calTEMP_dTEMPdt!(ocn::Ocean; Δt::Float64)

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

function calSALT_dSALTdt!(ocn::Ocean; Δt::Float64)

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