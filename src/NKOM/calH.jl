function calH_dHdT!(ocn::Ocean; Δt::Float64)

    @loop_hor ocn i j let
        
        tmp_H = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1]) * ρc
        ocn.H[i, j], ocn.dHdt[i, j] = tmp_H, (tmp_H - ocn.H[i, j]) / Δt

    end
end

function calH!(ocn::Ocean)

    @loop_hor ocn i j let
        ocn.H[i, j] = OC_getIntegratedTemperature(ocn, i, j; target_z = ocn.cols.zs[i, j][ocn.Nz[i, j]+1]) * ρc
    end
end
