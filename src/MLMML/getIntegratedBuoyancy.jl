function OC_getIntegratedTemperature(
    occ      :: OceanColumnCollection,
    i        :: Integer,
    j        :: Integer;
    target_z :: Float64 = NaN,
)

    if isnan(target_z)
        target_z = occ.zs_vw[occ.Nz[i, j]+1]
    end

    return getIntegratedQuantity(
        zs       = occ.zs_vw[i, j],
        qs       = occ.Ts_vw[i, j],
        q_ML     = occ.T_ML[i, j],
        h_ML     = occ.h_ML[i, j],
        Nz       = occ.Nz[i, j],
        target_z = target_z,
    )
end




function OC_getIntegratedSalinity(
    occ      :: OceanColumnCollection,
    i        :: Integer,
    j        :: Integer;
    target_z :: Float64 = NaN,
)

    if isnan(target_z)
        target_z = occ.zs_vw[occ.Nz[i, j]+1]
    end

    return getIntegratedQuantity(
        zs       = occ.zs_vw[i, j],
        qs       = occ.Ss_vw[i, j],
        q_ML     = occ.S_ML[i, j],
        h_ML     = occ.h_ML[i, j],
        Nz       = occ.Nz[i, j],
        target_z = target_z,
    )
end



function OC_getIntegratedBuoyancy(
    occ      :: OceanColumnCollection,
    i        :: Integer,
    j        :: Integer;
    target_z :: Float64 = NaN,
)

    if isnan(target_z)
        target_z = occ.zs_vw[i, j][occ.Nz[i, j]+1]
    end

    return getIntegratedQuantity(
        zs       = occ.zs_vw[i, j],
        qs       = occ.bs_vw[i, j],
        q_ML     = occ.b_ML[i, j],
        h_ML     = occ.h_ML[i, j],
        Nz       = occ.Nz[i, j],
        target_z = target_z,
    )
end


function getIntegratedQuantity(;
    zs       :: AbstractArray{Float64,1},
    qs       :: AbstractArray{Float64,1},
    q_ML     :: Float64,
    h_ML     :: Float64,
    Nz       :: Integer,
    target_z :: Float64,
)

    if target_z < zs[Nz+1]
        throw(ErrorException(
            format("target_z ({:f}) cannot be deeper than the minimum of zs ({:f}).", target_z, zs[Nz])
        ))
    end


    # Integrate mixed layer
    if -target_z < h_ML
        return q_ML * ( - target_z )
    end

    sum_q = 0.0
    sum_q += h_ML * q_ML


    # Test if entire ocean column is mixed layer
    FLDO = getFLDO(zs=zs, h_ML=h_ML, Nz=Nz)
    if FLDO == -1
        return sum_q
    end

    # Integrate FLDO
    if target_z > zs[FLDO+1]
        sum_q += qs[FLDO] * ( (-h_ML) - target_z)
        return sum_q
    end
    
    sum_q += qs[FLDO] * ( (-h_ML) - zs[FLDO+1]) 

    # Integrate rest layers
    if FLDO < Nz
        for i = FLDO+1 : Nz
            if target_z < zs[i+1]
                sum_q += qs[i] * (zs[i] - zs[i+1])
            else
                sum_q += qs[i] * (zs[i] - target_z)
                return sum_q
            end
        end
    else
        return sum_q
    end

end

