function OC_doConvectiveAdjustment!(
        occ :: OceanColumnCollection,
        i   :: Integer,
        j   :: Integer,
    )

    if_adjust = false

    if occ.bs[i, j, 1] < occ.bs[i, j, 2]

        if_adjust = true

        mixed_T = (occ.hs[1] * occ.Ts[i, j, 1] + occ.hs[2] * occ.Ts[i, j, 2]) / (-occ.zs[3])
        mixed_S = (occ.hs[1] * occ.Ss[i, j, 1] + occ.hs[2] * occ.Ss[i, j, 2]) / (-occ.zs[3])

        occ.Ts[i, j, :] .= mixed_T
        occ.Ss[i, j, :] .= mixed_S 

        OC_updateB!(occ, i, j)

    end 

    return if_adjust
end


