function setHQ(;
    occ :: OceanColumnCollection,
    h_ML:: Array{Float64, 2},
    Q_ML:: Array{Float64, 2},
    t   :: Array{Float64, 1},
)

    interpolatePeriodic!(
        time1 :: t,
        data1 :: h_ML,
        time2 :: occ.t,
        data2 :: occ.h_ML 
    )

    interpolatePeriodic!(
        time1 :: t,
        data1 :: Q_ML,
        time2 :: occ.t,
        data2 :: occ.Q_ML 
    )

    Δtt = 2.0 * occ.Δt
    for i = 1:occ.N_ocs
        occ.we[i, :] = (circshift(occ.h_ML[i,:], -1) - circshift(occ.h_ML[i,:], 1)) ./ Δtt
    end

end
