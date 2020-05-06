function calFLDOPartition!(
    m :: TmdModel
)

    @loop_hor m i j let 

        FLDO = m.state.FLDO[i, j]

        if FLDO == -1
            m.state.FLDO_ratio_top[i, j] = 0.0
            m.state.FLDO_ratio_bot[i, j] = 1.0
        else
            m.state.FLDO_ratio_top[i, j] = (m.state.h_ML[i, j] + m.env.z_bnd_av[FLDO, i, j]) / m.core.Δz_T[FLDO, i, j]
            m.state.FLDO_ratio_bot[i, j] = 1.0 - m.state.FLDO_ratio_top[i, j]
        end
    end
end
