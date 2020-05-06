function calBuoyancyPressure!(
    m :: TmdModel,
)
    @fast_extract m
   

    calFLDOPartition!(m)

    # First mix the FLDO 
    st.b_mixed .= st.b
    @loop_hor m i j let 
        Tmd.mixFLDO!(
            qs   = co.cols.b_mixed[i, j],
            zs   = co.cols.z_bnd_av[i, j],
            hs   = co.cols.Δz_T[i, j],
            q_ML = st.b_ML[i, j],
            FLDO = st.FLDO[i, j],
            FLDO_ratio_top = st.FLDO_ratio_top[i, j],
            FLDO_ratio_bot = st.FLDO_ratio_bot[i, j],
        )
    end

    # Calculate integrated buoyancy
    B_star = getSpace!(co.wksp, :T)
    @. B_star = co.Δz_T * st.b_mixed / 2

    #=
    println("b_mixed: ", st.b_mixed[1, 55, 45])
    println("Δz_T: ",   co.Δz_T[1, 55, 45])
    println("B_star: ", B_star[1, 55, 45])

    throw(ErrorException())
    =#

    # Calculate B_n
    @loop_hor m i j let

        B = co.cols.B[i, j]
        B[1] = B_star[1, i, j]

        for k=2:ev.Nz_av[i, j]
            B[k] = B[k-1] + B_star[k-1, i, j] + B_star[k, i, j]
        end
    end

end

