
function initialization!(
    m :: TmdModel
)

    @fast_extract m

    
    @loop_hor m i j let
        st.h_ML[i, j] = boundMLD(
            st.h_ML[i, j];
            h_ML_max=ev.h_ML_max[i, j],
            h_ML_min=ev.h_ML_min[i, j]
        )
    end

#=
    println(m.state.T_ML[35,20])
    println(m.state.S_ML[35,20])
    println(m.state.b_ML[35,20])
=#

    updateB!(m)
    updateFLDO!(m)
#=    println("!!!!!!!!!!!!!!!!!!!")
    println(m.state.T[:,35,20])
    println(m.state.S[:,35,20])
    println(m.state.b[:,35,20])
    println(m.state.FLDO[35,20])
    println(m.state.b_ML[35,20])
    println(m.state.T_ML[35,20])
    println(m.state.h_ML[35,20])=#
    if ev.convective_adjustment
        @loop_hor m i j let
            if OC_doConvectiveAdjustment!(m, i, j) && (i,j) == (35,20)
#            println("ADJUST")
            end
        end
    end

#=
    println("!!!!!!!!AFTER!!!!!!!!!")
    println(m.state.T[:,35,20])
    println(m.state.S[:,35,20])
    println(m.state.b[:,35,20])
=#

end
