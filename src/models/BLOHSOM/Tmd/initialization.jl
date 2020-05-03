
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

    updateB!(m)
    updateFLDO!(m)

    if ev.convective_adjustment
        @loop_hor m i j let
            OC_doConvectiveAdjustment!(m, i, j)
        end
    end

end
