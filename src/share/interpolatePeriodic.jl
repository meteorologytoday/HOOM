"""
    interpolatePeriodic!(;
        time1  :: AbstractArray{Float64, 1},
        data1  :: AbstractArray{Float64, 2},
        time2  :: AbstractArray{Float64, 1},
        data2  :: AbstractArray{Float64, 2},
        period :: Float64,
    )

# Description
This function interpolate a series of coarse periodic data into
finer periodic data. Result would be stored in `data2`.

`period` is the length of period of the input periodic data. 
`time1` and `time2` are assumed to be monotonically increaseing and
all of the number should be in the inveral [0, `period`].

"""
function interpolatePeriodic!(;
    time1  :: AbstractArray{Float64, 1},
    data1  :: AbstractArray{Float64, 2}, # first index: space; second index: time
    time2  :: AbstractArray{Float64, 1},
    data2  :: AbstractArray{Float64, 2}, # first index: space; second index: time
    period :: Float64,
)

    interpolate_mtx = _interpolateMatrix(time1, time2, period)

    for i = 1:size(data1)[1]
        data2[i, :] = interpolate_mtx * data1[i, :]
    end

end

function _interpolateMatrix(
    time1 :: AbstractArray{Float64, 1},
    time2 :: AbstractArray{Float64, 1},
    period :: Float64,
)

    N_ext = length(time1)+2
    time1_ext = zeros(Float64, N_ext)

    time1_ext[1]       = time1[end] - period

    if length(time1) > 1
        time1_ext[2:end-1] = time1[1:end]
    end

    time1_ext[end]     = time1[1] + period

    mtx = spzeros(Float64, length(time2), length(time1))
   
    if length(time1) == 1
        mtx .= 1.0
        return mtx
    end
    
    # The algorithm below assume length(time1) > 1
    i_2 = 1
    done = false
    for i = 1:N_ext-1

        ta = time1_ext[i]
        tb = time1_ext[i+1]
        while ta <= time2[i_2] && time2[i_2] < tb
            Δt = tb - ta
            ia = (i   ==     1  ) ? length(time1) : i-1
            ib = (i   == N_ext-1) ? 1             : i
            
            mtx[i_2, ia] = (tb - time2[i_2]) / Δt
            mtx[i_2, ib] = (time2[i_2] - ta) / Δt

            #println(tb - time2[i_2], "; ", time2[i_2] - ta, "; ", Δt, "; ", ia, " , ", ib)

            if i_2 < length(time2)
                i_2 += 1
            else
                done = true
                break
            end
        end
    
        if done
            break
        end
    end

    if !done 
        throw(ErrorException("Not all element of time2 is within [0, period]"))
    end

    return mtx
end
