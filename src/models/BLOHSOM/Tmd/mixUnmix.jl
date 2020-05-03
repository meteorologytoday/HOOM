function mixFLDO!(;
    qs   :: AbstractArray{Float64, 1},
    zs   :: AbstractArray{Float64, 1},
    hs   :: AbstractArray{Float64, 1},
    q_ML :: Float64,
    FLDO :: Integer,
    FLDO_ratio_top :: Float64,
    FLDO_ratio_bot :: Float64,
)
    
    Δq = 0.0
    
    if FLDO != -1
    
        Δq = q_ML - qs[FLDO]
        qs[FLDO] =  FLDO_ratio_top * q_ML + FLDO_ratio_bot * qs[FLDO]
        
    end
    
    return Δq
    
end

function unmixFLDOKeepDiff!(;
    qs   :: AbstractArray{Float64, 1},
    zs   :: AbstractArray{Float64, 1},
    hs   :: AbstractArray{Float64, 1},
    h_ML :: Float64,
    FLDO :: Integer,
    Nz   :: Integer,
    Δq   :: Float64,
#    verbose = false
)

    int_layer = 0
    integral = 0.0
    new_q_ML = 0.0

    if FLDO == -1

        for k = 1:Nz
            integral += hs[k] * qs[k]
        end 
        new_q_ML = integral / h_ML
        qs[1:Nz] .= new_q_ML

    else

#        verbose && println("FLDO = ", FLDO)
#        verbose && println("Δq: ", Δq)
        for k = 1:FLDO
            integral += hs[k] * qs[k]
#            verbose && println(k, "; hs: ", hs[k], "; qs: ", qs[k])
        end 
#        verbose && println("integral = ", integral)
#        verbose && println("h_ML = ", h_ML, "; zs[FLDO+1] = ", zs[FLDO+1])

        new_q_ML = (integral - Δq * (h_ML + zs[FLDO+1])) / ( - zs[FLDO+1] )
        qs[1:FLDO] .= new_q_ML
        qs[FLDO]    = new_q_ML - Δq


    end

    return new_q_ML

end



function remixML!(;
    qs   :: AbstractArray{Float64, 1},
    zs   :: AbstractArray{Float64, 1},
    hs   :: AbstractArray{Float64, 1},
    h_ML :: Float64,
    FLDO :: Integer,
    Nz   :: Integer,
)

    # no need to remix if h_ML is shallower than the first layer
    if FLDO == 1
        return qs[1]
    end

    int_layer = 0
    new_q_ML = 0.0

    if FLDO == -1
        int_layer = Nz
    else
        int_layer = FLDO - 1
        new_q_ML += ( h_ML + zs[FLDO] ) * qs[FLDO]
    end

    for k = 1:int_layer
        new_q_ML += hs[k] * qs[k]
    end 

    new_q_ML /= h_ML

    qs[1:int_layer] .= new_q_ML
    return new_q_ML

    return qs[1]

end
