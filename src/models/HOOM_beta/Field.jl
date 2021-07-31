mutable struct Field 

    _X       :: AbstractArray{Float64, 1}

    HMXL     :: AbstractArray{Float64, 2}
    TFLUX_DIV_implied      :: AbstractArray{Float64, 2}
    SFLUX_DIV_implied      :: AbstractArray{Float64, 2}

    qflx_T_correction :: AbstractArray{Float64, 3}
    qflx_S_correction :: AbstractArray{Float64, 3}
 
    # Advection related variables
    τx       :: AbstractArray{Float64, 1}
    τy       :: AbstractArray{Float64, 1}

    u        :: AbstractArray{Float64, 1}
    v        :: AbstractArray{Float64, 1}
    w        :: AbstractArray{Float64, 1}

end



mutable struct SugarView
    
    f :: Field

    TEMP     :: AbstractArray{Float64, 3}
    SALT     :: AbstractArray{Float64, 3}

end

