mutable struct TmdDiag
 
    intX             :: AbstractArray{Float64, 3} # Total X content, vertical integrated quantity
    dintXdt          :: AbstractArray{Float64, 3}

    intT             :: AbstractArray{Float64, 2} # Total heat content
    dintTdt          :: AbstractArray{Float64, 2} # Total heat content change rate
    intS             :: AbstractArray{Float64, 2} # Total salt
    dintSdt          :: AbstractArray{Float64, 2} # Total salt change rate

    qflx_X_correction :: AbstractArray{Float64, 3}
    qflx_T_correction :: AbstractArray{Float64, 2}
    qflx_S_correction :: AbstractArray{Float64, 2}

    XFLUX_DIV_implied :: AbstractArray{Float64, 3}
    TFLUX_DIV_implied :: AbstractArray{Float64, 2}
    SFLUX_DIV_implied :: AbstractArray{Float64, 2}


    function TmdDiag(
        env :: TmdEnv,
    )
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz
        NX = env.NX

        intX = zeros(Float64, Nx, Ny, NX)
        intT = view(intX, :, :, 1)
        intS = view(intX, :, :, 2)

        dintXdt = zeros(Float64, Nx, Ny, NX)
        dintTdt = view(dintXdt, :, :, 1)
        dintSdt = view(dintXdt, :, :, 2)
 
        qflx_X_correction = zeros(Float64, Nx, Ny, NX)
        qflx_T_correction = view(qflx_X_correction, :, :, 1)
        qflx_S_correction = view(qflx_X_correction, :, :, 2)

        XFLUX_DIV_implied = zeros(Float64, Nx, Ny, NX)
        TFLUX_DIV_implied = view(XFLUX_DIV_implied, :, :, 1)
        SFLUX_DIV_implied = view(XFLUX_DIV_implied, :, :, 2)

        return new(
            intX,
            dintXdt,

            intT,
            dintTdt,
            intS,
            dintSdt,

            qflx_X_correction,
            qflx_T_correction,
            qflx_S_correction,

            XFLUX_DIV_implied,
            TFLUX_DIV_implied,
            SFLUX_DIV_implied,
        )
    end

end


