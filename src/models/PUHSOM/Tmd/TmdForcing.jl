mutable struct TmdForcing
 
    X_wk   :: AbstractArray{Float64, 4}
    T_wk   :: AbstractArray{Float64, 3}
    S_wk   :: AbstractArray{Float64, 3}
 
    XSAS_wk   :: AbstractArray{Float64, 3}
    TSAS_wk   :: AbstractArray{Float64, 2}
    SSAS_wk   :: AbstractArray{Float64, 2}

    u_U    :: AbstractArray{Float64, 3}
    v_V    :: AbstractArray{Float64, 3}
    w_W    :: AbstractArray{Float64, 3}

    qflx2atm     :: AbstractArray{Float64, 2}
    qflx2atm_pos :: AbstractArray{Float64, 2}
    qflx2atm_neg :: AbstractArray{Float64, 2}

    qflx_X :: AbstractArray{Float64, 3}
    qflx_T :: AbstractArray{Float64, 2}
    qflx_S :: AbstractArray{Float64, 2}

    τx       :: AbstractArray{Float64, 2}
    τy       :: AbstractArray{Float64, 2}

    nswflx :: AbstractArray{Float64, 2}
    swflx  :: AbstractArray{Float64, 2}
    vsflx  :: AbstractArray{Float64, 2}
    ifrac  :: AbstractArray{Float64, 2}

    h_ML   :: AbstractArray{Float64, 2}


    function TmdForcing(
        env :: TmdEnv,
    )
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz
        NX = env.NX

        X_wk = zeros(Float64, Nz, Nx, Ny, NX)
        T_wk = view(X_wk, :, :, :, 1)
        S_wk = view(X_wk, :, :, :, 2)
 
        XSAS_wk = zeros(Float64, Nx, Ny, NX)
        TSAS_wk = view(XSAS_wk, :, :, 1)
        SSAS_wk = view(XSAS_wk, :, :, 2)
 

        u_U = zeros(Float64, Nz, Nx, Ny)
        v_V = zeros(Float64, Nz, Nx, Ny+1)
        w_W = zeros(Float64, Nz+1, Nx, Ny)
 
        qflx2atm = zeros(Float64, Nx, Ny)
        qflx2atm_pos = zeros(Float64, Nx, Ny)
        qflx2atm_neg = zeros(Float64, Nx, Ny)

        qflx_X  = zeros(Float64, Nx, Ny, NX)
        qflx_T = view(qflx_X, :, :, 1)
        qflx_S = view(qflx_X, :, :, 2)
 
        τx = zeros(Float64, Nx, Ny)
        τy = zeros(Float64, Nx, Ny)
        
        nswflx = zeros(Float64, Nx, Ny)
        swflx  = zeros(Float64, Nx, Ny)
        vsflx  = zeros(Float64, Nx, Ny)
        ifrac  = zeros(Float64, Nx, Ny)
        
        h_ML = zeros(Float64, Nx, Ny)
 
         return new(
            X_wk,
            T_wk,
            S_wk,
         
            XSAS_wk,
            TSAS_wk,
            SSAS_wk,

            u_U,
            v_V,
            w_W,

            qflx2atm,
            qflx2atm_pos,
            qflx2atm_neg,

            qflx_X,
            qflx_T,
            qflx_S,

            τx,
            τy,

            nswflx,
            swflx,
            vsflx,
            ifrac,

            h_ML,
        )
    end
end


