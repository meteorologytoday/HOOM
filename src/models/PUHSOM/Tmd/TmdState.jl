mutable struct TmdState
 
    X         :: AbstractArray{Float64, 4}
    T         :: AbstractArray{Float64, 3}
    S         :: AbstractArray{Float64, 3}
  
    X_ML     :: AbstractArray{Float64, 3}
    T_ML     :: AbstractArray{Float64, 2}
    S_ML     :: AbstractArray{Float64, 2}

    
    ΔX        :: AbstractArray{Float64, 3}
    ΔT        :: AbstractArray{Float64, 2}
    ΔS        :: AbstractArray{Float64, 2}

    dΔXdt     :: AbstractArray{Float64, 3}
    dΔTdt     :: AbstractArray{Float64, 2}
    dΔSdt     :: AbstractArray{Float64, 2}
    
    
    XSAS_wk   :: AbstractArray{Float64, 3}
    TSAS_wk   :: AbstractArray{Float64, 2}
    SSAS_wk   :: AbstractArray{Float64, 2}


    intX             :: AbstractArray{Float64, 3} # Total X content, vertical integrated quantity
    dintXdt          :: AbstractArray{Float64, 3}

    intT             :: AbstractArray{Float64, 2} # Total heat content
    dintTdt          :: AbstractArray{Float64, 2} # Total heat content change rate
    intS             :: AbstractArray{Float64, 2} # Total salt
    dintSdt          :: AbstractArray{Float64, 2} # Total salt change rate


    qflx_X_correction :: AbstractArray{Float64, 3}
    qflx_T_correction :: AbstractArray{Float64, 2}
    qflx_S_correction :: AbstractArray{Float64, 2}


    h_ML      :: AbstractArray{Float64, 2}
    h_MO      :: AbstractArray{Float64, 2}
    
    FLDO           :: AbstractArray{Int64,     2}
    FLDO_ratio_top :: AbstractArray{Float64,   2}
    FLDO_ratio_bot :: AbstractArray{Float64,   2}

    qflx2atm         :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes
    qflx2atm_pos     :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes
    qflx2atm_neg     :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    # Advection related variables
    τx       :: AbstractArray{Float64, 2}
    τy       :: AbstractArray{Float64, 2}

    u_U    :: AbstractArray{Float64, 3}
    v_V    :: AbstractArray{Float64, 3}
    w_W    :: AbstractArray{Float64, 3}


    b        :: AbstractArray{Float64, 3}
    b_ML     :: AbstractArray{Float64, 2}
 

    function TmdState(
        env :: TmdEnv,
    )
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz
        NX = env.NX

        X = zeros(Float64, Nz, Nx, Ny, NX)
        T = view(X, :, :, :, 1)
        S = view(X, :, :, :, 2)


        X_ML = zeros(Float64, Nx, Ny, NX)
        T_ML = view(X_ML, :, :, 1)
        S_ML = view(X_ML, :, :, 2)
        
        ΔX = zeros(Float64, Nx, Ny, NX)
        ΔT = view(ΔX, :, :, 1)
        ΔS = view(ΔX, :, :, 2)
        
        dΔXdt = zeros(Float64, Nx, Ny, NX)
        dΔTdt = view(dΔXdt, :, :, 1)
        dΔSdt = view(dΔXdt, :, :, 2)

    
        XSAS_wk = zeros(Float64, Nx, Ny, NX)
        TSAS_wk = view(XSAS_wk, :, :, 1)
        SSAS_wk = view(XSAS_wk, :, :, 2)
   
        intX = zeros(Float64, Nx, Ny, NX)
        intT = view(intX, :, :, 1)
        intS = view(intX, :, :, 2)

        dintXdt = zeros(Float64, Nx, Ny, NX)
        dintTdt = view(dintXdt, :, :, 1)
        dintSdt = view(dintXdt, :, :, 2)
 
        qflx_X_correction = zeros(Float64, Nx, Ny, NX)
        qflx_T_correction = view(qflx_X_correction, :, :, 1)
        qflx_S_correction = view(qflx_X_correction, :, :, 2)

           
        h_ML = zeros(Float64, Nx, Ny)
        h_MO = zeros(Float64, Nx, Ny)
        
        FLDO = zeros(Int64, Nx, Ny)
        FLDO_ratio_top = zeros(Float64, Nx, Ny)
        FLDO_ratio_bot = zeros(Float64, Nx, Ny)

        qflx2atm = zeros(Float64, Nx, Ny)
        qflx2atm_pos = zeros(Float64, Nx, Ny)
        qflx2atm_neg = zeros(Float64, Nx, Ny)


        τx = zeros(Float64, Nx, Ny)
        τy = zeros(Float64, Nx, Ny)

        u_U = zeros(Float64, Nz, Nx, Ny)
        v_V = zeros(Float64, Nz, Nx, Ny+1)
        w_W = zeros(Float64, Nz+1, Nx, Ny)

        b    = zeros(Float64, Nz, Nx, Ny)
        b_ML = zeros(Float64, Nx, Ny)


        return new(
            X,
            T,
            S,
           
            X_ML,
            T_ML,
            S_ML,
           
            ΔX,
            ΔT,
            ΔS,

            dΔXdt,
            dΔTdt,
            dΔSdt,
            
            XSAS_wk,
            TSAS_wk,
            SSAS_wk,

            intX,
            dintXdt,

            intT,
            dintTdt,
            intS,
            dintSdt,

            qflx_X_correction,
            qflx_T_correction,
            qflx_S_correction,

            h_ML,
            h_MO,
            
            FLDO,
            FLDO_ratio_top,
            FLDO_ratio_bot,

            qflx2atm,
            qflx2atm_pos,
            qflx2atm_neg,

            τx,
            τy,

            u_U,
            v_V,
            w_W,

            b,
            b_ML,
        )
    end

end


