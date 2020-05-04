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


    h_ML      :: AbstractArray{Float64, 2}
    h_MO      :: AbstractArray{Float64, 2}
    
    FLDO           :: AbstractArray{Int64,     2}
    FLDO_ratio_top :: AbstractArray{Float64,   2}
    FLDO_ratio_bot :: AbstractArray{Float64,   2}

    b        :: AbstractArray{Float64, 3}
    b_ML     :: AbstractArray{Float64, 2}
 
    b_mixed  :: AbstractArray{Float64, 3}
    B        :: AbstractArray{Float64, 3}
 
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

       
        h_ML = zeros(Float64, Nx, Ny)
        h_MO = zeros(Float64, Nx, Ny)
        
        FLDO = zeros(Int64, Nx, Ny)
        FLDO_ratio_top = zeros(Float64, Nx, Ny)
        FLDO_ratio_bot = zeros(Float64, Nx, Ny)

        b    = zeros(Float64, Nz, Nx, Ny)
        b_ML = zeros(Float64, Nx, Ny)
        
        b_mixed = zeros(Float64, Nz, Nx, Ny)
        B       = zeros(Float64, Nz, Nx, Ny)

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

            h_ML,
            h_MO,
            
            FLDO,
            FLDO_ratio_top,
            FLDO_ratio_bot,

            b,
            b_ML,

            b_mixed,
            B,
        )
    end

end


