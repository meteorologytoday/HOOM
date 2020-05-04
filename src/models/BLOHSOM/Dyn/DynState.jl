mutable struct DynState
    
    # total velocity, buoyancy
    u_total  :: AbstractArray{Float64, 3}
    v_total  :: AbstractArray{Float64, 3}
    B        :: AbstractArray{Float64, 3}

    X  :: AbstractArray{Float64, 4}   # all tracers, 1 is T, 2 is S
    T  :: AbstractArray{Float64, 3}   # dimensions in Nz_f grid
    S  :: AbstractArray{Float64, 3} 

    Φ  :: AbstractArray{Float64, 2}
    
    # barotropic component
    U  :: AbstractArray{Float64, 2}
    V  :: AbstractArray{Float64, 2}

    # Baroclinic component
    u  :: AbstractArray{Float64, 3}
    v  :: AbstractArray{Float64, 3}

    # Exception 
    τx  :: AbstractArray{Float64, 2}  # on T grid
    τy  :: AbstractArray{Float64, 2}  # on T grid

    # Excpetion 2
    G_u :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx+1, Ny, time)
    G_v :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx, Ny+1, time)


    function DynState(env, datakind)

        Nx = env.Nx
        Ny = env.Ny
        Nz_c = env.Nz_c
        Nz_f = env.Nz_f

        u_total =  zeros(Float64, Nx, Ny, Nz_c)
        v_total =  zeros(Float64, Nx, Ny+1, Nz_c)
        B       =  zeros(Float64, Nx, Ny, Nz_c)

        NX = 2 + env.NX_passive
        X =  zeros(Float64, Nx, Ny, Nz_f, NX)
        T = view(X, :, :, :, 1)
        S = view(X, :, :, :, 2)

        Φ =  zeros(Float64, Nx, Ny)
        
        U =  zeros(Float64, Nx, Ny  )
        V =  zeros(Float64, Nx, Ny+1)


        u =  zeros(Float64, Nx, Ny, Nz_c)
        v =  zeros(Float64, Nx, Ny+1, Nz_c)
        b =  zeros(Float64, Nx, Ny, Nz_c)

        τx = copy(Φ)
        τy = copy(Φ)

        G_u = zeros(Float64, Nx, Ny, Nz_c, 3)
        G_v = zeros(Float64, Nx, Ny+1, Nz_c, 3)

        return new(
            u_total, v_total, B,
            X, T, S,
            Φ, U, V,
               u, v,
            τx, τy,
            G_u, G_v,
        )
        
    end
end



