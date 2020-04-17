mutable struct State
    u_c  :: AbstractArray{Float64, 3}
    v_c  :: AbstractArray{Float64, 3}

    u_f  :: AbstractArray{Float64, 3}
    v_f  :: AbstractArray{Float64, 3}
    w_f  :: AbstractArray{Float64, 3}

    b  :: AbstractArray{Float64, 3}

    X  :: AbstractArray{Float64, 4}   # all tracers, 1 is T, 2 is S
    T  :: AbstractArray{Float64, 3}
    S  :: AbstractArray{Float64, 3} 

    # barotropic velocity
    U  :: AbstractArray{Float64, 2}
    V  :: AbstractArray{Float64, 2}
    η  :: AbstractArray{Float64, 2}

    function State(env, datakind)

        Nx = env.Nx
        Ny = env.Ny
        Nz_c = env.Nz_c
        Nz_f = env.Nz_f

        u_c = allocate(datakind, Float64, Nz_c, Nx+1, Ny)
        v_c = allocate(datakind, Float64, Nz_c, Nx, Ny+1)

        u_f = allocate(datakind, Float64, Nz_f, Nx+1, Ny)
        v_f = allocate(datakind, Float64, Nz_f, Nx, Ny+1)
        w_f = allocate(datakind, Float64, Nz_f+1, Nx, Ny)

        b = allocate(datakind, Float64, Nz_c, Nx, Ny)
        
        NX = 2 + env.NX_passive
        X = allocate(datakind, Float64, Nz_f, Nx, Ny, NX)
        T = view(X, :, :, :, 1)
        S = view(X, :, :, :, 2)

        U = allocate(datakind, Float64, Nx, Ny)
        V = allocate(datakind, Float64, Nx, Ny)
        η = allocate(datakind, Float64, Nx, Ny)

        return new(
            u_c, v_c, u_f, v_f, w_f,
            b, X, T, S, U, V, η,
        )
        
    end
end



