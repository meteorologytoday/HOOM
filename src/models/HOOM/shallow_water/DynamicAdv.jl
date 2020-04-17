mutable struct DynamicAdv    # Adam Bashford

    ASUM :: Union{AdvectionSpeedUpMatrix, Nothing}

    G_bt_u :: AbstractArray{Float64, 3}  # G term for barotropic mode.   (Nx+1, Ny, time)
    G_bt_v :: AbstractArray{Float64, 3}  # G term for barotropic mode.   (Nx, Ny+1, time)

    G_bc_u :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx+1, Ny, time)
    G_bc_v :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx, Ny+1, time)

    ∇b_x   :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nz_c, Nx+1, Ny)
    ∇b_y   :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nz_c, Nx, Ny+1)
    
    fu     :: AbstractArray{Float64, 3}   # Coriolis force for v             (Nz_c, Nx, Ny+1)
    fv     :: AbstractArray{Float64, 3}   # Coriolis force for u             (Nz_c, Nx+1, Ny)

    ∇u     :: AbstractArray{Float64, 4}   # Horizontal gradient of velocity  (Nz_c, Nx, Ny, XY)
    ∇v     :: AbstractArray{Float64, 4}   # Horizontal gradient of velocity  (Nz_c, Nx, Ny, XY)

    V∇u    :: AbstractArray{Float64, 3}   # Advection of velocity            (Nz_c, Nx+1, Ny)
    V∇v    :: AbstractArray{Float64, 3}   # Advection of velocity            (Nz_c, Nx, Ny+1)

    
    function DynamicAdv(env, state)
        
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz_c

        ASUM = AdvectionSpeedUpMatrix(;
                gi = env.gi,
                Nx = Nx,
                Ny = Ny,
                Nz_bone = Nz,
                Nz = env.Nz_av_c,
                mask3 = env.mask3_c,
                noflux_x_mask3 = env.noflux_x_mask3_c,
                noflux_y_mask3 = env.noflux_y_mask3_c,
                Δzs = env.Δz_c,
                hs  = env.H_c,
        )

        new(
            zeros(Float64,     Nx+1, Ny, 3),
            zeros(Float64, Nz, Nx, Ny+1, 3),

            zeros(Float64, Nz, Nx+1, Ny, 3),
            zeros(Float64, Nz, Nx, Ny+1, 3),

            zeros(Float64, Nz, Nx+1, Ny),
            zeros(Float64, Nz, Nx, Ny+1),

            zeros(Float64, Nz, Nx+1, Ny),
            zeros(Float64, Nz, Nx, Ny+1),

            zeros(Float64, Nz, Nx, Ny, 2),
            zeros(Float64, Nz, Nx, Ny, 2),
            
            zeros(Float64, Nz, Nx+1, Ny),
            zeros(Float64, Nz, Nx, Ny+1),
        )
    end
end


