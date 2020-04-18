mutable struct DynamicAdv    # Adam Bashford


    ASUM :: Union{AdvectionSpeedUpMatrix, Nothing}
    AVGM :: Union{AvgSpeedUpMatrix, Nothing}

    G_idx  :: NamedTuple

    G_bt_u :: AbstractArray{Float64, 3}  # G term for barotropic mode.   (Nx+1, Ny, time)
    G_bt_v :: AbstractArray{Float64, 3}  # G term for barotropic mode.   (Nx, Ny+1, time)

    G_bc_u :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx+1, Ny, time)
    G_bc_v :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx, Ny+1, time)

    ∂b∂x   :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nz_c, Nx+1, Ny)
    ∂b∂y   :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nz_c, Nx, Ny+1)
 
    ∂B∂x   :: AbstractArray{Float64, 2}   # Horizontal gradient of averaged buoyancy  (Nx+1, Ny)
    ∂B∂y   :: AbstractArray{Float64, 2}   # Horizontal gradient of averaged buoyancy  (Nx, Ny+1)
 
    ∂Φ∂x   :: AbstractArray{Float64, 2}   # Horizontal gradient of averaged buoyancy  (Nx+1, Ny)
    ∂Φ∂y   :: AbstractArray{Float64, 2}   # Horizontal gradient of averaged buoyancy  (Nx, Ny+1)
    
    fu     :: AbstractArray{Float64, 3}   # Coriolis force for v             (Nz_c, Nx, Ny+1)
    fv     :: AbstractArray{Float64, 3}   # Coriolis force for u             (Nz_c, Nx+1, Ny)

    ∂u∂x   :: AbstractArray{Float64, 3}   # on U grid (Nz_c, Nx, Ny)
    ∂u∂y   :: AbstractArray{Float64, 3}   # on U grid (Nz_c, Nx, Ny)
    ∂v∂x   :: AbstractArray{Float64, 3}   # on V grid (Nz_c, Nx, Ny+1)
    ∂v∂y   :: AbstractArray{Float64, 3}   # on V grid (Nz_c, Nx, Ny+1)

    u∂u∂x  :: AbstractArray{Float64, 3}   # on U grid (Nz_c, Nx, Ny)
    v∂u∂y  :: AbstractArray{Float64, 3}   # on U grid (Nz_c, Nx, Ny)
    u∂v∂x  :: AbstractArray{Float64, 3}   # on V grid (Nz_c, Nx, Ny+1)
    v∂v∂y  :: AbstractArray{Float64, 3}   # on V grid (Nz_c, Nx, Ny+1)

    DIV_Hu :: AbstractArray{Float64, 2}   # on T grid 

    τx_acc :: AbstractArray{Float64, 2}   # acceleration caused by stress of the first layer
    τy_acc :: AbstractArray{Float64, 2}   # acceleration caused by stress of the first layer

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

        AVGM = AvgSpeedUpMatrix(;
                Nx = Nx,
                Ny = Ny,
                Nz_f = env.Nz_f,
                Nz_c = env.Nz_c,
                mask2 = env.mask2,
                height_level_counts=env.height_level_counts,
                H_f = env.H_f,
                H_c = env.H_c,
        )

 

        #         now   Δt-ago  2Δt-ago
        G_idx = Dict(
            :now           =>  1,
            :one_Δt_ago    =>  2,
            :two_Δt_ago    =>  3,
        )

        new(

            G_idx,

            zeros(Float64,     Nx, Ny,   3), # Adam-Bashforth III 
            zeros(Float64, Nz, Nx, Ny+1, 3), # 3 past

            zeros(Float64, Nz, Nx, Ny, 3),
            zeros(Float64, Nz, Nx, Ny+1, 3),

            zeros(Float64, Nz, Nx, Ny),
            zeros(Float64, Nz, Nx, Ny+1),

            zeros(Float64, Nz, Nx, Ny+1),  # fu is on V grid
            zeros(Float64, Nz, Nx, Ny),    # fv is on U grid

            zeros(Float64, Nz, Nx, Ny),
            zeros(Float64, Nz, Nx, Ny),
            zeros(Float64, Nz, Nx, Ny+1),
            zeros(Float64, Nz, Nx, Ny+1),
 
            zeros(Float64, Nz, Nx, Ny),
            zeros(Float64, Nz, Nx, Ny),
            zeros(Float64, Nz, Nx, Ny+1),
            zeros(Float64, Nz, Nx, Ny+1),
            
            zeros(Float64, Nx, Ny),

            zeros(Float64, Nx, Ny),
            zeros(Float64, Nx, Ny+1),
        )
    end
end


