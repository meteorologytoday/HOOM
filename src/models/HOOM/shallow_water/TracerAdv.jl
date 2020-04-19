mutable struct TracerAdv # Leonard 1976

    ASUM :: Union{AdvectionSpeedUpMatrix, Nothing}

    GRAD_bnd_x :: AbstractArray{Float64, 3}
    GRAD_bnd_y :: AbstractArray{Float64, 3}
    GRAD_bnd_z :: AbstractArray{Float64, 3}

    CURV_x     :: AbstractArray{Float64, 3}
    CURV_y     :: AbstractArray{Float64, 3}
    CURV_z     :: AbstractArray{Float64, 3}

    XFLUX_DEN_x :: AbstractArray{Float64, 4}
    XFLUX_DEN_y :: AbstractArray{Float64, 4}
    XFLUX_DEN_z :: AbstractArray{Float64, 4}

    XFLUX_CONV    :: AbstractArray{Float64, 4}
    XFLUX_CONV_h  :: AbstractArray{Float64, 4}
    
    XFLUX_bot     :: AbstractArray{Float64, 3}
    
    div           :: AbstractArray{Float64, 3}

    
    workspace_heap :: AbstractArray{Float64, 4}
    workspaces     :: Any

    function TracerAdv(env)

        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz_f
        NX = env.NX

        GRAD_bnd_x   = zeros(Float64, Nz, Nx, Ny)
        GRAD_bnd_y   = zeros(Float64, Nz, Nx, Ny+1)
        GRAD_bnd_z   = zeros(Float64, Nz+1, Nx, Ny)

        CURV_x       = zeros(Float64, Nz, Nx, Ny)
        CURV_y       = zeros(Float64, Nz, Nx, Ny)
        CURV_z       = zeros(Float64, Nz, Nx, Ny)

        XFLUX_DEN_x  = zeros(Float64, Nz, Nx, Ny, NX)
        XFLUX_DEN_y  = zeros(Float64, Nz, Nx, Ny+1, NX)
        XFLUX_DEN_z  = zeros(Float64, Nz+1, Nx, Ny, NX)

        XFLUX_CONV   = zeros(Float64, Nz, Nx, Ny, NX)
        XFLUX_CONV_h = zeros(Float64, Nz, Nx, Ny, NX)
        
        XFLUX_bot    = zeros(Float64, Nx, Ny, NX)

        div     = zeros(Float64, Nz, Nx, Ny)

        ASUM = AdvectionSpeedUpMatrix(;
                        gi = env.gi,
                        Nx = Nx,
                        Ny = Ny,
                        Nz_bone = Nz,
                        Nz = env.Nz_av_f,
                        mask3 = env.mask3_f,
                        noflux_x_mask3 = env.noflux_x_mask3_f,
                        noflux_y_mask3 = env.noflux_y_mask3_f,
                        Δzs = env.Δz_f,
                        hs  = env.H_f,
        )


        workspace_heap = zeros(Float64, Nz, Nx, Ny, 3)
        workspaces = []
        for w=1:size(workspace_heap)[4]
           push!(workspaces, view(workspace_heap, :, :, :, w))
        end
        
        return new(
            ASUM,
            GRAD_bnd_x, GRAD_bnd_y, GRAD_bnd_z,
            CURV_x, CURV_y, CURV_z,
            XFLUX_DEN_x, XFLUX_DEN_y, XFLUX_DEN_z,
            XFLUX_CONV, XFLUX_CONV_h, XFLUX_bot,
            div,
            workspace_heap, workspaces,
        )
    end
end


