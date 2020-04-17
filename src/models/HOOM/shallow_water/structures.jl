
mutable struct Env
    
    #
    # This shallow water model is fixed height
    # and advect high vertical resultion while
    # horizontal is split into coarser grids
    #
    # c : coarse vertical grid
    # f : fine   vertical grid
    #
    # ----------------+---+---+
    # variable        | c | f |
    # ----------------+---+---+
    # u, v            | v | v |
    # b               | v |   |
    # T, S, X         |   | v |
    # H               | v |   |
    #
    # Variables are arranged in Arakawa-C grid
    # 
        
    gi :: DisplacedPoleCoordinate.GridInfo

    Nx :: Int64
    Ny :: Int64

    Nz_c :: Int64
    Nz_f :: Int64
    Nz_av_f :: AbstractArray{Int64,2}

    z_bnd_f :: AbstractArray{Float64, 1}
    height_level_counts :: AbstractArray{Int64, 1}
    
    H_c  :: AbstractArray{Float64, 1}
    
    H_f  :: AbstractArray{Float64, 3}
    Δz_f :: AbstractArray{Float64, 3}
   
    NX   :: Int64
    NX_passive :: Int64

    Dh   :: AbstractArray{Float64, 1}
    Dv   :: AbstractArray{Float64, 1}

    #f  :: AbstractArray{Float64, 2}
    #ϵ  :: AbstractArray{Float64, 2}
    f
    ϵ


    mask2           :: AbstractArray{Float64, 2}
    loop_idx        :: AbstractArray{Float64, 2}


    mask3_f          :: AbstractArray{Float64, 3}
    noflux_x_mask3_f :: AbstractArray{Float64, 3}
    noflux_y_mask3_f :: AbstractArray{Float64, 3}

    datakind         :: Symbol




    function Env(;
        gi,
        Nx                  :: Int64,
        Ny                  :: Int64,
        z_bnd_f             :: AbstractArray{Float64, 1},
        height_level_counts :: AbstractArray{Int64, 1},
        Dh                  :: AbstractArray{Float64, 1},
        Dv                  :: AbstractArray{Float64, 1},
        NX_passive          :: Int64 = 0,
        Nz_f_av             :: Union{AbstractArray{Int64, 3}, Nothing} = nothing,
        H_f                 :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
        Δz_f                :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
        mask3_f             :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
        noflux_x_mask3_f    :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
        noflux_y_mask3_f    :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
        datakind            :: Symbol,
        f,
        ϵ,
    )

        if ! (datakind in (:local, :shared))
            throw(ErrorException("Datakind not valid. Should be :local or :shared"))
        end

        z_bnd_f = copy(z_bnd_f)
        Nz_f = length(z_bnd_f) - 1

        height_level_counts = copy(height_level_counts)
        Nz_c =  length(height_level_counts)

        if sum(height_level_counts) != Nz_f
            throw(ErrorException("sum(height_level_counts) != length(z_bnd_f) - 1"))
        end

        H_f  = z_bnd_f[1:end-1] - z_bnd_f[2:end]
        H_c  = allocate(datakind, Float64, Nz_c)
       
        idx = 1
        for (k, cnt) in enumerate(height_level_counts)
            H_c[k] = sum(H_f[idx:idx+cnt-1])
            idx += cnt
        end

        # Assume no topography if mask3 is not provided
        if mask3_f == nothing
            H_f  = zeros(Float64, Nz_f, Nx, Ny)
            Δz_f = zeros(Float64, Nz_f - 1, Nx, Ny)
            mask3_f = ones(Nz_f, Nx, Ny)
            noflux_x_mask3_f = ones(Float64, Nz_f, Nx+1, Ny)
            noflux_y_mask3_f = ones(Float64, Nz_f, Nx, Ny+1)
            Nz_f_av = allocate(datakind, Int64, Nx, Ny)

            noflux_y_mask3_f[:, :,   1] .= 0.0
            noflux_y_mask3_f[:, :, end] .= 0.0
            


            Nz_f_av .= Nz_f
            H_f[:, 1, 1]  .= z_bnd_f[1:end-1] - z_bnd_f[2:end]
            Δz_f[:, 1, 1] .= (H_f[1:end-1, 1, 1] + H_f[2:end, 1, 1]) / 2.0
            for i=1:Nx, j=1:Ny
                H_f[:, i, j]  .= H_f[:, 1, 1]
                Δz_f[:, i, j] .= Δz_f[:, 1, 1]
            end
        end


        mask2 = mask3_f[1, :, :]
        loop_idx = zeros(Float64, 2, Int64(sum(mask2)))
        

        NX = NX_passive + 2

        if ! ( length(Dv) == length(Dh) == NX )
            throw(ErrorException("Diffusion coefficients do not match the number of tracers. T and S are included too."))
        end
        return new(
            gi,
            Nx, Ny, Nz_c, Nz_f, Nz_f_av,
            z_bnd_f, height_level_counts,
            H_c,
            H_f,
            Δz_f,
            NX, NX_passive,
            Dh, Dv,
            f, ϵ,
            mask2,
            loop_idx,
            mask3_f,
            noflux_x_mask3_f,
            noflux_y_mask3_f,
            datakind,
        ) 
    end
    

end

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
    function TracerAdv(env, state, datakind)

        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz_f
        NX = env.NX

        GRAD_bnd_x   = allocate(datakind, Float64, Nz, Nx+1, Ny)
        GRAD_bnd_y   = allocate(datakind, Float64, Nz, Nx, Ny+1)
        GRAD_bnd_z   = allocate(datakind, Float64, Nz+1, Nx, Ny)

        CURV_x       = allocate(datakind, Float64, Nz, Nx, Ny)
        CURV_y       = allocate(datakind, Float64, Nz, Nx, Ny)
        CURV_z       = allocate(datakind, Float64, Nz, Nx, Ny)

        XFLUX_DEN_x  = allocate(datakind, Float64, Nz, Nx+1, Ny, NX)
        XFLUX_DEN_y  = allocate(datakind, Float64, Nz, Nx, Ny+1, NX)
        XFLUX_DEN_z  = allocate(datakind, Float64, Nz+1, Nx, Ny, NX)

        XFLUX_CONV   = allocate(datakind, Float64, Nz, Nx, Ny, NX)
        XFLUX_CONV_h = allocate(datakind, Float64, Nz, Nx, Ny, NX)
        
        XFLUX_bot    = allocate(datakind, Float64, Nx, Ny, NX)

        div     = allocate(datakind, Float64, Nz, Nx, Ny)

        if env.datakind != :shared

            println(size(env.noflux_x_mask3_f))
            println(Nx,",",Ny)
            println(env.Nz_f)
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

        else

            ASUM = nothing

        end

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

mutable struct DynamicAdv    # Adam Bashford
    
end


mutable struct Model

    env     :: Env
    state   :: State
    tcr_adv :: TracerAdv
    dyn_adv :: DynamicAdv

    function Model(;
        gi,
        z_bnd_f, height_level_counts,
        f, ϵ,
        Dh,
        Dv,
        NX_passive = 0,
        Nz_f_av=nothing, H_f=nothing, Δz_f=nothing,
        mask3_f=nothing, noflux_x_mask3_f=nothing, noflux_y_mask3_f=nothing,
    )

      
        datakind=:local 

        env = Env(;
            gi = gi,
            Nx = gi.Nx,
            Ny = gi.Ny,
            z_bnd_f = z_bnd_f, 
            height_level_counts = height_level_counts,
            NX_passive = NX_passive,
            Dh=Dh,
            Dv=Dv,
            Nz_f_av = Nz_f_av,
            H_f = H_f,
            Δz_f = Δz_f,
            mask3_f = mask3_f,
            noflux_x_mask3_f = noflux_x_mask3_f,
            noflux_y_mask3_f = noflux_y_mask3_f,
            datakind=datakind,
            f=f,
            ϵ=ϵ,
        )

        state = State(env, datakind)
        tcr_adv = TracerAdv(env, state, datakind)
        dyn_adv = DynamicAdv()
        return new(
            env,
            state,
            tcr_adv,
            dyn_adv,
        )

    end
end


