
module ShallowWaterModule

    function allocate(datakind::Symbol, dtype::DataType, dims... ; func=Main.zeros)
        if datakind == :local
            return func(dtype, dims...)
        elseif datakind == :shared
            return SharedArray{dtype}(dims...)
        else
            ErrorException("Unknown kind: " * string(datakind)) |> throw
        end
    end


    mutable struct TracerAdv

        ASUM :: Union{AdvectionSpeedUpMatrix, Nothing}
        
        workspace1    :: AbstractArray{Float64, 3}
        workspace2    :: AbstractArray{Float64, 3}
        workspace3    :: AbstractArray{Float64, 3}

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
        
        div           :: AbstractArray{Float64, 3}

        function TracerAdv(;
            gi             :: DisplacedPoleCoordinate.GridInfo,
            Nx             :: Int64,
            Ny             :: Int64,
            Nz             :: Int64,
            NX             :: Int64,  # number of tracers
            Nz_av          :: AbstractArray{Int64, 2},
            mask3          :: AbstractArray{Float64, 3},
            noflux_x_mask3 :: AbstractArray{Float64, 3},
            noflux_y_mask3 :: AbstractArray{Float64, 3},
            Δz             :: AbstractArray{Float64, 3},
            H              :: AbstractArray{Float64, 3},
            datakind       :: Symbol,
            create_mtx     :: Bool,
        )

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

            div     = allocate(datakind, Float64, Nz, Nx, Ny)

            if create_mtx
                ASUM = AdvectionSpeedUpMatrix(;
                                gi = gridinfo,
                                Nx = Nx,
                                Ny = Ny,
                                Nz_bone = Nz_f,
                                Nz = Nz_f_av,
                                mask3 = mask3,
                                noflux_x_mask3 = noflux_x_mask3,
                                noflux_y_mask3 = noflux_y_mask3,
                                Δzs = Δz_f,
                                hs  = H_f,
                )
            else
                ASUM = nothing
            end

            return new(
                ASUM,
                GRAD_bnd_x, GRAD_bnd_y, GRAD_bnd_z,
                CURV_bnd_x, CURV_bnd_y, CURV_bnd_z,
                XFLUX_DEN_x, XFLUX_DEN_y, XFLUX_DEN_z,
                XFLUX_CONV, XFLUX_CONV_h,
                div,
            )
        end
    end

    mutable struct State
        u_c  :: AbstractArray{Float64, 3}
        v_c  :: AbstractArray{Float64, 3}

        u_f  :: AbstractArray{Float64, 3}
        v_f  :: AbstractArray{Float64, 3}

        b  :: AbstractArray{Float64, 3}

        T  :: AbstractArray{Float64, 3} 
        S  :: AbstractArray{Float64, 3} 
        X  :: AbstractArray{Float64, 3}   # passive tracer

 
        # barotropic velocity
        U  :: AbstractArray{Float64, 2}
        V  :: AbstractArray{Float64, 2}
        η  :: AbstractArray{Float64, 2}


        function State(Nx, Ny, Nz_f, Nz_c, NX, datakind)

            u_c = allocate(datakind, Float64, Nx, Ny,   Nz_c)
            v_c = allocate(datakind, Float64, Nx, Ny+1, Nz_c)

            u_f = allocate(datakind, Float64, Nx, Ny,   Nz_f)
            v_f = allocate(datakind, Float64, Nx, Ny+1, Nz_f)

            b = allocate(datakind, Float64, Nx, Ny, Nz_c)
            
            T = allocate(datakind, Float64, Nx, Ny, Nz_f)
            S = allocate(datakind, Float64, Nx, Ny, Nz_f)
            X = allocate(datakind, Float64, Nx, Ny, Nz_f, NX)

            U = allocate(datakind, Float64, Nx, Ny)
            V = allocate(datakind, Float64, Nx, Ny)
            η = allocate(datakind, Float64, Nx, Ny)

            return new(
                u_c, v_c, u_f, v_f,
                b, T, S, X, U, V, η,
            )
            
        end
    end

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

        Nx :: Int64
        Ny :: Int64

        Nz_c :: Int64
        Nz_f :: Int64

        z_bnd_f :: AbstractArray{Float64, 1}
        height_level_counts :: AbstractArray{Int64, 1}
        
        H_c  :: AbstractArray{Float64, 1}
        
        f  :: AbstractArray{Float64, 2}
        ϵ  :: AbstractArray{Float64, 2}

        mask3_f          :: AbstractArray{Float64, 3}
        noflux_x_mask3_f :: AbstractArray{Float64, 3}
        noflux_y_mask3_f :: AbstractArray{Float64, 3}

        function Env(;
            Nx                  :: Int64,
            Ny                  :: Int64,
            Nz_c                :: Int64,
            z_bnd_f             :: AbstractArray{Float64, 1},
            height_level_counts :: AbstractArray{Int64, 1},
            NX                  :: Int64 = 0,
            Nz_f_av             :: Union{AbstractArray{Int64, 3}, Nothing} = nothing,
            H_f                 :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
            Δz_f                :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
            mask3_f             :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
            noflux_x_mask3_f    :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
            noflux_y_mask3_f    :: Union{AbstractArray{Float64, 3}, Nothing} = nothing,
            datakind            :: Symbol,
        )
            z_bnd_f = copy(z_bnd_f)
            Nz_f = length(z_bnd_f) - 1

            height_level_counts = copy(height_level_counts)

            if Nz_c != length(height_level_counts)
                throw(ErrorException("length(height_level_counts) != Nz_c"))
            end

            H_f  = z_bnd_f[1:end-1] - z_bnd_f[2:end]
            H_c  = allocate(datakind, Float64, Nz_c)
           
            idx = 1
            for (k, cnt) in enumerate(height_level_counts)
                H_c[k] = sum(H_f[idx:idx+cnt-1])
                idx += cnt
            end

            # Assume no topography if mask3 is not provided
            if mask3 == nothing
                H_f  = allocate(datakind, Float64, Nx, Ny, Nz_f)
                Δz_f = allocate(datakind, Float64, Nx, Ny, Nz_f-1)
                mask3_f = allocate(Float64, Nx, Ny, Nz_f; func=Main.ones)
                noflux_x_mask3_f = allocate(Float64, Nx, Ny,   Nz_f; func=Main.ones)
                noflux_y_mask3_f = allocate(Float64, Nx, Ny+1, Nz_f; func=Main.ones)
                Nz_f_av = allocate(datakind, Int64, Nx, Ny)

                Nz_f_av .= Nz_f
                H_f[:, 1, 1]  .= z_bnd_f[1:end-1] - z_bnd_f[2:end]
                Δz_f[:, 1, 1] .= (H_f[1:end-1, 1, 1] + H_f[2:end, 1, 1]) / 2.0
                for i=1:Nx, j=1:Ny
                    H_f[:, i, j]  .= H_f[:, 1, 1]
                    Δz_f[:, i, j] .= Δz_f[:, 1, 1]
                end
            end
     

            state = State(Nx, Ny, Nz_f, Nz_c, NX) 


            return new(
                Nx, Ny, Nz_c, Nz_f,
                z_bnd_f, height_level_counts,
                H_c,
                f, ϵ,
                mask3_f,
                noflux_x_mask3_f,
                noflux_y_mask3_f,
            ) 
        end
        

    end

    mutable struct Model
        env     :: Env
        state   :: State
        tcr_adv :: TracerAdv
        dyn_adv :: DynamicAdv

        function Model(;
            Nx, Ny, Nz_c, z_bnd_f, height_level_counts,
            NX = 0,
            Nz_f_av=nothing, H_f=nothing, Δz_f=nothing,
            mask3=nothing, noflux_x_mask3=nothing, noflux_y_mask3=nothing,
            f, ϵ,
            datakind, 
        )
            
            env = Env(;
                Nx = Nx,
                Ny = Ny,
                Nz_c = Nz_c, 
                z_bnd_f = z_bnd_f, 
                height_level_counts = height_level_counts,
                NX = NX,
                Nz_f_av = Nz_f_av,
                H_f = H_f,
                Δz_f = Δz_f,
                mask3 = mask3,
                noflux_x_mask3 = noflux_x_mask3,
                noflux_y_mask3 = noflux_y_mask3,
                datakind,
            )

            state = State(Nx, Ny, Nz_f, Nz_c, NX, datakind)
            
            tcr_adv = TracerAdv(;
                gi = gi,
                Nx = Nx,
                Ny = Ny,
                Nz = Nz_f,
                NX = 2 + NX,  # T, S and the rest
                Nz_av = Nz_f_av
                mask3 = mask3
                
            ) 

            return new(
                env,
                state,
                tcr_adv,
                dyn_adv,
            )



        end
    end


    mutable struct ABIIIObj
        gi :: GridInfo
        
    
        function ABIIIObj(
            ShallowWater
        )
                        
        end
    end


    function ABIII!(
        o :: ABIIIObj
    )

    end
end
