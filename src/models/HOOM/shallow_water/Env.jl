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


