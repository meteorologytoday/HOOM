mutable struct DynEnv
    
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
        
    gi :: PolelikeCoordinate.GridInfo

    Δt :: Float64
    Dh :: Float64

    Nx :: Int64
    Ny :: Int64
    Nz_c :: Int64
    Nz_f :: Int64
    
    z_bnd_f :: AbstractArray{Float64, 1}
    height_level_counts :: AbstractArray{Int64, 1}
    
    H_c  :: AbstractArray{Float64, 1}
    H_f  :: AbstractArray{Float64, 1}
    H_total :: Float64
    Φ_total :: Float64
   
    NX   :: Int64
    NX_passive :: Int64

    mask :: AbstractArray{Float64, 2}     # where mixed layer model is active 




    function DynEnv(;
        gi                  :: PolelikeCoordinate.GridInfo,
        Δt                  :: Float64,
        Dh                  :: Float64,
        Nx                  :: Int64,
        Ny                  :: Int64,
        z_bnd_f             :: AbstractArray{Float64, 1},
        height_level_counts :: AbstractArray{Int64, 1},
        NX_passive          :: Int64 = 0,
        mask                :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    )
 
        z_bnd_f = copy(z_bnd_f)
        Nz_f = length(z_bnd_f) - 1

        H_f, H_c = _helper_calLayerThickness(z_bnd_f, height_level_counts)

        Nz_c = length(H_c)

        # If mask is not provided
        if mask == nothing
            mask = ones(Nx, Ny)
        end

        NX = NX_passive + 2

        H_total = sum(H_f)
        Φ_total = g * H_total
        
        println("H_total: ", H_total)
        println("Φ_total: ", Φ_total)

        return new(
            gi,
            Δt,
            Dh,
            Nx, Ny, Nz_c, Nz_f, 
            z_bnd_f,
            height_level_counts,
            H_c,
            H_f,
            H_total,
            Φ_total,
            NX, NX_passive,
            mask,
        ) 
    end
    

end


