mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    Nz       :: Integer           # Number of layers
    
    zs       :: Array{Float64, 1} # Position of (N+1) grid points
    K        :: Float64           # Diffusion coes

    mask     :: Array{Float64, 2}
    mask_idx :: Any

    b_ML     :: Array{Float64, 2}
    h_ML     :: Array{Float64, 2}

    bs       :: Array{Float64, 3}
    FLDO     :: Array{Float64, 2}

    # Derived quantities
    N_ocs  :: Integer           # Number of columns
    hs     :: Array{Float64, 1} # Thickness of layers
    Δzs    :: Array{Float64, 1} # Δz between layers

    # workspace
    wksp   :: Workspace

    function OceanColumnCollection(;
        Nx     :: Integer,
        Ny     :: Integer,
        Nz     :: Integer,
        zs     :: Array{Float64, 1},
        bs     :: Array{Float64, 1},
        K      :: Float64,
        b_ML   :: Float64,
        h_ML   :: Float64,
        FLDO   :: Integer,
        mask   :: Union{Array{Float64, 2}, Nothing} = nothing,
    )

        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0

        if mask == nothing
            mask = zeros(Float64, Nx, Ny)
            mask .+= 1.0
        end

        mask_idx = (mask .== 0.0)

        h_ML = zeros(Float64, Nx, Ny)
        b_ML = zeros(Float64, Nx, Ny)

        bs   = zeros(Float64, Nx, Ny, Nz)
        FLDO = zeros(Float64, Nx, Ny)

        wksp = Workspace(Nx, Ny)

        return new(
            Nx, Ny, Nz,
            zs, K,
            mask, mask_idx,
            b_ML, h_ML, bs, FLDO,
            Nx * Ny, hs, Δzs,
            wksp
        )
    end

end


