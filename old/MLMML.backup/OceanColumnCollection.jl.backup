mutable struct OceanColumnCollection
    N_ocs    :: Integer           # Number of columns

    Nx    :: Integer           # Number of columns in i direction
    Ny    :: Integer           # Number of columns in j direction

    N        :: Integer           # Number of layers
    zs       :: Array{Float64, 1} # Position of (N+1) grid points
    K        :: Float64           # Diffusion coes
    mask     :: Array{Float64}
    mask_idx :: Any

    # Derived quantities
    hs     :: Array{Float64, 1} # Thickness of layers
    Δzs    :: Array{Float64, 1} # Δz between layers

    
    ocs    :: Array{MLMML.OceanColumn, 2}

    # workspace
    wksp   :: Workspace

    function OceanColumnCollection(;
        Nx     :: Integer,
        Ny     :: Integer,
        N      :: Integer,
        zs     :: Array{Float64, 1},
        bs     :: Array{Float64, 1},
        K      :: Float64,
        b_ML   :: Float64,
        h_ML   :: Float64,
        FLDO   :: Integer,
        mask   :: Union{Array{Float64, 2}, Nothing} = nothing,
    )

        N_ocs = Nx * Ny

        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0

        if mask == nothing
            mask = zeros(Float64, Nx, Ny)
            mask .+= 1.0
        end

        mask_idx = (mask .== 0.0)
        ocs = Array{MLMML.OceanColumn}(undef, Nx, Ny)

        for i=1:Nx, j=1:Ny

            if mask[i, j] == 0.0
                continue
            end

            ocs[i, j] = MLMML.OceanColumn(
                N = N,
                zs = zs,
                bs = bs,
                K  = K,
                b_ML = b_ML,
                h_ML = h_ML,
                FLDO = FLDO,
                hs = hs,
                Δzs = Δzs
            )
        end 

        wksp = Workspace(Nx, Ny)

        return new(N_ocs, Nx, Ny, N, zs, K, mask, mask_idx, hs, Δzs, ocs, wksp)
    end

end


