mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    Nz       :: Integer           # Number of layers
    
    zs       :: AbstractArray{Float64, 1} # Position of (N+1) grid points
    K        :: Float64           # Diffusion coes

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    sst      :: AbstractArray{Float64, 2}
    b_ML     :: AbstractArray{Float64, 2}
    h_ML     :: AbstractArray{Float64, 2}

    bs       :: AbstractArray{Float64, 3}
    FLDO     :: AbstractArray{Int64, 2}
    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    # Derived quantities
    N_ocs  :: Integer           # Number of columns
    hs     :: AbstractArray{Float64, 1} # Thickness of layers
    Δzs    :: AbstractArray{Float64, 1} # Δz between layers

    function OceanColumnCollection(;
        Nx     :: Integer,
        Ny     :: Integer,
        zs     :: AbstractArray{Float64, 1},
        bs     :: Union{AbstractArray{Float64, 1}, Nothing} = nothing,
        K      :: Float64 = 0.0,
        b_ML   :: Float64 = 0.0,
        h_ML   :: Float64 = 0.0,
        mask   :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    )

        Nz = length(zs) - 1

        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0

        if mask == nothing
            mask = zeros(Float64, Nx, Ny)
            mask .+= 1.0
        end
        mask_idx = (mask .== 0.0)

        _h_ML     = zeros(Float64, Nx, Ny)
        _b_ML     = zeros(Float64, Nx, Ny)

        _bs       = zeros(Float64, Nx, Ny, Nz)
        _FLDO     = zeros(Int64, Nx, Ny)
        qflx2atm  = zeros(Float64, Nx, Ny)
        sst       = zeros(Float64, Nx, Ny)

        _h_ML .= h_ML
        _b_ML .= b_ML
        _FLDO .= getFLDO(zs=zs, h_ML=h_ML)
        sst   .= b2T(b_ML)
       
        if bs != nothing 
            for i=1:Nx, j=1:Ny
                _bs[i, j, :] = bs
            end
        end


        occ = new(
            Nx, Ny, Nz,
            zs, K,
            mask, mask_idx,
            sst, _b_ML, _h_ML, _bs, _FLDO, qflx2atm,
            Nx * Ny, hs, Δzs
        )

        return occ
    end

end

function makeBlankOceanColumnCollection(
    Nx   :: Integer,
    Ny   :: Integer,
    zs   :: AbstractArray{Float64, 1};
    mask :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)

    return OceanColumnCollection(;
        Nx   = Nx,
        Ny   = Ny,
        zs   = zs,
        h_ML = h_ML_min,
        mask = mask,
    )
end


function makeBasicOceanColumnCollection(
    Nx      :: Integer,
    Ny      :: Integer,
    zs      :: AbstractArray{Float64, 1};
    b_slope :: Float64 = 30.0 / 5000.0 * g * α,
    b_ML    :: Float64 = 1.0,
    h_ML    :: Float64 = h_ML_min,
    Δb      :: Float64 = 0.0,
    K       :: Float64 = 1e-5,
    mask :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)
    bs = zeros(Float64, length(zs)-1)
    for i = 1:length(bs)
        z = (zs[i] + zs[i+1]) / 2.0
        if z > -h_ML
            bs[i] = b_ML
        else
            bs[i] = b_ML - Δb - b_slope * (-z - h_ML)
        end
    end

    return OceanColumnCollection(;
        Nx   = Nx,
        Ny   = Ny,
        zs   = zs,
        bs   = bs,
        K    = K,
        b_ML = b_ML,
        h_ML = h_ML,
        mask = mask,
    )
end
