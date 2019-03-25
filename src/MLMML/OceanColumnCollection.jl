mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    Nz       :: Integer           # Number of layers
    
    zs       :: AbstractArray{Float64, 1} # Position of (N+1) grid points
    K_T      :: Float64           # Diffusion coe of temperature
    K_S      :: Float64           # Diffusion coe of salinity

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    b_ML     :: AbstractArray{Float64, 2}
    S_ML     :: AbstractArray{Float64, 2}
    T_ML     :: AbstractArray{Float64, 2}
    h_ML     :: AbstractArray{Float64, 2}

    bs       :: AbstractArray{Float64, 3}
    Ss       :: AbstractArray{Float64, 3}
    Ts       :: AbstractArray{Float64, 3}
    FLDO     :: AbstractArray{Int64, 2}
    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    h_ML_min :: Float64
    h_ML_max :: Float64
    we_max   :: Float64

    # Derived quantities
    N_ocs  :: Integer           # Number of columns
    hs     :: AbstractArray{Float64, 1} # Thickness of layers
    Δzs    :: AbstractArray{Float64, 1} # Δz between layers
    H      :: Integer

    function OceanColumnCollection(;
        Nx     :: Integer,
        Ny     :: Integer,
        zs     :: AbstractArray{Float64, 1},
        Ss     :: Union{AbstractArray{Float64, 1}, Nothing} = nothing,
        Ts     :: Union{AbstractArray{Float64, 1}, Nothing} = nothing,
        K_T    :: Float64,
        K_S    :: Float64,
        S_ML   :: Float64,
        T_ML   :: Float64,
        h_ML   :: Float64,
        h_ML_min :: Float64,
        h_ML_max :: Float64,
        we_max   :: Float64,
        mask   :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    )

        if -h_ML_max < zs[end]
            throw(ErrorException("h_ML_max should not be equal or greater than ocean max depth."))
        end

        zs = copy(zs)
        Nz = length(zs) - 1

        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0

        if mask == nothing
            mask = SharedArray{Float64}(Nx, Ny)
            mask .+= 1.0
        else
            mask = copy(mask)
        end

        mask_idx = (mask .== 0.0)


        _b_ML     = SharedArray{Float64}(Nx, Ny)
        _S_ML     = SharedArray{Float64}(Nx, Ny)
        _T_ML     = SharedArray{Float64}(Nx, Ny)
        _h_ML     = SharedArray{Float64}(Nx, Ny)

        _bs       = SharedArray{Float64}(Nx, Ny, Nz)
        _Ss       = SharedArray{Float64}(Nx, Ny, Nz)
        _Ts       = SharedArray{Float64}(Nx, Ny, Nz)
        _FLDO     = zeros(Int64, Nx, Ny)
        qflx2atm  = SharedArray{Float64}(Nx, Ny)

        _h_ML .= h_ML
        _S_ML .= S_ML
        _T_ML .= T_ML
       
        if Ts != nothing 
            for i=1:Nx, j=1:Ny
                _Ts[i, j, :] = Ts
            end
        end

        if Ss != nothing 
            for i=1:Nx, j=1:Ny
                _Ss[i, j, :] = Ss
            end
        end

        occ = new(
            Nx, Ny, Nz,
            zs, K_T, K_S,
            mask, mask_idx,
            _b_ML, _S_ML, _T_ML, _h_ML,
            _bs,   _Ss,   _Ts,
            _FLDO, qflx2atm,
            h_ML_min, h_ML_max, we_max,
            Nx * Ny, hs, Δzs, zs[1] - zs[end]
        )

        updateB!(occ)
        updateFLDO!(occ)
        

        return occ
    end

end

function copyOCC!(fr_occ::OceanColumnCollection, to_occ::OceanColumnCollection)

    if (fr_occ.Nx, fr_occ.Ny, fr_occ.Nz) != (to_occ.Nx, to_occ.Ny, to_occ.Nz)
        throw(ErrorException("These two OceanColumnCollection have different dimensions."))
    end
    
    to_occ.K_T = fr_occ.K_T
    to_occ.K_S = fr_occ.K_S
    to_occ.h_ML_max = fr_occ.h_ML_max
    to_occ.h_ML_min = fr_occ.h_ML_min
    to_occ.we_max = fr_occ.we_max

    to_occ.b_ML[:, :]      = fr_occ.b_ML
    to_occ.T_ML[:, :]      = fr_occ.T_ML
    to_occ.S_ML[:, :]      = fr_occ.S_ML
    to_occ.h_ML[:, :]      = fr_occ.h_ML
    to_occ.FLDO[:, :]      = fr_occ.FLDO
    to_occ.qflx2atm[:, :]  = fr_occ.qflx2atm
    to_occ.bs[:, :, :]     = fr_occ.bs
    to_occ.Ts[:, :, :]     = fr_occ.Ts
    to_occ.Ss[:, :, :]     = fr_occ.Ss
end


function copyOCC(occ::OceanColumnCollection)
    occ2 = makeBlankOceanColumnCollection(occ.Nx, occ.Ny, occ.zs; mask=mask)
    copyOCC!(occ, occ2)
    
    return occ2
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
        h_ML = -zs[2],
        K_T  = 0.0,
        K_S  = 0.0,
        S_ML = 0.0,
        T_ML = 0.0,
        h_ML_min = -zs[2],
        h_ML_max = -zs[end-1],
        we_max   = 0.0,
        mask = mask,
    )
end


function makeBasicOceanColumnCollection(
    Nx      :: Integer,
    Ny      :: Integer,
    zs      :: AbstractArray{Float64, 1};
    T_slope :: Float64,
    S_slope :: Float64,
    T_ML    :: Float64,
    S_ML    :: Float64,
    h_ML    :: Float64,
    ΔT      :: Float64,
    ΔS      :: Float64,
    K_T     :: Float64,
    K_S     :: Float64,
    h_ML_min:: Float64,
    h_ML_max:: Float64,
    we_max  :: Float64,
    mask :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)
    Ts = zeros(Float64, length(zs)-1)
    Ss = zeros(Float64, length(zs)-1)
    for i = 1:length(Ts)
        z = (zs[i] + zs[i+1]) / 2.0
        if z > -h_ML
            Ts[i] = T_ML
            Ss[i] = S_ML
        else
            Ts[i] = T_ML - ΔT - T_slope * (-z - h_ML)
            Ss[i] = S_ML - ΔS - S_slope * (-z - h_ML)
        end
    end

    return OceanColumnCollection(;
        Nx   = Nx,
        Ny   = Ny,
        zs   = zs,
        Ts   = Ts,
        Ss   = Ss,
        K_T  = K_T,
        K_S  = K_S,
        T_ML = T_ML,
        S_ML = S_ML,
        h_ML = h_ML,
        h_ML_min = h_ML_min,        
        h_ML_max = h_ML_max,
        we_max = we_max,
        mask = mask,
    )
end
