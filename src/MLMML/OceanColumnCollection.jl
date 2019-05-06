mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    Nz_bone  :: Integer           # Number of layers  in k direction
    
    zs_bone  :: AbstractArray{Float64, 1} # Unmasked zs bone
    topo     :: AbstractArray{Float64, 2} # Depth of the topography. Negative value if it is underwater
    zs       :: AbstractArray{Float64, 3} # Actuall zs coordinate masked by topo
    Nz       :: AbstractArray{Float64, 2} # Number of layers that is active

    K_T      :: Float64           # Diffusion coe of temperature
    K_S      :: Float64           # Diffusion coe of salinity

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    b_ML     :: AbstractArray{Float64, 2}
    T_ML     :: AbstractArray{Float64, 2}
    S_ML     :: AbstractArray{Float64, 2}
    h_ML     :: AbstractArray{Float64, 2}

    bs       :: AbstractArray{Float64, 3}
    Ts       :: AbstractArray{Float64, 3}
    Ss       :: AbstractArray{Float64, 3}
    FLDO     :: AbstractArray{Int64, 2}
    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    h_ML_min :: Float64
    h_ML_max :: Float64
    we_max   :: Float64

    # Derived quantities
    N_ocs  :: Integer           # Number of columns
    hs     :: AbstractArray{Float64, 3} # Thickness of layers
    Δzs    :: AbstractArray{Float64, 3} # Δz between layers

    function OceanColumnCollection(;
        Nx       :: Integer,
        Ny       :: Integer,
        zs_bone  :: AbstractArray{Float64, 1},
        Ts       :: Union{AbstractArray{Float64, 1}, Nothing} = nothing,
        Ss       :: Union{AbstractArray{Float64, 1}, Nothing} = nothing,
        K_T      :: Float64,
        K_S      :: Float64,
        T_ML     :: Float64,
        S_ML     :: Float64,
        h_ML     :: Float64,
        h_ML_min :: Float64,
        h_ML_max :: Float64,
        we_max   :: Float64,
        mask     :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
        topo     :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    )

        if - h_ML_max < zs[end]
            throw(ErrorException("h_ML_max should not be equal or greater than ocean max depth."))
        end

        # Dealing with z coordinate
        zs_bone = copy(zs_bone)
        Nz_bone = length(zs_bone) - 1

        zs  = SharedArray{Float64}(Nx, Ny, Nz_bone + 1)
        Nz  = SharedArray{Float64}(Nx, Ny)
        hs  = SharedArray{Float64}(Nx, Ny, Nz_bone)
        Δzs = SharedArray{Float64}(Nx, Ny, Nz_bone - 1)

        zs  .= NaN
        Nz  .= NaN
        hs  .= NaN
        Δzs .= NaN

        if topo == nothing
            topo = zeros(Float64, Nx, Ny)
            topo .= zs_bone[end]
        else
            topo = copy(topo)
        end


        for i=1:Nx, j=1:Ny

            # Determine Nz
            for k=2:length(zs_bone)
                if zs_bone[k] <= topo[i, j]
                    Nz[i, j] = k-1
                    break
                end
            end

            # Construct vertical coordinate
            zs[i, j, 1:Nz] = zs_bone[1:Nz]
            zs[i, j, Nz+1] = topo[i, j]

            # Construct thickness of each layer
            hs[i, j, :]  = zs[i, j, 1:end-1] - zs[i, j, 2:end]
            Δzs[i, j, :] = (hs[i, j, 1:end-1] + hs[i, j, 2:end]) / 2.0
            
        end


        if mask == nothing
            mask = SharedArray{Float64}(Nx, Ny)
            mask .+= 1.0
        else
            mask = copy(mask)
        end

        mask_idx = (mask .== 0.0)

        _b_ML     = SharedArray{Float64}(Nx, Ny)
        _T_ML     = SharedArray{Float64}(Nx, Ny)
        _S_ML     = SharedArray{Float64}(Nx, Ny)
        _h_ML     = SharedArray{Float64}(Nx, Ny)

        _bs       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _Ts       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _Ss       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _FLDO     = zeros(Int64, Nx, Ny)
        qflx2atm  = SharedArray{Float64}(Nx, Ny)

        _h_ML .= h_ML
        _T_ML .= T_ML
        _S_ML .= S_ML
       
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
            Nx, Ny, Nz_bone,
            zs_bone, topo, zs, Nz,
            K_T, K_S,
            mask, mask_idx,
            _b_ML, _T_ML, _S_ML, _h_ML,
            _bs,   _Ts,   _Ss,
            _FLDO, qflx2atm,
            h_ML_min, h_ML_max, we_max,
            Nx * Ny, hs, Δzs,
        )

        updateB!(occ)
        updateFLDO!(occ)
        

        return occ
    end

end

function copyOCC!(fr_occ::OceanColumnCollection, to_occ::OceanColumnCollection)

    if (fr_occ.Nx, fr_occ.Ny, fr_occ.Nz_bone) != (to_occ.Nx, to_occ.Ny, to_occ.Nz_bone)
        throw(ErrorException("These two OceanColumnCollection have different dimensions."))
    end

    to_occ.zs_bone[:]       = fr_occ.zs_bone
    to_occ.topo[:, :]       = fr_occ.topo
    to_occ.zs[:, :, :]      = fr_occ.zs
    to_occ.Nz[:, :]         = fr_occ.Nz
  
    to_occ.K_T              = fr_occ.K_T
    to_occ.K_S              = fr_occ.K_S

    to_occ.mask[:, :]       = fr_occ.mask
    to_occ.mask_idx[:, :]   = fr_occ.mask_idx

    to_occ.b_ML[:, :]       = fr_occ.b_ML
    to_occ.T_ML[:, :]       = fr_occ.T_ML
    to_occ.S_ML[:, :]       = fr_occ.S_ML
    to_occ.h_ML[:, :]       = fr_occ.h_ML

    to_occ.bs[:, :, :]      = fr_occ.bs
    to_occ.Ss[:, :, :]      = fr_occ.Ss
    to_occ.Ts[:, :, :]      = fr_occ.Ts
    to_occ.FLDO[:, :]       = fr_occ.FLDO
    to_occ.qflx2atm[:, :]   = fr_occ.qflx2atm

    to_occ.h_ML_min         = fr_occ.h_ML_min
    to_occ.h_ML_max         = fr_occ.h_ML_max
    to_occ.we_max           = fr_occ.we_max

    to_occ.N_ocs            = fr_occ.N_ocs
    to_occ.hs[:, :, :]      = fr_occ.hs
    to_occ.Δzs[:, :, :]     = fr_occ.Δzs
end

function copyOCC(occ::OceanColumnCollection)
    occ2 = makeBlankOceanColumnCollection(occ.Nx, occ.Ny, occ.zs; mask=mask)
    copyOCC!(occ, occ2)
    
    return occ2
end


function makeBlankOceanColumnCollection(
    Nx      :: Integer,
    Ny      :: Integer,
    zs_bone :: AbstractArray{Float64, 1};
    mask    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    topo    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)

    return OceanColumnCollection(;
        Nx       = Nx,
        Ny       = Ny,
        zs_bone  = zs_bone,
        Ss       = nothing,
        Ts       = nothing,
        K_T      = 0.0,
        K_S      = 0.0,
        S_ML     = 0.0,
        T_ML     = 0.0,
        h_ML     = -zs[2],
        h_ML_min = -zs[2],
        h_ML_max = -zs[end-1],
        we_max   = 0.0,
        mask     = mask,
        topo     = topo,
    )
end


function makeBasicOceanColumnCollection(
    Nx      :: Integer,
    Ny      :: Integer,
    zs_bone :: AbstractArray{Float64, 1};
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
    mask    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    topo    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)

    
    Ts = zeros(Float64, length(zs_bone)-1)
    Ss = zeros(Float64, length(zs_bone)-1)
    for i = 1:length(Ts)
        z = (zs_bone[i] + zs_bone[i+1]) / 2.0
        if z > -h_ML
            Ts[i] = T_ML
            Ss[i] = S_ML
        else
            Ts[i] = T_ML - ΔT - T_slope * (-z - h_ML)
            Ss[i] = S_ML - ΔS - S_slope * (-z - h_ML)
        end
    end

    return OceanColumnCollection(;
        Nx        = Nx,
        Ny        = Ny,
        zs_bone   = zs_bone,
        Ts        = Ts,
        Ss        = Ss,
        K_T       = K_T,
        K_S       = K_S,
        T_ML      = T_ML,
        S_ML      = S_ML,
        h_ML      = h_ML,
        h_ML_min  = h_ML_min,        
        h_ML_max  = h_ML_max,
        we_max    = we_max,
        mask      = mask,
        topo      = topo,
    )
end
