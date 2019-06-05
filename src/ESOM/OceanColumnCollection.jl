mutable struct Workspace
    τx        :: AbstractArray{Float64, 2}
    τy        :: AbstractArray{Float64, 2}
    M1x        :: AbstractArray{Float64, 2}
    M1y        :: AbstractArray{Float64, 2}
    M1x_T1     :: AbstractArray{Float64, 2}
    M1y_T1     :: AbstractArray{Float64, 2}
    M1x_T2     :: AbstractArray{Float64, 2}
    M1y_T2     :: AbstractArray{Float64, 2}
    M1x_S1     :: AbstractArray{Float64, 2}
    M1y_S1     :: AbstractArray{Float64, 2}
    M1x_S2     :: AbstractArray{Float64, 2}
    M1y_S2     :: AbstractArray{Float64, 2}
    div_M1T1   :: AbstractArray{Float64, 2}
    div_M1T2   :: AbstractArray{Float64, 2}
    div_M1S1   :: AbstractArray{Float64, 2}
    div_M1S2   :: AbstractArray{Float64, 2}
    div_M1     :: AbstractArray{Float64, 2}

    function Workspace(
        Nx :: Integer,
        Ny :: Integer,
    )
        τx      = SharedArray{Float64}(Nx, Ny)
        τy      = SharedArray{Float64}(Nx, Ny)
        M1x      = SharedArray{Float64}(Nx, Ny)
        M1y      = SharedArray{Float64}(Nx, Ny)
        M1x_T1   = SharedArray{Float64}(Nx, Ny)
        M1y_T1   = SharedArray{Float64}(Nx, Ny)
        M1x_T2   = SharedArray{Float64}(Nx, Ny)
        M1y_T2   = SharedArray{Float64}(Nx, Ny)
        M1x_S1   = SharedArray{Float64}(Nx, Ny)
        M1y_S1   = SharedArray{Float64}(Nx, Ny)
        M1x_S2   = SharedArray{Float64}(Nx, Ny)
        M1y_S2   = SharedArray{Float64}(Nx, Ny)
        div_M1T1 = SharedArray{Float64}(Nx, Ny)
        div_M1T2 = SharedArray{Float64}(Nx, Ny)
        div_M1S1 = SharedArray{Float64}(Nx, Ny)
        div_M1S2 = SharedArray{Float64}(Nx, Ny)
        div_M1   = SharedArray{Float64}(Nx, Ny)
        
        return new(
            τx, τy,
            M1x, M1y,
            M1x_T1, M1y_T1, M1x_T2, M1y_T2,
            M1x_S1, M1y_S1, M1x_S2, M1y_S2,
            div_M1T1, div_M1T2,
            div_M1S1, div_M1S2,
            div_M1,
        )
    end
end


mutable struct OceanColumnCollection

    gi       :: DisplacedPoleCoordinate.GridInfo
    gi_file  :: AbstractString

    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    
    hs     :: AbstractArray{Float64, 1} # Thickness of layers

    topo     :: AbstractArray{Float64, 2} # Depth of the topography. Negative value if it is underwater

    Kh_T     :: Float64           # Horizontal diffusion coe of temperature
    Kh_S     :: Float64           # Horizontal diffusion coe of salinity
    fs       :: AbstractArray{Float64, 2}
    ϵs       :: AbstractArray{Float64, 2}

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    bs       :: AbstractArray{Float64, 3}
    Ts       :: AbstractArray{Float64, 3}
    Ss       :: AbstractArray{Float64, 3}

    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes
    et_x     :: AbstractArray{Float64, 2} # Ekman transport in x direction
    et_y     :: AbstractArray{Float64, 2} # Ekman transport in y direction

    # Derived quantities
    zs       :: AbstractArray{Float64, 1}
    N_ocs    :: Integer           # Number of columns
    
    wksp     :: Workspace

    function OceanColumnCollection(;
        gridinfo_file :: AbstractString,
        Nx       :: Integer,
        Ny       :: Integer,
        hs       :: AbstractArray{Float64, 1},
        Ts       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        Ss       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        Kh_T     :: Float64,
        Kh_S     :: Float64,
        fs       :: Union{AbstractArray{Float64, 2}, Float64, Nothing} = nothing,
        ϵs       :: Union{AbstractArray{Float64, 2}, Float64},
        mask     :: Union{AbstractArray{Float64, 2}, Nothing},
        topo     :: Union{AbstractArray{Float64, 2}, Nothing},
    )

        if length(hs) != 2 && all(hs .> 0)
            throw(ErrorException("ESOM needs 2 thicknesses. Also all values in `hs` should be positive."))
        end

        _hs = copy(hs)
        _zs = zeros(Float64, 3)
        _zs[:] = [0.0, -hs[1], - (hs[1] + hs[2]) ]

        # ===== [BEGIN] topo, mask =====
       
        _topo = SharedArray{Float64}(Nx, Ny)
        _mask = SharedArray{Float64}(Nx, Ny)

        if topo == nothing
            _topo .= _zs[end]
        else
            _topo[:, :] = topo
        end
        
        # mask =>   inactive = 0, ocn = 1
        if mask == nothing
            _mask .+= 1.0
        else
            _mask[:, :] = mask
        end

        # mask_idx needs to be determined after topography is probably tuned.
        mask_idx = (_mask .== 1.0)

        # ===== [END] topo, mask =====

        # ===== [BEGIN] Set mask if not deep enough =====

        for i=1:Nx, j=1:Ny

            if _mask[i, j] == 0
                continue
            end

            if _topo[i, j] > _zs[2]
                _mask[i, j] = 0
            end
        end
        
        # ===== [END] Set mask if not deep enough =====
        

        # ===== [BEGIN] Column information =====

        _bs       = SharedArray{Float64}(Nx, Ny, 2)
        _Ts       = SharedArray{Float64}(Nx, Ny, 2)
        _Ss       = SharedArray{Float64}(Nx, Ny, 2)
        _qflx2atm = SharedArray{Float64}(Nx, Ny)
        _et_x     = SharedArray{Float64}(Nx, Ny)
        _et_y     = SharedArray{Float64}(Nx, Ny)

        if typeof(Ts) <: AbstractArray{Float64, 3}
            _Ts[:, :, :] = Ts
        elseif typeof(Ts) <: AbstractArray{Float64, 1}
            for i=1:Nx, j=1:Ny
                _Ts[i, j, :] = Ts
            end
        elseif typeof(Ts) <: Float64 
            _Ts .= Ts
        end

        if typeof(Ss) <: AbstractArray{Float64, 3}
            _Ss[:, :, :] = Ss
        elseif typeof(Ss) <: AbstractArray{Float64, 1}
            for i=1:Nx, j=1:Ny
                _Ss[i, j, :] = Ss
            end
        elseif typeof(Ss) <: Float64 
            _Ss .= Ss
        end

        # ===== [END] Column information =====




        # ===== [BEGIN] Mask out data =====

        # Clean up all variables
        for v in [_bs, _Ts, _Ss]
            for i=1:Nx, j=1:Ny
                if _mask[i, j] == 0
                    v[i, j, :] .= NaN
                end
            end
        end

        # ===== [END] Mask out data =====

        # ===== [BEG] GridInfo =====
        mi = ModelMap.MapInfo{Float64}(gridinfo_file)
        gridinfo = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv; angle_unit=:deg)

        # ===== [END] GridInfo =====

        # ===== [BEGIN] fs and ϵs =====

        _fs       = SharedArray{Float64}(Nx, Ny)
        _ϵs       = SharedArray{Float64}(Nx, Ny)

        if typeof(fs) <: AbstractArray{Float64, 2}
            _fs[:, :] = fs
        elseif typeof(fs) <: Float64 
            _fs .= fs
        elseif fs == nothing
           _fs[:, :] = 2 * Ωe * sin.(mi.yc * π / 180.0)
        end

        println("Ωe: ", Ωe)

        if typeof(ϵs) <: AbstractArray{Float64, 2}
            _ϵs[:, :] = ϵs
        elseif typeof(ϵs) <: Float64 
            _ϵs .= ϵs
        end

        # ===== [END] fs and ϵs =====

        occ = new(
            gridinfo,
            gridinfo_file,
            Nx, Ny, hs,
            _topo, Kh_T, Kh_S,
            _fs, _ϵs,
            _mask, mask_idx,
            _bs,   _Ts,   _Ss,
            _qflx2atm, _et_x, _et_y,
            _zs, Nx * Ny, Workspace(Nx, Ny)
        )

        updateB!(occ)

        return occ
    end

end

#=
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
    to_occ.Ts[:, :, :]      = fr_occ.Ts
    to_occ.Ss[:, :, :]      = fr_occ.Ss
    to_occ.FLDO[:, :]       = fr_occ.FLDO
    to_occ.qflx2atm[:, :]   = fr_occ.qflx2atm

    to_occ.h_ML_min[:, :]   = fr_occ.h_ML_min
    to_occ.h_ML_max[:, :]   = fr_occ.h_ML_max
    to_occ.we_max           = fr_occ.we_max

    to_occ.Ts_clim_relax_time = fr_occ.Ts_clim_relax_time
    to_occ.Ss_clim_relax_time = fr_occ.Ss_clim_relax_time

    to_occ.Ts_clim[:, :, :] = fr_occ.Ts_clim
    to_occ.Ss_clim[:, :, :] = fr_occ.Ss_clim


    to_occ.N_ocs            = fr_occ.N_ocs
    to_occ.hs[:, :, :]      = fr_occ.hs
    to_occ.Δzs[:, :, :]     = fr_occ.Δzs

end
=#

#=
function copyOCC(occ::OceanColumnCollection)
    occ2 = makeBlankOceanColumnCollection(occ.Nx, occ.Ny, occ.zs; mask=mask)
    copyOCC!(occ, occ2)
    
    return occ2
end
=#
#=
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
        Ts       = nothing,
        Ss       = nothing,
        K_T      = 0.0,
        K_S      = 0.0,
        T_ML     = 0.0,
        S_ML     = 0.0,
        h_ML     = -zs_bone[2],
        h_ML_min = -zs_bone[2],
        h_ML_max = -zs_bone[end-1],
        we_max   = 0.0,
        mask     = mask,
        topo     = topo,
    )
end
=#
