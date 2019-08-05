mutable struct OceanColumnCollection

    Nx       :: Integer           # Number of columns in x direction
    Ny       :: Integer           # Number of columns in y direction
    Nz       :: Integer           # Number of columns in z direction
    
    hs     :: AbstractArray{Float64, 1} # Thickness of layers

    topo     :: AbstractArray{Float64, 2} # Depth of the topography. Negative value if it is underwater

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    bs       :: AbstractArray{Float64, 3}
    Ts       :: AbstractArray{Float64, 3}
    Ss       :: AbstractArray{Float64, 3}
    T_ML     :: AbstractArray{Float64, 2}
    S_ML     :: AbstractArray{Float64, 2}
    FLDO     :: AbstractArray{Int64, 2}
    h_ML_min :: AbstractArray{Float64, 2}
    h_ML_max :: AbstractArray{Float64, 2}

    # Climatology states
    Ts_clim_relax_time :: Union{Float64, Nothing}
    Ss_clim_relax_time :: Union{Float64, Nothing}

    Ts_clim  :: Union{AbstractArray{Float64, 3}, Nothing}
    Ss_clim  :: Union{AbstractArray{Float64, 3}, Nothing}



    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    # Derived quantities
    zs       :: AbstractArray{Float64, 1}
    N_ocs    :: Integer           # Number of columns
    
    function OceanColumnCollection(;
        Nx       :: Integer,
        Ny       :: Integer,
        hs       :: AbstractArray{Float64, 1},
        T_ML     :: Union{AbstractArray{Float64, 2}, Float64},
        S_ML     :: Union{AbstractArray{Float64, 2}, Float64},
        Ts       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        Ss       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        mask     :: Union{AbstractArray{Float64, 2}, Nothing},
        topo     :: Union{AbstractArray{Float64, 2}, Nothing},

        h_ML_min :: Union{AbstractArray{Float64, 2}, Float64},
        h_ML_max :: Union{AbstractArray{Float64, 2}, Float64},
        Ts_clim_relax_time :: Union{Float64, Nothing},
        Ss_clim_relax_time :: Union{Float64, Nothing},
        Ts_clim  :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Nothing},
        Ss_clim  :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Nothing},
    )

        if length(hs) != 2 && all(hs .> 0)
            throw(ErrorException("EntSOM needs more than 2 layers. Also all values in `hs` should be positive."))
        end

        Nz = length(hs) 

        _hs = copy(hs)
        _zs = zeros(Float64, Nz+1)
        _zs[1] = 0.0

        for i = 2:length(_zs)
            _zs[i] = _zs[i-1] - _hs[i-1]
        end

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
        

        # ===== [BEGIN] Column information =====

        _bs       = SharedArray{Float64}(Nx, Ny, 2)
        _Ts       = SharedArray{Float64}(Nx, Ny, 2)
        _Ss       = SharedArray{Float64}(Nx, Ny, 2)
        _qflx2atm = SharedArray{Float64}(Nx, Ny)

        _qflx2atm .= 0.0


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

        # ===== [BEGIN] Climatology =====

        # TODO: Need to detect whether all
        #       climatology data points are
        #       valid or not.

        if Ts_clim == nothing

            _Ts_clim = nothing

        else
            
            _Ts_clim = SharedArray{Float64}(Nx, Ny, Nz_bone)
            
            if typeof(Ts_clim) <: AbstractArray{Float64, 3}

                _Ts_clim[:, :, :] = Ts_clim

            elseif typeof(Ts_clim) <: AbstractArray{Float64, 1}

                for i=1:Nx, j=1:Ny
                    _Ts_clim[i, j, :] = Ts_clim
                end

            end

        end 


        if Ss_clim == nothing
            
            _Ss_clim = nothing

        else
            
            _Ss_clim = SharedArray{Float64}(Nx, Ny, Nz_bone)
            
            if typeof(Ss_clim) <: AbstractArray{Float64, 3}

                _Ss_clim[:, :, :] = Ss_clim

            elseif typeof(Ss_clim) <: AbstractArray{Float64, 1}

                for i=1:Nx, j=1:Ny
                    _Ss_clim[i, j, :] = Ss_clim
                end

            end

        end 

        # ===== [END] Climatology =====

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
