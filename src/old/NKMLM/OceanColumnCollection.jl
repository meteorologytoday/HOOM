mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    Nz_bone  :: Integer           # Number of layers  in k direction
    
    zs_bone  :: AbstractArray{Float64, 1} # Unmasked zs bone
    topo     :: AbstractArray{Float64, 2} # Depth of the topography. Negative value if it is underwater
    zs       :: AbstractArray{Float64, 3} # Actuall zs coordinate masked by topo
    Nz       :: AbstractArray{  Int64, 2} # Number of layers that is active

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
    FLDO     :: AbstractArray{  Int64, 2}
    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes

    h_ML_min :: AbstractArray{Float64, 2}
    h_ML_max :: AbstractArray{Float64, 2}
    we_max   :: Float64

    # Climatology states
    Ts_clim_relax_time :: Union{Float64, Nothing}
    Ss_clim_relax_time :: Union{Float64, Nothing}

    Ts_clim  :: Union{AbstractArray{Float64, 3}, Nothing}
    Ss_clim  :: Union{AbstractArray{Float64, 3}, Nothing}

    # Derived quantities
    N_ocs  :: Integer           # Number of columns
    hs     :: AbstractArray{Float64, 3} # Thickness of layers
    Δzs    :: AbstractArray{Float64, 3} # Δz between layers

    # 1D Views to make clean code
    zs_vw :: Any 
    bs_vw :: Any
    Ts_vw :: Any
    Ss_vw :: Any
    Ts_clim_vw :: Any
    Ss_clim_vw :: Any

    function OceanColumnCollection(;
        Nx       :: Integer,
        Ny       :: Integer,
        zs_bone  :: AbstractArray{Float64, 1},
        Ts       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        Ss       :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Float64},
        K_T      :: Float64,
        K_S      :: Float64,
        T_ML     :: Union{AbstractArray{Float64, 2}, Float64},
        S_ML     :: Union{AbstractArray{Float64, 2}, Float64},
        h_ML     :: Union{AbstractArray{Float64, 2}, Float64, Nothing},
        h_ML_min :: Union{AbstractArray{Float64, 2}, Float64},
        h_ML_max :: Union{AbstractArray{Float64, 2}, Float64},
        we_max   :: Float64,
        Ts_clim_relax_time :: Union{Float64, Nothing},
        Ss_clim_relax_time :: Union{Float64, Nothing},
        Ts_clim  :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Nothing},
        Ss_clim  :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1}, Nothing},
        mask     :: Union{AbstractArray{Float64, 2}, Nothing},
        topo     :: Union{AbstractArray{Float64, 2}, Nothing},
    )

        # ===== [BEGIN] topo, mask, h_ML_min, h_ML_max =====
        # Min/max of ML is tricky because it cannot be
        # deeper than the bottom boundary
        # Also, in real data topo can be 0 and not masked out
       
        _topo = SharedArray{Float64}(Nx, Ny)
        _h_ML_min = SharedArray{Float64}(Nx, Ny)
        _h_ML_max = SharedArray{Float64}(Nx, Ny)
        _mask = SharedArray{Float64}(Nx, Ny)

        if topo == nothing
            _topo .= zs_bone[end]
        else
            _topo[:, :] = topo
        end
       
        # mask =>   lnd = 0, ocn = 1
        if mask == nothing
            _mask .+= 1.0
        else
            _mask[:, :] = mask
        end

        mask_idx = (_mask .== 1.0)


        if typeof(h_ML_min) <: AbstractArray{Float64, 2}
            _h_ML_min[:, :] = h_ML_min
        elseif typeof(h_ML_min) <: Float64
            _h_ML_min .= h_ML_min
        end

        if typeof(h_ML_max) <: AbstractArray{Float64, 2}
            _h_ML_max[:, :] = h_ML_max
        elseif typeof(h_ML_max) <: Float64
            _h_ML_max .= h_ML_max
        end

        # Detect and fix h_ML_{max,min}
        coord_max = -zs_bone[end]
        for i=1:Nx, j=1:Ny

            if _mask[i, j] == 0
                _h_ML_max[i, j] = NaN
                _h_ML_min[i, j] = NaN
                continue
            end

            hmax = _h_ML_max[i, j]
            hmin = _h_ML_min[i, j]
            hbot = - _topo[i, j]

            if hbot < 0
                throw(ErrorException(format("Topography is negative at idx ({:d}, {:d})", i, j)))
            elseif hbot == 0
                #throw(ErrorException(format("Topography is zero at idx ({:d}, {:d})", i, j)))
                println(format("Topography is zero at idx ({:d}, {:d})", i, j))
            end

            if hmax < hmin
                throw(ErrorException(format("h_ML_max must ≥ h_ML_min. Problem happens at idx ({:d}, {:d})", i, j)))
            end

            if coord_max < hmin
                println(format("Point ({},{}) got h_min {:.2f} which is larger than coord_max {}. Tune h_ML_min to match it.", i, j, hmin, coord_max))
                hmin = coord_max
            end

            if coord_max < hmax
                println(format("Point ({},{}) got h_max {:.2f} which is larger than coord_max {}. Tune h_ML_max to match it.", i, j, hmax, coord_max))
                hmax = coord_max
            end
 


            if hmin > hbot
                #println(mask[i,j]) 
                println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_min {}. Tune h_ML_min/max to depth.", i, j, hbot, hmin))
                hbot = hmin
            end

            if hmax > hbot
                println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_max {}. Tune the h_ML_max to depth.", i, j, hbot, hmax))
                hmax = hbot
            end

            _h_ML_min[i, j] = hmin
            _h_ML_max[i, j] = hmax
            _topo[i, j]     = -hbot

        end
      
        
        # ===== [END] topo, mask, h_ML_min, h_ML_max =====

        # ===== [BEGIN] z coordinate =====
        zs_bone = copy(zs_bone)
        Nz_bone = length(zs_bone) - 1

        Nz   = SharedArray{Int64}(Nx, Ny)
        zs   = SharedArray{Float64}(Nx, Ny, Nz_bone + 1)
        hs   = SharedArray{Float64}(Nx, Ny, Nz_bone    )
        Δzs  = SharedArray{Float64}(Nx, Ny, Nz_bone - 1)

        zs  .= NaN
        Nz  .= 0
        hs  .= NaN
        Δzs .= NaN

        for i=1:Nx, j=1:Ny

            if _mask[i, j] == 0
                continue
            end

            # Determine Nz

            # Default is that topo is deeper than
            # the bottom of zs_bone
            _Nz = Nz_bone
            for k=2:length(zs_bone)
                if zs_bone[k] <= _topo[i, j]
                    _Nz = k-1
                    #println(format("This topo gets: zs_bone[{:d}] = {:f}, _topo[{:d},{:d}]={:f}", k, zs_bone[k], i, j, _topo[i,j]))
                    break
                end
            end

            Nz[i, j] = _Nz

            # Construct vertical coordinate
            zs[i, j, 1:_Nz] = zs_bone[1:_Nz]

            zs[i, j, _Nz+1] = max(_topo[i, j], zs_bone[_Nz+1])

            # Construct thickness of each layer
            hs[ i, j, 1:_Nz]  = zs[i, j, 1:_Nz] - zs[i, j, 2:_Nz+1]
            Δzs[i, j, 1:_Nz-1] = (hs[i, j, 1:_Nz-1] + hs[i, j, 2:_Nz]) / 2.0
            
        end
        
        # ===== [END] z coordinate =====

        # ===== [BEGIN] Column information =====

        _b_ML     = SharedArray{Float64}(Nx, Ny)
        _T_ML     = SharedArray{Float64}(Nx, Ny)
        _S_ML     = SharedArray{Float64}(Nx, Ny)
        _h_ML     = SharedArray{Float64}(Nx, Ny)

        _bs       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _Ts       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _Ss       = SharedArray{Float64}(Nx, Ny, Nz_bone)
        _FLDO     = SharedArray{Int64}(Nx, Ny)
        qflx2atm  = SharedArray{Float64}(Nx, Ny)


        if typeof(h_ML) <: AbstractArray{Float64, 2}
            _h_ML[:, :] = h_ML
        elseif typeof(h_ML) <: Float64
            _h_ML .= h_ML
        elseif h_ML == nothing
            _h_ML .= h_ML_min
        end

        if typeof(T_ML) <: AbstractArray{Float64, 2}
            _T_ML[:, :] = T_ML
        elseif typeof(T_ML) <: Float64
            _T_ML .= T_ML
        end

        if typeof(S_ML) <: AbstractArray{Float64, 2}
            _S_ML[:, :] = S_ML
        elseif typeof(S_ML) <: Float64
            _S_ML .= S_ML
        end

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

        # ===== [BEGIN] Construct Views =====
        zs_vw = Array{SubArray}(undef, Nx, Ny)
        bs_vw = Array{SubArray}(undef, Nx, Ny)
        Ts_vw = Array{SubArray}(undef, Nx, Ny)
        Ss_vw = Array{SubArray}(undef, Nx, Ny)

        for i=1:Nx, j=1:Ny
            zs_vw[i, j]      = view(zs,  i, j, :)
            bs_vw[i, j]      = view(_bs, i, j, :)
            Ts_vw[i, j]      = view(_Ts, i, j, :)
            Ss_vw[i, j]      = view(_Ss, i, j, :)
        end

        Ts_clim_vw = nothing
        Ss_clim_vw = nothing

        if Ts_clim != nothing
            Ts_clim_vw = Array{SubArray}(undef, Nx, Ny)
            for i=1:Nx, j=1:Ny
                Ts_clim_vw[i, j] = view(_Ts_clim, i, j, :)
            end
        end
 
        if Ss_clim != nothing
            Ss_clim_vw = Array{SubArray}(undef, Nx, Ny)
            for i=1:Nx, j=1:Ny
                Ss_clim_vw[i, j] = view(_Ss_clim, i, j, :)
            end
        end
     
        # ===== [END] Construct Views =====

        # ===== [BEGIN] Mask out data =====

        mask3 = zeros(Int64, Nx, Ny, Nz_bone)
        mask3 .= 1

        # Clean up all variables
        for i=1:Nx, j=1:Ny
            mask3[i, j, Nz[i, j] + 1:end] .= 0 
        end

        println("sum of mask3: ", sum(mask3))

        mask3_lnd_idx = (mask3  .== 0)
        mask2_lnd_idx = (_mask  .== 0)
        
        for v in [_bs, _Ts, _Ss, _Ts_clim, _Ss_clim]
            if v == nothing
                continue
            end

            v[mask3_lnd_idx] .= NaN
        end 

        for v in [_b_ML, _T_ML, _S_ML, _h_ML, _h_ML_min, _h_ML_max]
            v[mask2_lnd_idx] .= NaN
        end 


        # ===== [END] Mask out data

        # ===== [BEGIN] check integrity =====

        # Check topography, h_ML_min/max and zs
        for i=1:Nx, j=1:Ny
            if _mask[i, j] == 0
                continue
            end

            if ! (-_h_ML_min[i, j] >= - _h_ML_max[i, j] >= zs[i, j, Nz[i, j] + 1] >= _topo[i, j])
                println("idx: (", i, ", ", j, ")")
                println("h_ML_min: ", _h_ML_min[i, j])
                println("h_ML_max: ", _h_ML_max[i, j])
                println("z_deepest: ", zs[i, j, Nz[i, j] + 1])
                println("topo: ", _topo[i, j])
                ErrorException("Relative relation is wrong") |> throw
            end
        end

        
        # Check if there is any hole in climatology 
        
        mask3_idx = (mask3 .== 1)
        valid_grids = sum(mask3_idx)
        total_data  = Nx * Ny * Nz_bone

        println("Total  data count: ", total_data)
        println("Valid  data count: ", valid_grids)
        println("Masked data count: ", total_data - valid_grids)
        
        if sum(isfinite.(_Ss)) != valid_grids
            throw(ErrorException("Salinity data has holes"))
        end
 
        if sum(isfinite.(_Ts)) != valid_grids
            throw(ErrorException("Temperature data has holes"))
        end
 
        if _Ts_clim != nothing && sum(isfinite.(_Ts_clim)) != valid_grids
            throw(ErrorException("Temperature climatology has holes"))
        end
 
        if _Ss_clim != nothing && sum(isfinite.(_Ss_clim)) != valid_grids
            throw(ErrorException("Salinity climatology has holes"))
        end


        # Check if h_ML_min h_ML_max is negative
        if any(_h_ML_min[mask_idx] .<= 0)
            throw(ErrorException("h_ML_min should always be positive (cannot be zero or negative)"))
        end

        if any(_h_ML_max[mask_idx] .< 0)
            throw(ErrorException("h_ML_max should always be non-negative"))
        end 
        
        # ===== [END] check integrity =====


        # ===== [BEGIN] Mask out values below topo =====
        # ===== [END] Mask out values below topo =====


        occ = new(
            Nx, Ny, Nz_bone,
            zs_bone, _topo, zs, Nz,
            K_T, K_S,
            _mask, mask_idx,
            _b_ML, _T_ML, _S_ML, _h_ML,
            _bs,   _Ts,   _Ss,
            _FLDO, qflx2atm,
            _h_ML_min, _h_ML_max, we_max,
            Ts_clim_relax_time, Ss_clim_relax_time,
            _Ts_clim, _Ss_clim,
            Nx * Ny, hs, Δzs,
            zs_vw, bs_vw, Ts_vw, Ss_vw, Ts_clim_vw, Ss_clim_vw,
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



function makeBasicOceanColumnCollection(;
    Nx      :: Integer,
    Ny      :: Integer,
    zs_bone :: AbstractArray{Float64, 1},
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
    Ts_clim :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1},  Nothing} = nothing,
    Ss_clim :: Union{AbstractArray{Float64, 3}, AbstractArray{Float64, 1},  Nothing} = nothing,
    Ts_clim_relax_time :: Union{Float64, Nothing} = nothing,
    Ss_clim_relax_time :: Union{Float64, Nothing} = nothing,

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
        Nx                 = Nx,
        Ny                 = Ny,
        zs_bone            = zs_bone,
        Ts                 = Ts,
        Ss                 = Ss,
        K_T                = K_T,
        K_S                = K_S,
        T_ML               = T_ML,
        S_ML               = S_ML,
        h_ML               = h_ML,
        h_ML_min           = h_ML_min,        
        h_ML_max           = h_ML_max,
        we_max             = we_max,
        Ts_clim_relax_time = Ts_clim_relax_time,
        Ss_clim_relax_time = Ss_clim_relax_time,
        Ts_clim            = Ts_clim,
        Ss_clim            = Ss_clim,
        mask               = mask,
        topo               = topo,
    )
end
