mutable struct Core

    amo :: Union{AdvancedMatrixOperators, Nothing}

    wksp :: Workspace

    function Core(
        ev :: Env;
    )
    
    println("Test if slab AMO is there")
    @time amo = AdvancedMatrixOperators(;
        gd = ev.gd_slab,
        mask_T = env.mask_T,
    )


        # Build Advection Matrix

        # Build Radiation Matrix
        swflx_factor_W = (1.0 - ev.R) * exp.(ev.gd.z_W / ζ2)
        lwflx_factor_W =        ev.R  * exp.(ev.gd.z_W / ζ1)

        # Bottom absorbs everything
        swflx_factor_W[end, :, :] .= 0.0
        lwflx_factor_W[end, :, :] .= 0.0

        function build!(id_mtx, idx)
            local result
            rows = size(id_mtx)[1]
            
            # using transpose speeds up by 100 times 
            tp = transpose(id_mtx) |> sparse
            result = transpose(tp[:, view(idx, :)]) |> sparse
            dropzeros!(result)
            idx .= 0 # clean so that debug is easir when some girds are not assigned
            return result
        end

        # Build the matrix to broadcast sW to W grid
        num_sW = reshape( collect(1:amo_slab.bmo.W_pts), amo_slab.bmo.W_dim...)
        W = repeat(num_sW, outer=(ev.Nz+1, 1, 1))
        W_broadcast_sW = build!(amo_slab.bmo.W_I_W, W)

        swflx_conv = - amo.T_DIVz_W * spdiagm(0 => view(swflx_factor_W, :)) * W_broadcast_sW 
        lwflx_conv = - amo.T_DIVz_W * spdiagm(0 => view(lwflx_factor_W, :)) * W_broadcast_sW 


        # ===== [BEGIN] topo, mask, h_ML_min, h_ML_max =====
        # Min/max of ML is tricky because it cannot be
        # deeper than the bottom boundary
        # Also, in real data topo can be 0 and not masked out
       
        _topo        = allocate(datakind, Float64, Nx, Ny)
        _h_ML_min    = allocate(datakind, Float64, Nx, Ny)
        _h_ML_max    = allocate(datakind, Float64, Nx, Ny)



        if topo == nothing
            _topo .= zs_bone[end]
        else
            _topo[:, :] = topo
        end
       
        # Arrage like (2, cnt) instead of (cnt, 2) to
        # enhance speed through memory cache
        valid_idx = allocate(datakind, Int64, 2, sum(mask_idx))
        
        let k = 1
            for idx in CartesianIndices((Nx, Ny))
                if _mask[idx] == 1.0
                    valid_idx[1, k] = idx[1]
                    valid_idx[2, k] = idx[2]

                    k += 1
                end
            end

            if k != size(valid_idx)[2] + 1
                throw(ErrorException("Initialization error making `valid_idx`"))
            end
        end


        #println("################ h_ML_min: ", h_ML_min)


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
                verbose && println(format("Point ({},{}) got h_min {:.2f} which is larger than coord_max {}. Tune h_ML_min to match it.", i, j, hmin, coord_max))
                hmin = coord_max
            end

            if coord_max < hmax
                verbose && println(format("Point ({},{}) got h_max {:.2f} which is larger than coord_max {}. Tune h_ML_max to match it.", i, j, hmax, coord_max))
                hmax = coord_max
            end
 


            if hmin > hbot
                verbose && println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_min {}. Tune h_ML_min/max to depth.", i, j, hbot, hmin))
                hbot = hmin
            end

            if hmax > hbot
                verbose && println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_max {}. Tune the h_ML_max to depth.", i, j, hbot, hmax))
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

        Nz   = allocate(datakind, Int64, Nx, Ny)
        zs   = allocate(datakind, Float64, Nz_bone + 1, Nx, Ny)
        hs   = allocate(datakind, Float64, Nz_bone    , Nx, Ny)
        Δzs  = allocate(datakind, Float64, Nz_bone + 1, Nx, Ny)

        zs  .= 1.0
        Nz  .= 0
        hs  .= 1.0
        Δzs .= 1.0

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
            zs[1:_Nz, i, j] = zs_bone[1:_Nz]

            zs[_Nz+1, i, j] = max(_topo[i, j], zs_bone[_Nz+1])

            # Construct thickness of each layer
            hs[ 1:_Nz,   i, j] = zs[1:_Nz, i, j] - zs[2:_Nz+1, i, j]
            Δzs[2:_Nz, i, j] = (hs[1:_Nz-1, i, j] + hs[2:_Nz, i, j]) / 2.0
            Δzs[1, i, j] = Δzs[2, i, j]
            Δzs[end, i, j] = Δzs[end-1, i, j]
            
           
        end
        
        # ===== [END] z coordinate =====

        # ===== [BEGIN] Column information =====

        _b_ML     = allocate(datakind, Float64, Nx, Ny)
        _T_ML     = allocate(datakind, Float64, Nx, Ny)
        _S_ML     = allocate(datakind, Float64, Nx, Ny)
        _ΔT       = allocate(datakind, Float64, Nx, Ny)
        _ΔS       = allocate(datakind, Float64, Nx, Ny)
        _dΔTdt    = allocate(datakind, Float64, Nx, Ny)
        _dΔSdt    = allocate(datakind, Float64, Nx, Ny)
        _h_ML     = allocate(datakind, Float64, Nx, Ny)
        _h_MO     = allocate(datakind, Float64, Nx, Ny)
        _fric_u   = allocate(datakind, Float64, Nx, Ny)
        _dTdt_ent    = allocate(datakind, Float64, Nx, Ny)
        _dSdt_ent    = allocate(datakind, Float64, Nx, Ny)

        _TSAS_clim   = allocate(datakind, Float64, Nx, Ny)
        _SSAS_clim   = allocate(datakind, Float64, Nx, Ny)
        _TFLUX_bot       = allocate(datakind, Float64, Nx, Ny)
        _SFLUX_bot       = allocate(datakind, Float64, Nx, Ny)
        _SFLUX_top       = allocate(datakind, Float64, Nx, Ny)
        _TFLUX_DIV_implied      = allocate(datakind, Float64, Nx, Ny)
        _SFLUX_DIV_implied      = allocate(datakind, Float64, Nx, Ny)

        _bs       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _Ts       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _Ss       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _FLDO     = allocate(datakind, Int64, Nx, Ny)
        _FLDO_ratio_top = allocate(datakind, Float64, Nx, Ny)
        _FLDO_ratio_bot = allocate(datakind, Float64, Nx, Ny)
        _qflx2atm  = allocate(datakind, Float64, Nx, Ny)
        _qflx2atm_pos = allocate(datakind, Float64, Nx, Ny)
        _qflx2atm_neg = allocate(datakind, Float64, Nx, Ny)
        _TEMP         = allocate(datakind, Float64, Nx, Ny)
        _dTEMPdt      = allocate(datakind, Float64, Nx, Ny)
        _SALT      = allocate(datakind, Float64, Nx, Ny)
        _dSALTdt   = allocate(datakind, Float64, Nx, Ny)
        _qflx_T_correction = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _qflx_S_correction = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _seaice_nudge_energy = allocate(datakind, Float64, Nx, Ny)

        if typeof(h_ML) <: AbstractArray{Float64, 2}
            _h_ML[:, :] = h_ML
        elseif typeof(h_ML) <: Float64
            _h_ML .= h_ML
        elseif h_ML == nothing
            _h_ML .= h_ML_min
        end

        # Need to constraint h_ML
        for i=1:Nx, j=1:Ny

            (_mask[i, j] == 0.0) && continue
            
            _h_ML[i, j] = boundMLD(_h_ML[i, j]; h_ML_max=_h_ML_max[i, j], h_ML_min=_h_ML_min[i, j])

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
            _Ts[:, :, :] = toZXY(Ts, arrange)
        elseif typeof(Ts) <: AbstractArray{Float64, 1}
            for i=1:Nx, j=1:Ny
                _Ts[:, i, j] = Ts
            end
        elseif typeof(Ts) <: Float64 
            _Ts .= Ts
        end

        if typeof(Ss) <: AbstractArray{Float64, 3}
            _Ss[:, :, :] = toZXY(Ss, arrange)
        elseif typeof(Ss) <: AbstractArray{Float64, 1}
            for i=1:Nx, j=1:Ny
                _Ss[:, i, j] = Ss
            end
        elseif typeof(Ss) <: Float64 
            _Ss .= Ss
        end

        _Ts_mixed = copy(_Ts)
        _Ss_mixed = copy(_Ts)

        # ===== [END] Column information =====

        # ===== [BEGIN] Radiation =====

        _rad_decay_coes  = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _rad_absorp_coes = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        for i=1:Nx, j=1:Ny

            if _mask[i, j] == 0.0
                continue
            end

            for k=1:Nz[i, j]
                _rad_decay_coes[k, i, j]  = exp(zs[k, i, j] / ζ)         # From surface to top of the layer
                _rad_absorp_coes[k, i, j] = 1.0 - exp(- hs[k, i, j] / ζ)
            end

            # Since we assume the bottome of ocean absorbs anything
            _rad_absorp_coes[Nz[i, j], i, j] = 1.0
        end


        # ===== [END] Radiation =====




        # ===== [BEGIN] fs and ϵs =====

        _fs       = allocate(datakind, Float64, Nx, Ny)
        _ϵs       = allocate(datakind, Float64, Nx, Ny)

        if typeof(fs) <: AbstractArray{Float64, 2}
            _fs[:, :] = fs
        elseif typeof(fs) <: Float64 
            _fs .= fs
        elseif fs == nothing
           _fs[:, :] = 2 * Ωe * sin.(gridinfo.c_lat)
        end

        if typeof(ϵs) <: AbstractArray{Float64, 2}
            _ϵs[:, :] = ϵs
        elseif typeof(ϵs) <: Float64 
            _ϵs .= ϵs
        end

        # ===== [END] fs and ϵs =====

        # ===== [BEGIN] Advection variables =====

        _τx      = allocate(datakind, Float64, Nx, Ny)
        _τy      = allocate(datakind, Float64, Nx, Ny)
        _u       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _v       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _w       = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        _u_bnd   = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _v_bnd   = allocate(datakind, Float64, Nz_bone, Nx, Ny+1)
        _w_bnd   = allocate(datakind, Float64, Nz_bone+1, Nx, Ny)

        _GRAD_bnd_x   = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _GRAD_bnd_y   = allocate(datakind, Float64, Nz_bone, Nx, Ny+1)
        _GRAD_bnd_z   = allocate(datakind, Float64, Nz_bone+1, Nx, Ny)

        _CURV_x       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _CURV_y       = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _CURV_z       = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        _TFLUX_DEN_x  = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _TFLUX_DEN_y  = allocate(datakind, Float64, Nz_bone, Nx, Ny+1)
        _TFLUX_DEN_z  = allocate(datakind, Float64, Nz_bone+1, Nx, Ny)

        _SFLUX_DEN_x  = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _SFLUX_DEN_y  = allocate(datakind, Float64, Nz_bone, Nx, Ny+1)
        _SFLUX_DEN_z  = allocate(datakind, Float64, Nz_bone+1, Nx, Ny)


        _div     = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        _TFLUX_CONV   = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _TFLUX_CONV_h = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _SFLUX_CONV   = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _SFLUX_CONV_h = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        _∇∇T     = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        _∇∇S     = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        # ===== [END] Advection variables =====

        # ===== [BEGIN] Climatology =====

        _Ts_clim = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        if Ts_clim == nothing

        else
            
            if typeof(Ts_clim) <: AbstractArray{Float64, 3}

                _Ts_clim[:, :, :] = toZXY(Ts_clim, arrange)

            elseif typeof(Ts_clim) <: AbstractArray{Float64, 1}

                for i=1:Nx, j=1:Ny
                    _Ts_clim[:, i, j] = Ts_clim
                end

            end

        end 


        _Ss_clim = allocate(datakind, Float64, Nz_bone, Nx, Ny)
        if Ss_clim == nothing
            
        else
            
            if typeof(Ss_clim) <: AbstractArray{Float64, 3}

                _Ss_clim[:, :, :] = toZXY(Ss_clim, arrange)

            elseif typeof(Ss_clim) <: AbstractArray{Float64, 1}

                for i=1:Nx, j=1:Ny
                    _Ss_clim[:, i, j] = Ss_clim
                end

            end

        end 

        # ===== [END] Climatology =====


        # ===== [BEGIN] Mask out data =====

        _mask3        = allocate(datakind, Float64, Nz_bone, Nx, Ny)

        _mask3 .= 1.0
        # Clean up all variables
        for i=1:Nx, j=1:Ny
            _mask3[Nz[i, j] + 1:end, i, j] .= 0.0
        end

        #println("sum of _mask3: ", sum(_mask3))

        mask3_lnd_idx = (_mask3  .== 0.0)
        mask2_lnd_idx = (_mask  .== 0.0)

        for v in [_bs, _Ts, _Ss, _Ts_clim, _Ss_clim]
            if v == nothing
                continue
            end

            v[mask3_lnd_idx] .= -999
        end 

        for v in [_b_ML, _T_ML, _S_ML, _h_ML, _h_ML_min, _h_ML_max]
            v[mask2_lnd_idx] .= -999
        end 

        # ===== [END] Mask out data

        # ==== [BEGIN] Determine nomotionmask3 ====
        _nomotionmask3 = copy(_mask3)
        for i=1:Nx, j=1:Ny

            if _topo[i, j] > -300.0
                _nomotionmask3[:, i, j] .= 0.0
                continue
            end

            if _mask[i, j] == 1.0
                _nomotionmask3[ Nz[i, j], i, j ] = 0.0
            end
        
        end
        # ==== [END] Determine nomotionmask3 ====


        # ===== [BEGIN] check integrity =====

        # Check topography, h_ML_min/max and zs
        for i=1:Nx, j=1:Ny
            if _mask[i, j] == 0
                continue
            end

            if ! (-_h_ML_min[i, j] >= - _h_ML_max[i, j] >= zs[Nz[i, j] + 1, i, j] >= _topo[i, j])
                println("idx: (", i, ", ", j, ")")
                println("h_ML_min: ", _h_ML_min[i, j])
                println("h_ML_max: ", _h_ML_max[i, j])
                println("z_deepest: ", zs[Nz[i, j] + 1, i, j])
                println("topo: ", _topo[i, j])
                ErrorException("Relative relation is wrong") |> throw
            end
        end

        
        # Check if there is any hole in climatology 
        
        mask3_idx = (_mask3 .== 1)
        valid_grids = sum(mask3_idx)
        total_data  = Nx * Ny * Nz_bone

        #println("Total  data count: ", total_data)
        #println("Valid  data count: ", valid_grids)
        #println("Masked data count: ", total_data - valid_grids)

        #=        
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
        =#

        # Check if h_ML_min h_ML_max is negative
        if any(_h_ML_min[mask_idx] .<= 0)
            throw(ErrorException("h_ML_min should always be positive (cannot be zero or negative)"))
        end

        if any(_h_ML_max[mask_idx] .< 0)
            throw(ErrorException("h_ML_max should always be non-negative"))
        end 
        
        # ===== [END] check integrity =====

        cols = nothing
        lays = nothing

#        if id != 0

            # ===== [BEGIN] Construct Views of Lays =====
            lays = ((
                hs      = Array{SubArray}(undef, Nz_bone),
                bs      = Array{SubArray}(undef, Nz_bone),
                Ts      = Array{SubArray}(undef, Nz_bone),
                Ss      = Array{SubArray}(undef, Nz_bone),
                u       = Array{SubArray}(undef, Nz_bone),
                v       = Array{SubArray}(undef, Nz_bone),
                div     = Array{SubArray}(undef, Nz_bone),
                ∇∇T     = Array{SubArray}(undef, Nz_bone),
                ∇∇S     = Array{SubArray}(undef, Nz_bone),
                mask3   = Array{SubArray}(undef, Nz_bone),
            ))
     
            
            for k=1:Nz_bone
                lays.hs[k]      = view(hs, k, :, :)
                lays.bs[k]      = view(_bs, k, :, :)
                lays.Ts[k]      = view(_Ts, k, :, :)
                lays.Ss[k]      = view(_Ss, k, :, :)
                lays.u[k]       = view(_u, k, :, :)
                lays.v[k]       = view(_v, k, :, :)
                lays.div[k]     = view(_div, k, :, :)
                lays.∇∇T[k]     = view(_∇∇T, k, :, :)
                lays.∇∇S[k]     = view(_∇∇S, k, :, :)
                lays.mask3[k]   = view(_mask3, k, :, :)
            end

            # ===== [END] Construct Views of Lays =====

            # ===== [BEGIN] Construct Views of Cols =====
            cols = ((
                zs  = Array{SubArray}(undef, Nx, Ny),
                Δzs = Array{SubArray}(undef, Nx, Ny),
                hs = Array{SubArray}(undef, Nx, Ny),
                w  = Array{SubArray}(undef, Nx, Ny),
                w_bnd = Array{SubArray}(undef, Nx, Ny),
                bs = Array{SubArray}(undef, Nx, Ny),
                Ts = Array{SubArray}(undef, Nx, Ny),
                Ss = Array{SubArray}(undef, Nx, Ny),
                Ts_mixed = Array{SubArray}(undef, Nx, Ny),
                Ss_mixed = Array{SubArray}(undef, Nx, Ny),
                rad_decay_coes  = Array{SubArray}(undef, Nx, Ny),
                rad_absorp_coes = Array{SubArray}(undef, Nx, Ny),
                Ts_clim = (Ts_clim == nothing) ? nothing : Array{SubArray}(undef, Nx, Ny),
                Ss_clim = (Ss_clim == nothing) ? nothing : Array{SubArray}(undef, Nx, Ny),
                qflx_T_correction = Array{SubArray}(undef, Nx, Ny),
                qflx_S_correction = Array{SubArray}(undef, Nx, Ny),
            ))
            
            for i=1:Nx, j=1:Ny
                cols.zs[i, j]              = view(zs,  :, i, j)
                cols.Δzs[i, j]             = view(Δzs, :, i, j)
                cols.hs[i, j]              = view(hs,  :, i, j)
                cols.w[i, j]               = view(_w,  :, i, j)
                cols.w_bnd[i, j]           = view(_w_bnd,  :, i, j)
                cols.bs[i, j]              = view(_bs, :, i, j)
                cols.Ts[i, j]              = view(_Ts, :, i, j)
                cols.Ss[i, j]              = view(_Ss, :, i, j)
                cols.Ts_mixed[i, j]        = view(_Ts_mixed, :, i, j)
                cols.Ss_mixed[i, j]        = view(_Ss_mixed, :, i, j)
                cols.rad_decay_coes[i, j]  = view(_rad_decay_coes,  :, i, j)
                cols.rad_absorp_coes[i, j] = view(_rad_absorp_coes, :, i, j)
                cols.qflx_T_correction[i, j] = view(_qflx_T_correction, :, i, j)
                cols.qflx_S_correction[i, j] = view(_qflx_S_correction, :, i, j)
            end

            if Ts_clim != nothing
                for i=1:Nx, j=1:Ny
                    cols.Ts_clim[i, j] = view(_Ts_clim, :, i, j)
                end
            end
 
            if Ss_clim != nothing
                for i=1:Nx, j=1:Ny
                    cols.Ss_clim[i, j] = view(_Ss_clim, :, i, j)
                end
            end

#        end
        # ===== [END] Construct Views of Cols =====

        # ===== [BEGIN] Making speed-up matrix

        if true || id != 0
            ASUM = AdvectionSpeedUpMatrix(;
                gi    = gridinfo,
                Nz    = Nz_bone,
                Nz_av = Nz,
                mask3 = _mask3,
                Δz_W  = Δzs,
                Δz_T  = hs,
                nomotionmask3 = _nomotionmask3,
            )
        else
            ASUM = nothing
        end
        # ===== [END] Making speed-up matrix

        deep_ocn_correction_start_layer = getLayerFromDepth(;
            zs = zs_bone,
            z  = -1.0,
            Nz = length(zs_bone) - 1
        )

        diag = Dict(
            :total_heat => zeros(Float64, 1),
            :total_salt => zeros(Float64, 1),
            :total_heat_budget_residue => zeros(Float64, 1),
            :total_salt_budget_residue => zeros(Float64, 1),
        )

        # ===== [BEGIN] Q-flux finding forcing file =====
        qflx_finding_cdm = nothing
        if id == 0
           
            # do nothing 
            
        else
            if sub_yrng == nothing
                thorw(ErrorException("Init worker ocean, sub_yrng must be provided."))
            end

            if qflx_finding_forcing_file != ""

                println("Load qflx file: ", qflx_finding_forcing_file)
                qflx_finding_cdm = CyclicDataManager(;
                    filename     = qflx_finding_forcing_file,
                    varname_time = "time", 
                    varnames     = ["TEMP", "SALT"],
                    beg_time     = 0.0,
                    cyc_time     = 365.0,
                    spatial_rng  = (:, sub_yrng, :),  # in the shape of the file
                    xyz2zxy      = true,
                )
            end
               
        end


        # ===== [END] Q-flux finding forcing file =====


        ocn = new(
            id,
            gridinfo,
            gridinfo_file,
            mi,
            Nx, Ny, Nz_bone,
            zs_bone, _topo, zs, Nz,
            K_v, Dh_T, Dv_T, Dh_S, Dv_S,
            _fs, _ϵs,
            _mask3, _nomotionmask3,
            _mask, mask_idx, valid_idx,
            _b_ML, _T_ML, _S_ML, _ΔT, _ΔS, _dΔTdt, _dΔSdt,
            _h_ML, _h_MO, _fric_u, _dTdt_ent, _dSdt_ent,
            _TSAS_clim, _SSAS_clim,
            _TFLUX_bot, _SFLUX_bot, _SFLUX_top,
            _TFLUX_DIV_implied, _SFLUX_DIV_implied,
            _bs,   _Ts,   _Ss,
            _Ts_mixed,   _Ss_mixed,
            _FLDO, _FLDO_ratio_top, _FLDO_ratio_bot,
            _qflx2atm, _qflx2atm_pos, _qflx2atm_neg,
            _TEMP, _dTEMPdt, _SALT, _dSALTdt,
            _qflx_T_correction, _qflx_S_correction,
            _seaice_nudge_energy,
            _h_ML_min, _h_ML_max, we_max,
            _τx, _τy,
            _u, _v, _w,
            _u_bnd, _v_bnd, _w_bnd,
            _GRAD_bnd_x, _GRAD_bnd_y, _GRAD_bnd_z,
            _CURV_x, _CURV_y, _CURV_z,
            _TFLUX_DEN_x, _TFLUX_DEN_y, _TFLUX_DEN_z,
            _SFLUX_DEN_x, _SFLUX_DEN_y, _SFLUX_DEN_z,
            _div,
            _TFLUX_CONV, _TFLUX_CONV_h,
            _SFLUX_CONV, _SFLUX_CONV_h,
            _∇∇T, _∇∇S,
            R, ζ,
            _rad_decay_coes, _rad_absorp_coes,
            Ts_clim_relax_time, Ss_clim_relax_time,
            _Ts_clim, _Ss_clim,
            Nx * Ny, hs, Δzs,
            ( in_flds == nothing ) ? InputFields(datakind, Nx, Ny) : in_flds,
            lays,
            cols,
            ( id == 0 ) ? nothing : AccumulativeVariables(Nx, Ny, Nz_bone),
            ASUM,
            Workspace(Nx=Nx, Ny=Ny, Nz=Nz_bone, shape=:zxy),
            deep_ocn_correction_start_layer,
            qflx_finding_cdm,
            diag,
        )

        updateB!(ocn)
        updateFLDO!(ocn)

        if do_convective_adjustment
            @loop_hor ocn i j let
                OC_doConvectiveAdjustment!(ocn, i, j)
            end
        end

        calTEMP!(ocn)

        return ocn
    end

end


