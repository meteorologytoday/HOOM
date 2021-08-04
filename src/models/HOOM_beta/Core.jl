mutable struct Core

    gf        :: Union{PolelikeCoordinate.GridFile, Nothing}

    gd        :: PolelikeCoordinate.Grid
    gd_slab   :: PolelikeCoordinate.Grid

    mask_sT   :: AbstractArray{Float64, 3}
    mask_T    :: AbstractArray{Float64, 3}
 
    amo_slab  :: Union{AdvancedMatrixOperators, Nothing}
    amo       :: Union{AdvancedMatrixOperators, Nothing}

    mtx       :: Dict

    vd        :: VerticalDiffusion

    cdatam    :: Union{CyclicDataManager, Nothing}

    function Core(
        ev :: Env,
    )

        cfg = ev.config

        gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
            cfg[:domain_file];
            R  = 6371229.0,
            Ω  = 2π / (86400 / (1 + 365/365)),
        )


        gd      = PolelikeCoordinate.genGrid(gf, cfg[:z_w] ; sub_yrng=ev.sub_yrng) 
        gd_slab = PolelikeCoordinate.genGrid(gf, [0, -1.0]; sub_yrng=ev.sub_yrng) 

        mask_sT = reshape(gf.mask[:, ev.sub_yrng], 1, ev.Nx, ev.Ny)
        mask_T  = repeat( mask_sT, outer=(gd.Nz, 1, 1) )
 
        @time amo_slab = AdvancedMatrixOperators(;
            gd = gd_slab,
            mask_T = mask_sT,
        )

        @time amo = AdvancedMatrixOperators(;
            gd = gd,
            mask_T = mask_T,
        )



        # Build Advection Matrix


        function build!(id_mtx, idx)
            local result
#            rows = size(id_mtx)[1]
            
            # using transpose speeds up by 100 times 
            tp = transpose(id_mtx) |> sparse
            result = transpose(tp[:, view(idx, :)]) |> sparse
            dropzeros!(result)

            idx .= 0 # clean so that debug is easir when some girds are not assigned
            return result
        end

        # Build the matrix to broadcast sW to W grid
        # Notice that we cannot use W_pts because gd_slab has two W layers 
        num_sT = reshape( collect(1:amo_slab.bmo.T_pts), amo_slab.bmo.T_dim...)
        mapping_T = repeat(num_sT, outer=(ev.Nz+1, 1, 1))
        W_broadcast_sT = build!(amo_slab.bmo.T_I_T, mapping_T)

        # Build Radiation Matrix
        swflx_factor_W =  cfg[:rad_R]  * exp.(gd.z_W / cfg[:rad_ζ1]) + (1.0 - cfg[:rad_R]) * exp.(gd.z_W / cfg[:rad_ζ2])
        swflx_factor_W[end, :, :] .= 0.0 # Bottom absorbs everything

        nswflx_factor_W = 0.0 * gd.z_W
        nswflx_factor_W[1, :, :] .= 1.0



        # f and ϵ matrices
        f_sT = 2 * gd.Ω * sin.(gd_slab.ϕ_T)
        ϵ_sT = f_sT * 0 .+ cfg[:ϵ]
        D_sT = f_sT.^2 + ϵ_sT.^2
        invD_sT = D_sT.^(-1.0)

        mtx = Dict(
            :T_swflxConv_sT  => - amo.T_DIVz_W * spdiagm(0 => view(swflx_factor_W, :)) * W_broadcast_sT,
            :T_nswflxConv_sT => - amo.T_DIVz_W * spdiagm(0 => view(nswflx_factor_W, :)) * W_broadcast_sT,
            :invD_sT         => invD_sT,
            :f_sT            => f_sT,
            :ϵ_sT            => ϵ_sT,
            :D_sT            => D_sT,
        ) 

        vd = VerticalDiffusion(amo; K_iso=cfg[:Ks_V], K_cva=cfg[:Ks_V_cva])


        cdata_varnames = []

        if cfg[:MLD_scheme] == :datastream
            push!(cdata_varnames, "HBLT")
        end

        if cfg[:Qflx] == :on
            push!(cdata_varnames, "Qflx_T")
            push!(cdata_varnames, "Qflx_S")
        end
        
        if cfg[:weak_restoring] == :on || cfg[:Qflx_finding] == :on
            push!(cdata_varnames, "TEMP")
            push!(cdata_varnames, "SALT")
        end

        if length(cdata_varnames) == 0
            cdatam = nothing
        else
            
            if cfg[:cdata_file] == ""
                throw(ErrorException("Some config require cyclic data forcing file"))
            else
                cdatam = CyclicDataManager(;
                    filename     = cfg[:cdata_file],
                    varname_time = "time", 
                    varnames     = cdata_varnames,
                    beg_time     = 0.0,
                    cyc_time     = 365 * 86400.0,
                    sub_yrng     = ev.sub_yrng,
                )
            end
        end

        return new(
            gf,
            gd,
            gd_slab,

            mask_sT,
            mask_T,

            amo_slab,
            amo,

            mtx,    

            vd,

            cdatam,
        )
    end

end


