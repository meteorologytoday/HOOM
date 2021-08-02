mutable struct Core

    gf          :: Union{PolelikeCoordinate.GridFile, Nothing}


    gd       :: PolelikeCoordinate.Grid
    gd_slab  :: PolelikeCoordinate.Grid

    mask_sT :: AbstractArray{Float64, 3}
    mask_T :: AbstractArray{Float64, 3}
 
  
    amo      :: Union{AdvancedMatrixOperators, Nothing}
    amo_slab :: Union{AdvancedMatrixOperators, Nothing}

    function Core(
        ev :: Env;
    )

        gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
            ev.gf_filename;
            R  = 6371229.0,
            Ω  = 2π / (86400 / (1 + 365/365)),
        )


        gd      = PolelikeCoordinate.genGrid(gf, ev.z_w ; sub_yrng=ev.sub_yrng) 
        gd_slab = PolelikeCoordinate.genGrid(gf, [0, -1.0]; sub_yrng=ev.sub_yrng) 

        mask_sT = reshape(gf.mask[:, ev.sub_yrng], 1, :, ev.Ny)
        mask_T  = repeat( mask_sT, outer=(gd.Nz, 1, 1) )
 
        writeLog("Mask shape sT {:s}", string(size(mask_sT)); force=true)       
        writeLog("Mask shape T  {:s}", string(size(mask_T)); force=true)       

        writeLog("gd      shape: {:d} {:d} {:d}", gd.Nz, gd.Nx, gd.Ny; force=true)       
        writeLog("gd_slab shape: {:d} {:d} {:d}", gd_slab.Nz, gd_slab.Nx, gd_slab.Ny; force=true)       
        @time amo_slab = AdvancedMatrixOperators(;
            gd = gd_slab,
            mask_T = mask_sT,
        )

        @time amo = AdvancedMatrixOperators(;
            gd = gd,
            mask_T = mask_T,
        )



        # Build Advection Matrix

        # Build Radiation Matrix
        swflx_factor_W = (1.0 - ev.R) * exp.(gd.z_W / ev.ζ2)
        lwflx_factor_W =        ev.R  * exp.(gd.z_W / ev.ζ1)

        # Bottom absorbs everything
        swflx_factor_W[end, :, :] .= 0.0
        lwflx_factor_W[end, :, :] .= 0.0

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

        swflx_conv = - amo.T_DIVz_W * spdiagm(0 => view(swflx_factor_W, :)) * W_broadcast_sT
        lwflx_conv = - amo.T_DIVz_W * spdiagm(0 => view(lwflx_factor_W, :)) * W_broadcast_sT 

        return new(
            gf,
            gd,
            gd_slab,

            mask_sT,
            mask_T,
            amo,
        )
    end

end


