mutable struct Core

    amo :: Union{AdvancedMatrixOperators, Nothing}

    function Core(
        ev :: Env;
    )
    
    println("Test if slab AMO is there")
    @time amo = AdvancedMatrixOperators(;
        gd = ev.gd_slab,
        mask_T = ev.mask_sT,
    )


        # Build Advection Matrix

        # Build Radiation Matrix
        swflx_factor_W = (1.0 - ev.R) * exp.(ev.gd.z_W / ev.ζ2)
        lwflx_factor_W =        ev.R  * exp.(ev.gd.z_W / ev.ζ1)

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

        return new(
            amo,
        )
    end

end


