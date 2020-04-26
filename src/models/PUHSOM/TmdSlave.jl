mutable struct TmdSlave

    model        :: Tmd.TmdModel

    ocn_env      :: OcnEnv
    shared_data  :: SharedData

    function TmdSlave(
        ocn_env      :: OcnEnv,
        shared_data  :: SharedData,
        y_split_info :: YSplitInfo,
    )
       
        Nx = ocn_env.Nx
        Ny = ocn_env.Ny
        Nz = ocn_env.Nz_f
        NX = ocn_env.NX

        mi = ModelMap.MapInfo{Float64}(ocn_env.hrgrid_file)
        gi = PolelikeCoordinate.CurvilinearSphericalGridInfo(;
            R=Re,
            Ω=Ωe,
            Nx=mi.nx,
            Ny=mi.ny,
            c_lon=mi.xc,
            c_lat=mi.yc,
            vs_lon=mi.xv,
            vs_lat=mi.yv,
            area=mi.area,
            angle_unit=:deg,
            sub_yrng=y_split_info.pull_fr_rng,
        )
       
        # ===== [BEGIN] X_wr =====
        
        X_wr = zeros(Float64, Nz, Nx, Ny, NX)
        X_wr .= NaN
        for (x, t) in enumerate(ocn_env.t_X_wr)
            if ! isnan(t)
                Dataset(ocn_env.X_wr_file[x], "r") do ds
                    X_wr[:, :, :, x] = PermutedDimsArray( nomissing( ds[ocn_env.X_wr_varname][:, :, :] ), (3, 1, 2))
                end    
            end
        end
        # ===== [END] X_wr =====

 
        model = Tmd.TmdModel(
            gi      = gi,
            Δt      = ocn_env.Δt,
            z_bnd   = ocn_env.z_bnd_f,
            topo    = ocn_env.topo,
            mask2   = ocn_env.mask2,
            Kh_X    = ocn_env.Kh_X,
            Kv_X    = ocn_env.Kv_X,
            we_max  = ocn_env.we_max,
            R       = ocn_env.R,
            ζ       = ocn_env.ζ,
            MLT_rng = ocn_env.MLT_rng,
            NX_passive = ocn_env. NX_passive,
            t_X_wr     = ocn_env.t_X_wr,
            X_wr       = X_wr,
            radiation_scheme = Symbol(ocn_env.radiation_scheme),
            convective_adjustment = ocn_env.convective_adjustment,
        )         

        return new(
            model, 
            ocn_env,
            shared_data,
        )

    end

end
