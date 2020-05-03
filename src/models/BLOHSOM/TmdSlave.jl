mutable struct TmdSlave

    model        :: Tmd.TmdModel
    ocn_env      :: OcnEnv
    shared_data  :: SharedData
    
    data_exchanger :: DataExchanger
    buffer_data  :: Dict

    y_split_info :: YSplitInfo

    va       :: Union{VerticalAverager, Nothing}

    function TmdSlave(
        ocn_env      :: OcnEnv,
        shared_data  :: SharedData,
        y_split_info :: YSplitInfo,
    )
       
        """ WARNING !!!!! Ny should be length(sub_yrng) """


        Nx = ocn_env.Nx
        Ny = length(y_split_info.pull_fr_rng)
        Nz_f = ocn_env.Nz_f
        Nz_c = ocn_env.Nz_c
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
        
        X_wr = zeros(Float64, Nz_f, Nx, Ny, NX)
        X_wr .= NaN
        for (x, t) in enumerate(ocn_env.t_X_wr)
            if ! isnan(t)
                Dataset(ocn_env.X_wr_file[x], "r") do ds
                    X_wr[:, :, :, x] = PermutedDimsArray( nomissing( ds[ocn_env.X_wr_varname][:, y_split_info.pull_fr_rng, :] ), (3, 1, 2))
                end    
            end
        end
        # ===== [END] X_wr =====


        # ===== [BEGIN] Making vertical averager =====
        va = VerticalAverager(
            z_bnd_f             = ocn_env.z_bnd_f,
            height_level_counts = ocn_env.height_level_counts,
        )
        
        # ===== [END] Making vertical averager =====
 
        model = Tmd.TmdModel(
            gi      = gi,
            Δt      = ocn_env.Δt / ocn_env.substep_tmd,
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


        data_exchanger = DataExchanger([
           :FR_DYN, :TO_DYN, :BND, :TO_MAS, :FR_MAS, 
        ])

        buffer_data = Dict(
            :u_c => zeros(Float64, Nz_c, Nx, Ny  ),
            :v_c => zeros(Float64, Nz_c, Nx, Ny+1),
            :b_c => zeros(Float64, Nz_c, Nx, Ny  ),
        )



        return new(
            model, 
            ocn_env,
            shared_data,
            data_exchanger,
            buffer_data,
            y_split_info,
            va,
        )

    end

end

