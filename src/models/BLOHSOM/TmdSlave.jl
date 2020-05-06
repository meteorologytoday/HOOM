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

        gi = PolelikeCoordinate.genGridInfo(
            ocn_env.hrgrid,
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
            gi       = gi,
            Δt       = ocn_env.Δt,
            substeps = ocn_env.substeps_tmd,
            z_bnd    = ocn_env.z_bnd_f,
            topo     = ocn_env.topo,
            mask2    = ocn_env.mask2,
            Kh_X     = ocn_env.Kh_X,
            Kv_X     = ocn_env.Kv_X,
            we_max   = ocn_env.we_max,
            R        = ocn_env.R,
            ζ        = ocn_env.ζ,
            MLT_rng  = ocn_env.MLT_rng,
            NX_passive = ocn_env. NX_passive,
            t_X_wr     = ocn_env.t_X_wr,
            X_wr       = X_wr,
            MLT_scheme = ocn_env.MLT_scheme,
            radiation_scheme = Symbol(ocn_env.radiation_scheme),
            convective_adjustment = ocn_env.convective_adjustment,
            use_Qflux     = ocn_env.use_Qflux,
            finding_Qflux = ocn_env.finding_Qflux,
        )         


        data_exchanger = DataExchanger([
           :FR_DYN, :TO_DYN, :BND, :TO_MAS, :FR_MAS, 
        ])

        buffer_data = Dict(
            :u_total_c => zeros(Float64, Nx, Ny,   Nz_c),
            :v_total_c => zeros(Float64, Nx, Ny+1, Nz_c),
            :B_c => zeros(Float64, Nx, Ny,   Nz_c),
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

