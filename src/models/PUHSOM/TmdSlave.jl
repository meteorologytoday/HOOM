mutable struct TmdSlave

    model        :: Tmd.TmdModel
    ocn_env      :: OcnEnv
    shared_data  :: SharedData
    
    data_exchanger :: DataExchanger
    buffer_data  :: Dict

    y_split_info :: YSplitInfo

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


        data_exchanger = DataExchanger([
           :FR_DYN, :TO_DYN, :BND, :TO_MAS, 
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
        )

    end

end

function setupBinding!(
    slave :: TmdSlave,
)

    println("TmdSlave setupBinding...")

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    du_there = sd.data_units

    bd = slave.buffer_data
    s  = m.state

    ysi = slave.y_split_info

    bindings = (
        ([:FR_DYN], DataUnit(:u_c, :cU, :zxy, bd[:u_c], false), :u_c),
        ([:FR_DYN], DataUnit(:v_c, :cV, :zxy, bd[:v_c], false), :v_c),
        ([:TO_DYN], DataUnit(:b_c, :cT, :zxy, bd[:b_c], false), :b_c),
        
        # T, S, FLDO and such
        ([:BND, :TO_MAS], DataUnit(:X,    :fT, :zxy, s.X,    true), :X   ),
        ([:BND, :TO_MAS], DataUnit(:X_ML, :sT, :xy , s.X_ML, true), :X_ML),
    )

    println("createBinding..")

    for (group_labels, here, there_key) in bindings
       
        println("Doing : ", here.id, "; ", du_there[there_key].id) 
        here_yrng  = Colon()
        if here.grid in (:fV, :cV, :sV)
            there_yrng = ysi.pull_fr_rng[1]:(ysi.pull_fr_rng[end]+1)
        else
            there_yrng = ysi.pull_fr_rng
        end

        createBinding!(
            de,
            here,
            du_there[there_key],
            here_yrng,
            there_yrng,
            labels = group_labels,
        )
    end

    println("done.")
end
