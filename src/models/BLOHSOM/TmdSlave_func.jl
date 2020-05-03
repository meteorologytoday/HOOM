function stepModel!(
    slave :: TmdSlave,
)

    Tmd.stepModel!(slave.model)

end

function setupBinding!(
    slave :: TmdSlave,
)

    #println("TmdSlave setupBinding...")

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    du_there = sd.data_units

    bd = slave.buffer_data
    s  = m.state
    f  = m.forcing

    ysi = slave.y_split_info

    hasXdim = true
    noXdim  = false

    bindings = (
        ((:u_c, :cU, :zxy, bd[:u_c], noXdim), :u_c),
        ((:v_c, :cV, :zxy, bd[:v_c], noXdim), :v_c),
        ((:b_c, :cT, :zxy, bd[:b_c], noXdim), :b_c),
        
        # T, S, FLDO and such
        ((:X,    :fT, :zxy, s.X,    hasXdim), :X   ),
        ((:X_ML, :sT, :xy , s.X_ML, hasXdim), :X_ML),
        ((:h_ML, :sT, :xy , s.h_ML, noXdim),  :h_ML),
        ((:FLDO, :sT, :xy , s.FLDO, noXdim),  :FLDO),

        # forcings
        ((:SWFLX,  :sT, :xy, f.swflx,  noXdim), :SWFLX),
        ((:NSWFLX, :sT, :xy, f.nswflx, noXdim), :NSWFLX),

    )

    groups = Dict(
        :FR_DYN => (:u_c, :v_c,),
        :TO_DYN => (:b_c,),
        :BND    => (:X, :X_ML,),
        :FR_MAS => (:SWFLX, :NSWFLX),
        :TO_MAS => (:X, :X_ML, :h_ML, :FLDO),
    )
    
    #println("createBinding..")

    for (here_args, there_key) in bindings
       
        here = DataUnit(here_args...)
        #println("Doing : ", here.id, "; ", du_there[there_key].id) 
        
        here_pull_yrng  = Colon()

        if here.grid in (:fV, :cV, :sV)
            there_pull_yrng = ysi.pull_fr_rng[1]:(ysi.pull_fr_rng[end]+1)
            here_push_yrng  = ysi.push_fr_rng[1]:(ysi.push_fr_rng[end]+1)
            there_push_yrng = ysi.push_to_rng[1]:(ysi.push_to_rng[end]+1)
        else
            there_pull_yrng  = ysi.pull_fr_rng
            here_push_yrng   = ysi.push_fr_rng
            there_push_yrng  = ysi.push_to_rng
        end

        createBinding!(
            de,
            here.id,
            here,
            du_there[there_key],
            here_pull_yrng,
            there_pull_yrng,
            here_push_yrng,
            there_push_yrng,
        )
    end

    for (label, du_ids) in groups
        for id in du_ids
            addToGroup!(de, id, label)
        end
    end

    #println("done.")
end



function calCoarseBuoyancy!(
    slave :: TmdSlave,
)

    println("TODO: somehow need to calculate the avg buoyancy in the grid box where h_ML is in.")

    #mixXinsideFLDO!(state.b)

    calAverage_f2c!(
        slave.va,
        slave.model.state.b,
        slave.buffer_data[:b_c],
        shape=:zxy,
    )

end

function projVelocity!(
    slave :: TmdSlave,
)

    projVertical_c2f!(
        slave.va,
        slave.buffer_data[:u_c],
        slave.model.forcing.u_U,
    )
 
    projVertical_c2f!(
        slave.va,
        slave.buffer_data[:v_c],
        slave.model.forcing.v_V,
    )
 
end
