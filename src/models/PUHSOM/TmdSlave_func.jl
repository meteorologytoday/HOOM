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
            here,
            du_there[there_key],
            here_pull_yrng,
            there_pull_yrng,
            here_push_yrng,
            there_push_yrng,
            labels = group_labels,
        )
    end

    println("done.")
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
        slave.model.state.u_U,
    )
 
    projVertical_c2f!(
        slave.va,
        slave.buffer_data[:v_c],
        slave.model.state.v_V,
    )
 
end
