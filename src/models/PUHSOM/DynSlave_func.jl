function setupBinding!(
    slave :: DynSlave,
)

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    s  = m.state
    du_there = sd.data_units

    bindings = (
        (DataUnit(:u_c, :cU, :xyz, s.u_c, false), :u_c),
        (DataUnit(:v_c, :cV, :xyz, s.v_c, false), :v_c),
        (DataUnit(:b_c, :cT, :xyz, s.b_c, false), :b_c),
        (DataUnit(:Φ  , :sT, :xy , s.Φ,   false), :Φ  ),
    )

    groups = Dict(
        :FR_TMD => (:b_c,),
        :TO_TMD => (:u_c, :v_c,),
        :TO_MAS => (:Φ,),
        :TEST   => (:Φ,),
    )
 
    for (here, there_key) in bindings
       
        createBinding!(
            de,
            here.id,
            here,
            du_there[there_key],
            Colon(),
            Colon(),
            Colon(),
            Colon(),
        )
    end

    for (label, du_ids) in groups
        for id in du_ids
            addToGroup!(de, id, label)
        end
    end


end
