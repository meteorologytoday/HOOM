function setupBinding!(
    slave :: DynSlave,
)

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    s  = m.state
    du_there = sd.data_units

    bindings = (
        ([:TO_TMD], DataUnit(:u_c, :cU, :xyz, s.u_c, false), :u_c),
        ([:TO_TMD], DataUnit(:v_c, :cV, :xyz, s.v_c, false), :v_c),
        ([:FR_TMD], DataUnit(:b_c, :cT, :xyz, s.b_c, false), :b_c),
        ([:TO_MAS, :TEST], DataUnit(:Φ  , :sT, :xy , s.Φ,   false), :Φ  ),
    )

    for (group_labels, here, there_key) in bindings
       
        println("Doing : ", here.id, "; ", du_there[there_key].id) 

        createBinding!(
            de,
            here,
            du_there[there_key],
            Colon(),
            Colon(),
            Colon(),
            Colon(),
            labels = group_labels,
        )
    end
end
