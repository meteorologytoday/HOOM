function setupBinding!(
    slave :: DynSlave,
)

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    s  = m.state
    c  = m.core
    du_there = sd.data_units

    bindings = (
        (DataUnit(:u_total, :cU, :xyz, s.u_total, false), :u_total_c),
        (DataUnit(:v_total, :cV, :xyz, s.v_total, false), :v_total_c),
        (DataUnit(:B,       :cT, :xyz, s.B,       false), :B_c      ),
        (DataUnit(:Φ,       :sT, :xy , s.Φ,       false), :Φ        ),
        (DataUnit(:∂B∂x,    :cU, :xyz, c.∂B∂x,    false), :∂B∂x     ),
        (DataUnit(:∂B∂y,    :cV, :xyz, c.∂B∂y,    false), :∂B∂y     ),
    )
    #=
    groups = Dict(
        :FR_TMD => (:b_c,),
        :TO_TMD => (:u_c, :v_c,),
        :TO_MAS => (:Φ,),
        :TEST   => (:Φ,),
    )
    =#
    for (here, there_key) in bindings
      
        # The name of the binding uses
        # the key in the shared_data
 
        createBinding!(
            de,
            there_key,
            here,
            du_there[there_key],
            Colon(),
            Colon(),
            Colon(),
            Colon(),
        )
    end

    #=
    for (label, du_ids) in groups
        for id in du_ids
            addToGroup!(de, id, label)
        end
    end
    =#

end
