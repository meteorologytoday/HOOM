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
        ((:u_total_c, :cU, :xyz, bd[:u_total_c], noXdim), :u_total_c),
        ((:v_total_c, :cV, :xyz, bd[:v_total_c], noXdim), :v_total_c),
        ((:B_c,       :cT, :xyz, bd[:B_c],       noXdim), :B_c),


        
        # T, S, FLDO and such
        ((:X,    :fT, :zxy, s.X,    hasXdim), :X   ),
        ((:X_ML, :sT, :xy , s.X_ML, hasXdim), :X_ML),
        ((:h_ML, :sT, :xy , s.h_ML, noXdim),  :h_ML),
        ((:FLDO, :sT, :xy , s.FLDO, noXdim),  :FLDO),
        ((:b,    :fT, :zxy, s.b,    noXdim),  :b   ),
        ((:b_ML, :sT, :xy,  s.b_ML, noXdim),  :b_ML),
        ((:B,    :fT, :zxy, s.B,    noXdim),  :B   ),

        # forcings
        ((:SWFLX,  :sT, :xy,  f.swflx,  noXdim), :SWFLX),
        ((:NSWFLX, :sT, :xy,  f.nswflx, noXdim), :NSWFLX),
        ((:u_U,    :fU, :zxy, f.u_U,  noXdim),   :u_U ),
        ((:v_V,    :fV, :zxy, f.v_V,  noXdim),   :v_V ),
        ((:w_W,    :fW, :zxy, f.w_W,  noXdim),   :w_W ),

    )

    #=
    groups = Dict(
        :FR_DYN => (:u_c, :v_c,),
        :TO_DYN => (:b_c,),
        :BND    => (:X, :X_ML,),
        :FR_MAS => (:SWFLX, :NSWFLX),
        :TO_MAS => (:X, :X_ML, :h_ML, :FLDO),
    )
    =#
    #println("createBinding..")

    for (here_args, there_key) in bindings
       
        here = DataUnit(here_args...)
        println("Doing : ", here.id, "; ", du_there[there_key].id, "; size: ", size(here.data)) 
        
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
            there_key,
            here,
            du_there[there_key],
            here_pull_yrng,
            there_pull_yrng,
            here_push_yrng,
            there_push_yrng,
        )
    end

    #=
    for (label, du_ids) in groups
        for id in du_ids
            addToGroup!(de, id, label)
        end
    end
    =#

    #println("done.")
end



function calCoarseBuoyancyPressure!(
    slave :: TmdSlave,
)

    calAverage_f2c!(
        slave.va,
        PermutedDimsArray(slave.model.state.B, (2, 3, 1)),
        slave.buffer_data[:B_c],
    )

    #println("Sum of state.b: ", sum(slave.model.state.b))
    #println("Sum of b_c: ", sum(slave.buffer_data[:b_c]))

end

function projVelocity!(
    slave :: TmdSlave,
)

    projVertical_c2f!(
        slave.va,
        slave.buffer_data[:u_total_c],
        PermutedDimsArray(slave.model.forcing.u_U, (2, 3, 1)),
    )
 
    projVertical_c2f!(
        slave.va,
        slave.buffer_data[:v_total_c],
        PermutedDimsArray(slave.model.forcing.v_V, (2, 3, 1)),
    )
 
end
