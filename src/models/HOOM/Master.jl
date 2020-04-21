function init!(
    restart_file,
    ocn_env,
    ocn_state,
)


    parallization_info = decide_partition_of_cores(ocn_env)
    data_manager = create_datamanager(ocn_state)
    
    dyn_slaves = create_dynslave(parallization_info, data_manager, env)
    tcr_slaves = create_tcrslave(p_info, data_manager, env)
    mld_slaves = create_mldslave(p_info, data_manager, env)
    
    slaves_init!(data_manager)

    if restart_file != nothing
        loadRestart(restart_file, data_manager)
    end
    
    syncOcean(from = :master, to=:all, data_manager)
end

function run!(
    model,
    write_restart,
)

    load ... 
    substep_dyn
    substep_tcr
    substep_mld

    # Currently mld_core does not need info from
    # dyn_core so we do not need to pass dyn fields
    # to mld core
    for t=1:substep_dyn
        step_dyn(dyn_slaves)
    end
    syncOcean(from = :dyn_slave, to = :tcr_slave, data_manager)
 
    # this involves passing tracer through boundaries
    # so need to sync every time after it evolves
    for t=1:substep_tcr
        step_tcr(tcr_slaves)
        if t != substep_tcr
            syncOcean(from = :tcr_slave, to = :tcr_slave, data_manager)
        else
            syncOcean(from = :tcr_slave, to = :mld_slave, data_manager)
        end
    end

    # Supposely MLD dynamics changes dyn and tcr fields vertically
    # so it only sync by the end of simulation and sync to
    # all other components
    for t=1:substep_tcr
        step_mld(mld_slaves)
    end
    syncOcean(from = :mld_slave, to = :all, data_manager)
   
    if write_restart
        writeRestart(
            dyn_slave,
            tcr_slave,
            mld_slave, 
        )
    end
     
end





function loadData!()
end
