function init!(
    restart_file,
)

    ocn_env = readRestart(restart_file)
    
    ocn_state          = OcnState(ocn_env)
    shared_data       = create_datamanager(ocn_state)
    parallization_info = decide_partition_of_cores(ocn_env)

    
    dyn_slaves = create_dynslave(parallization_info, shared_data, env)
    tcr_slaves = create_tcrslave(p_info, shared_data, env)
    mld_slaves = create_mldslave(p_info, shared_data, env)
    
    slaves_init!(shared_data)

    if restart_file != nothing
        loadRestart(restart_file, shared_data)
    end
    
    syncOcean(from = :master, to=:all, shared_data)
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
    syncOcean(from = :dyn_slave, to = :tcr_slave, shared_data)
 
    # this involves passing tracer through boundaries
    # so need to sync every time after it evolves
    for t=1:substep_tcr
        step_tcr(tcr_slaves)
        if t != substep_tcr
            syncOcean(from = :tcr_slave, to = :tcr_slave, shared_data)
        else
            syncOcean(from = :tcr_slave, to = :mld_slave, shared_data)
        end
    end

    # Supposely MLD dynamics changes dyn and tcr fields vertically
    # so it only sync by the end of simulation and sync to
    # all other components
    for t=1:substep_mld
        step_mld(mld_slaves)
    end
    syncOcean(from = :mld_slave, to = :all, shared_data)
   
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
