function init!(
    ocn_env
)


    parallization_info = decide_partition_of_cores(ocn_env)
    data_manager = create_datamanager()
    
    dyn_slaves = create_dynslave(data_manager, env)
    tcr_slaves = create_tcrslave(data_manager, env)
    mld_slaves = create_mldslave(data_manager, env)
    
    slaves_init(la)
    
    
 

end






function loadData!()
end
