mutable struct Master

    dyn_slave       # Non parallizable
    tcr_slaves      # parallizable
    mld_slaves      # parallizable

    CPU_id          # CPU ids
    
    data_manager    # manage underlying data exchange
    data_recorder   # history file tool
   
    ocn
 
    function Master(
        
    )


        return new(
            dynamic_core    
        )
    end

end


mutable struct ParallizationInfo
end

mutable struct TcrSlave
    data_manager
    core 
end

mutable struct DynSlave
    data_manager 
    core
end

mutable struct MLDSlave
    data_manager 
    core
end


mutable struct DataManager
    binding
    data
end

mutable struct OceanState
    dyn_state
    tcr_state
    mld_state
end

