
#
# Important concepts
#
#
#
# "States" are what is managed by data_manager to sync between processes
# Anything that does not have to share should not be in "states"
#
# DataManager must manage all variable in "States". 

mutable struct Master

    dyn_slave       # Non parallizable
    tcr_slaves      # parallizable
    mld_slaves      # parallizable

    pid_id          # Process ids
    
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

