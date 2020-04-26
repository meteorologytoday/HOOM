function init!(
    model :: Model;
    ocn_env = Union{String, OceanEnv},
)

    # TODO
    ocn_env = (typeof(ocn_env) == String) ? loadOcnEnv(env) : ocn_env

    shared_data    = SharedData(ocn_env)
    job_dist_info  = JobDistributionInfo(ocn_env; overlap=2)
    
    registerSharedData!(shared_data)
    
    # Potential redundant: seems like shared_data is already doing this
    #ocn_state      = OcnState(shared_data)

    @sync let
        @spawnat job_dist_info.dyn_slave_pid let
            include("DynSlave.jl")
            global dyn_slave = DynSlave(ocn_env, shared_data)
        end

        for pid in job_dist_info.tmd_slave_pids
            @spawnat pid let
                include("TmdSlave.jl")
                global tmd_slave = TmdSlave(ocn_env, shared_data)
            end
        end
    end

    println("Stop here")
    readline()

    


   
    
    slaves_init!(shared_data)

    if restart_file != nothing
        loadRestart(restart_file, shared_data)
    end
    
    syncOcean(from = :master, to=:all, shared_data)
end


#=
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

    # Supposedly MLD dynamics changes dyn and tcr fields vertically
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

=#

function registerSharedData!(sd::SharedData)

    # F
    descs = (

        # These are intensively used
        # in tcr_mld_core. zxy will be
        # much memory friendly
        (:T,     :fT, :zxy, Float64),
        (:S,     :fT, :zxy, Float64),
        (:b,     :fT, :zxy, Float64),
        (:T_ML,  :sT, :xy,  Float64),
        (:S_ML,  :sT, :xy,  Float64),
        (:h_ML,  :sT, :xy,  Float64),
        (:FLDO,  :sT, :xy,    Int64),
        
        # These are used by dyn_core 
        (:u_c,   :cU, :xyz, Float64),
        (:v_c,   :cV, :xyz, Float64),
        (:Phi,   :sT, :xy,  Float64),

        # Forcings and return fluxes to coupler
        (:SWFLX,   :sT, :xy,  Float64),
        (:NSWFLX,  :sT, :xy,  Float64),
        (:TAUX,    :sT, :xy,  Float64),
        (:TAUY,    :sT, :xy,  Float64),
        (:IFRAC,   :sT, :xy,  Float64),
        (:FRWFLX,  :sT, :xy,  Float64),
        (:VSFLX,   :sT, :xy,  Float64),
        (:QFLX_T,  :sT, :xy,  Float64),
        (:QFLX_S,  :sT, :xy,  Float64),
        (:T_CLIM,  :sT, :xy,  Float64),
        (:S_CLIM,  :sT, :xy,  Float64),
        (:MLD,     :sT, :xy,  Float64),
    ) 


    for (id, grid, shape, dtype) in descs
        regVariable!(sd, id, grid, shape, dtype)
    end




end


function loadData!()
end
