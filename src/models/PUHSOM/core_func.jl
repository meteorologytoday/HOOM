function init!(
    ocn_env = Union{String, OceanEnv},
)

    # TODO
    ocn_env = (typeof(ocn_env) == String) ? loadOcnEnv(env) : ocn_env

    shared_data    = SharedData(ocn_env)
    job_dist_info  = JobDistributionInfo(ocn_env; overlap=2)

    model = Model(
        ocn_env,
        shared_data,
        job_dist_info,
    )

    println("Register Shared Data") 
    registerSharedData!(model)
    
    # Potential redundant: seems like shared_data is already doing this
    #ocn_state      = OcnState(shared_data)


    println("Creating slaves on nodes")
    @sync let

        @spawnat job_dist_info.dyn_slave_pid let
            global dyn_slave = PUHSOM.DynSlave(ocn_env, shared_data)

            PUHSOM.setupBinding!(dyn_slave)

        end

        for (p, pid) in enumerate(job_dist_info.tmd_slave_pids)
            @spawnat pid let
                global tmd_slave = PUHSOM.TmdSlave(
                    ocn_env,
                    shared_data,
                    job_dist_info.y_split_infos[p],
                )

                PUHSOM.setupBinding!(tmd_slave)
            end
        end
    end

    println("Slave created and data exchanger is set.")
    println("TODO: Skip read restart file for now.")

    #=
    if restart_file != nothing
        loadRestart(restart_file, shared_data)
    end
    =#

    syncTmd!(model, :TO_DYN, :S2M)
    return model

end


function stepModel!(
    model :: Model,
    write_restart :: Bool,
)

    env = model.env

    # Sync
    syncDyn!(model, :FR_TMD, :M2S)
   
     
    # Currently tmd_core does not need info from
    # dyn_core so we do not need to pass dyn fields
    # to mld core
#    for t=1:env.substep_dyn
#        Dyn.stepModel!(dyn_slaves)
#    end
 
#=
    @sync @spawnat model.job_dist_info.dyn_slave_pid let
            for t=1:env.substep_dyn
                Dyn.stepModel!(dyn_slave.model)
            end
    end
=# 
    
    # Sending updated velocity to tmd
    syncDyn!(model, :TO_TMD, :S2M)
    syncTmd!(model, :FR_DYN, :M2S)

    ##### tmd slave should distribute u,v to fine grids here #####
    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid projVelocity!(tmd_slave)
    end

    # this involves passing tracer through boundaries
    # so need to sync every time after it evolves
    for t=1:env.substep_tmd

        #step_tmd(tcr_slaves)

        syncTmd!(model, :BND, :S2M)
        syncTmd!(model, :BND, :M2S)

    end
    
    ##### tmd slave should calcaulte b of coarse grid #####
    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid let
            PUHSOM.calCoarseBuoyancy!(tmd_slave)
        end
    end

    syncTmd!(model, :TO_MAS, :S2M)
    syncDyn!(model, :TO_MAS, :S2M)

    #= 
    if write_restart
        writeRestart(
            dyn_slave,
            tcr_slave,
            mld_slave, 
        )
    end
    =#
end


function registerSharedData!(model::Model)

    descs_X = (
        (:X,     :fT, :zxy, Float64),
        (:X_ML,  :sT, :xy,  Float64),
    )
 
    descs_noX = (

        (:h_ML,  :sT, :xy,  Float64),
        (:FLDO,  :sT, :xy,    Int64),
        
        # These are used by dyn_core 
        (:u_c,   :cU, :xyz, Float64),
        (:v_c,   :cV, :xyz, Float64),
        (:b_c,   :cT, :xyz, Float64),
        (:Φ  ,   :sT, :xy,  Float64),

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
        (:MLT,     :sT, :xy,  Float64),
    ) 


    for (id, grid, shape, dtype) in descs_X
        regVariable!(model.shared_data, model.env, id, grid, shape, dtype, has_Xdim=true)
    end

    for (id, grid, shape, dtype) in descs_noX
        regVariable!(model.shared_data, model.env, id, grid, shape, dtype, has_Xdim=false)
    end

end


push_pull_relation = Dict(
    :S2M => :PUSH,
    :M2S => :PULL,
)
function syncTmd!(
    model :: Model,
    group_label :: Symbol,
    direction :: Symbol,
)

    direction = push_pull_relation[direction]

    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid let
            PUHSOM.syncData!(tmd_slave.data_exchanger, group_label, direction)
        end
    end
end

function syncDyn!(
    model:: Model,
    group_label :: Symbol,
    direction :: Symbol,
)

    direction = push_pull_relation[direction]

    @sync let
        @spawnat model.job_dist_info.dyn_slave_pid let
            PUHSOM.syncData!(dyn_slave.data_exchanger, group_label, direction)
        end
    end
end

function loadData!()
end
