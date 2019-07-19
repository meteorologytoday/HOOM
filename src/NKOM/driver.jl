# used by master
mutable struct SubOCCInfo
    sub_Nxs    :: AbstractArray{Int64, 1}
    beg_x_idxs :: AbstractArray{Int64, 1}
end

function makeSubOCC(
    occ :: OceanColumnCollection,
    block_id   :: Integer,
    beg_x_idx  :: Integer,
    sub_Nx     :: Integer,
)

    global master_occ = occ
    global rng2 = [beg_x_idx:beg_x_idx+sub_Nx-1, Colon()]
    global rng3 = [beg_x_idx:beg_x_idx+sub_Nx-1, Colon(), Colon()]
   
 
    global worker_occ = OceanColumnCollection(
        id             = block_id,
        gridinfo_file  = nothing,
        Nx             = sub_Nx,
        Ny             = occ.Ny,
        zs_bone        = occ.zs_bone,
        Ts             = occ.Ts[rng3...],
        Ss             = occ.Ss[rng3...],
        K_T            = occ.K_T,
        K_S            = occ.K_S,
        T_ML           = occ.T_ML[rng2...],
        S_ML           = occ.S_ML[rng2...],
        h_ML           = occ.h_ML[rng2...],
        h_ML_min       = occ.h_ML_min[rng2...],
        h_ML_max       = occ.h_ML_max[rng2...],
        we_max         = occ.we_max,
        Ts_clim_relax_time = occ.Ts_clim_relax_time,
        Ss_clim_relax_time = occ.Ss_clim_relax_time,
        Ts_clim        = ( occ.Ts_clim != nothing ) ? occ.Ts_clim[rng3...] : nothing,
        Ss_clim        = ( occ.Ss_clim != nothing ) ? occ.Ss_clim[rng3...] : nothing,
        mask           = occ.mask[rng2...],
        topo           = occ.topo[rng2...],
        fs             = occ.fs[rng2...],
        ϵs             = occ.ϵs[rng2...],
        in_flds        = SubInputFields(occ.in_flds, rng2...)
    )


end

function init(
    occ::OceanColumnCollection,
)

    println("Number of workers: ", nworkers())

    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    sub_Nx = ceil(Integer, occ.Nx / nworkers())
    sub_Nxs = [sub_Nx for block_id = 1:nworkers()]
    sub_Nxs[end] = occ.Nx - (nworkers()-1) * sub_Nx

    beg_x_idxs = [sub_Nx * (block_id - 1) + 1 for block_id = 1:nworkers()]

    global sub_occ_info = SubOCCInfo(sub_Nxs, beg_x_idxs)

    @sync for (i, p) in enumerate(workers())
        
        # We have P processors, N workers, N blocks
        # Block ids are numbered from 1 to N
        @spawnat p makeSubOCC(occ, i, beg_x_idxs[i], sub_Nxs[i])

    end
end

function run!(
    occ    :: OceanColumnCollection;
    use_qflx      :: Bool,
    use_h_ML      :: Bool,
    Δt            :: Float64,
    diffusion_Δt  :: Float64,
    relaxation_Δt :: Float64,
    do_diffusion  :: Bool = true, 
    do_relaxation :: Bool = true, 
    do_convadjust :: Bool = true, 
)

    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))
    
    for (i, p) in enumerate(workers())

        @spawnat p NKOM.stepOceanColumnCollection!(
            worker_occ,
            use_qflx      = use_qflx,
            use_h_ML      = use_h_ML,
            Δt            = Δt,
            diffusion_Δt  = diffusion_Δt,
            relaxation_Δt = relaxation_Δt,
            do_diffusion  = do_diffusion,
            do_relaxation = do_relaxation,
            do_convadjust = do_convadjust,
        )

    end
end

function syncToMaster(occ::OceanColumnCollection)

    global master_occ
    
    (occ.id != 0) || throw(ErrorException("`id` should not be 0 (master)."))

    master_occ.Ts[rng3...] = occ.Ts
    master_occ.Ss[rng3...] = occ.Ss

    master_occ.FLDO[rng2...] = occ.FLDO
    master_occ.T_ML[rng2...] = occ.T_ML 
    master_occ.S_ML[rng2...] = occ.S_ML
    master_occ.h_ML[rng2...] = occ.h_ML

    #master_occ.T_ML[rng2...] .+= 10.0 * rand()
    #println( pointer_from_objref(master_occ.T_ML))
end


function sync!(
    occ :: OceanColumnCollection;
)
    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    @sync for (i, p) in enumerate(workers())
        @spawnat p NKOM.syncToMaster(worker_occ)
    end


end
