# used by master
mutable struct SubOCCInfo
    sub_Nys    :: AbstractArray{Int64, 1}
    beg_y_idxs :: AbstractArray{Int64, 1}
end

function makeSubOCC(
    occ :: OceanColumnCollection,
    block_id   :: Integer,
    beg_y_idx  :: Integer,
    sub_Ny     :: Integer,
)

    global master_occ = occ
    global rng2 = [Colon(), beg_y_idx:beg_y_idx+sub_Ny-1]
    global rng3 = [Colon(), Colon(), beg_y_idx:beg_y_idx+sub_Ny-1]
    global master_sub_in_flds = SubInputFields(occ.in_flds, rng2...)
   
    global worker_occ = OceanColumnCollection(
        id             = block_id,
        gridinfo_file  = nothing,
        Nx             = occ.Nx,
        Ny             = sub_Ny,
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
        R              = occ.R,
        ζ              = occ.ζ,
        Ts_clim_relax_time = occ.Ts_clim_relax_time,
        Ss_clim_relax_time = occ.Ss_clim_relax_time,
        Ts_clim        = ( occ.Ts_clim != nothing ) ? occ.Ts_clim[rng3...] : nothing,
        Ss_clim        = ( occ.Ss_clim != nothing ) ? occ.Ss_clim[rng3...] : nothing,
        mask           = occ.mask[rng2...],
        topo           = occ.topo[rng2...],
        fs             = occ.fs[rng2...],
        ϵs             = occ.ϵs[rng2...],
        in_flds        = InputFields(:local, occ.Nx, sub_Ny),
        arrange        = :zxy,
    )


end

function init(
    occ::OceanColumnCollection,
)

    println("Number of workers: ", nworkers())

    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    sub_Ny = ceil(Integer, occ.Nx / nworkers())
    sub_Nys = [sub_Ny for block_id = 1:nworkers()]
    sub_Nys[end] = occ.Ny - (nworkers()-1) * sub_Ny

    beg_y_idxs = [sub_Ny * (block_id - 1) + 1 for block_id = 1:nworkers()]

    global sub_occ_info = SubOCCInfo(sub_Nys, beg_y_idxs)

    @sync for (i, p) in enumerate(workers())
        
        # We have P processors, N workers, N blocks
        # Block ids are numbered from 1 to N
        @spawnat p makeSubOCC(occ, i, beg_y_idxs[i], sub_Nys[i])

    end
end

function run!(
    occ    :: OceanColumnCollection;
    kwargs... 
)

    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    futures= []

    @sync for (i, p) in enumerate(workers())

        @spawnat p copyfrom!(worker_occ.in_flds, master_sub_in_flds)
        @spawnat p NKOM.stepOceanColumnCollection!(
            worker_occ;
            kwargs...
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
    master_occ.h_MO[rng2...] = occ.h_MO
    master_occ.fric_u[rng2...] = occ.fric_u

    master_occ.qflx2atm[rng2...] = occ.qflx2atm


end


function sync!(
    occ :: OceanColumnCollection;
)
    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    @sync for (i, p) in enumerate(workers())
        @spawnat p NKOM.syncToMaster(worker_occ)
    end


end
