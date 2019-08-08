mutable struct SubOcean
    mission           :: Symbol
    rng2              :: AbstractArray
    rng3              :: AbstractArray
    master_in_flds    :: InputFields
    master_ocn        :: Ocean
    worker_ocn        :: Ocean
end

function makeSubOcean(
    master_ocn   :: Ocean,
    mission      :: Symbol,
    block_id     :: Integer,
    beg_idx      :: Integer,
    sub_N        :: Integer,
)

    global rng2, rng3

    if mission == :vt
        rng2 = [Colon(), beg_idx:beg_idx+sub_N-1]
        rng3 = [Colon(), Colon(), beg_idx:beg_idx+sub_N-1]
        Ny = sub_N
    elseif mission == :hz
        rng2 = [Colon(), Colon()]
        rng3 = [beg_idx:beg_idx+sub_N-1, Colon(), Colon()]
        Ny = master_ocn.Ny
    else
        throw(ErrorException("[makeSubOcean] Unknown symbol for `kind`: ", string(kind)))
    end

    return SubOcean(
        mission,
        rng2,
        rng3,
        SubInputFields(master_ocn.in_flds, rng2...), 
        master_ocn,
        Ocean(
            id             = block_id,
            gridinfo_file  = nothing,
            Nx             = master_ocn.Nx,
            Ny             = Ny,
            zs_bone        = master_ocn.zs_bone,
            Ts             = master_ocn.Ts[rng3...],
            Ss             = master_ocn.Ss[rng3...],
            K_T            = master_ocn.K_T,
            K_S            = master_ocn.K_S,
            T_ML           = master_ocn.T_ML[rng2...],
            S_ML           = master_ocn.S_ML[rng2...],
            h_ML           = master_ocn.h_ML[rng2...],
            h_ML_min       = master_ocn.h_ML_min[rng2...],
            h_ML_max       = master_ocn.h_ML_max[rng2...],
            we_max         = master_ocn.we_max,
            R              = master_ocn.R,
            ζ              = master_ocn.ζ,
            Ts_clim_relax_time = master_ocn.Ts_clim_relax_time,
            Ss_clim_relax_time = master_ocn.Ss_clim_relax_time,
            Ts_clim        = ( master_ocn.Ts_clim != nothing ) ? master_ocn.Ts_clim[rng3...] : nothing,
            Ss_clim        = ( master_ocn.Ss_clim != nothing ) ? master_ocn.Ss_clim[rng3...] : nothing,
            mask           = master_ocn.mask[rng2...],
            topo           = master_ocn.topo[rng2...],
            fs             = master_ocn.fs[rng2...],
            ϵs             = master_ocn.ϵs[rng2...],
            in_flds        = InputFields(:local, master_ocn.Nx, Ny),
            arrange        = :zxy,
        )
    )

end

function syncToMaster(subocn::SubOcean)

#    (subocn.woocc.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
    
    master_ocn.Ts[rng3...] = worker_ocn.Ts
    master_ocn.Ss[rng3...] = worker_ocn.Ss

    master_ocn.FLDO[rng2...] = worker_ocn.FLDO
    master_ocn.T_ML[rng2...] = worker_ocn.T_ML 
    master_ocn.S_ML[rng2...] = worker_ocn.S_ML
    master_ocn.h_ML[rng2...] = worker_ocn.h_ML
    master_ocn.h_MO[rng2...] = worker_ocn.h_MO
    master_ocn.fric_u[rng2...] = worker_ocn.fric_u

    master_ocn.qflx2atm[rng2...] = worker_ocn.qflx2atm


end

function syncFromMaster!(subocn::SubOcean)

#    (subocn.woocc.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
    
    # View is to avoid array allocation
    worker_ocn.Ts[:] = view( master_ocn.Ts, rng3...)
    worker_ocn.Ss[:] = view( master_ocn.Ss, rng3...)
    
    worker_ocn.FLDO[:] = view( master_ocn.FLDO, rng2...)
    worker_ocn.T_ML[:] = view( master_ocn.T_ML, rng2...)
    worker_ocn.S_ML[:] = view( master_ocn.S_ML, rng2...)
    worker_ocn.h_ML[:] = view( master_ocn.h_ML, rng2...)

    copyfrom!(worker_ocn.in_flds, subocn.master_in_flds)
end



function init(ocn::Ocean)

    println("Number of workers: ", nworkers())

    (ocn.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(ocn.id)))


    # Sub ocean cols
    sub_Ny = ceil(Integer, ocn.Nx / nworkers())
    sub_Nys = [sub_Ny for block_id = 1:nworkers()]
    sub_Nys[end] = ocn.Ny - (nworkers()-1) * sub_Ny
    beg_y_idxs = [sub_Ny * (block_id - 1) + 1 for block_id = 1:nworkers()]

    # Sub ocean lays
    sub_Nz = ceil(Integer, ocn.Nz_bone / nworkers())
    sub_Nzs = [sub_Nz for block_id = 1:nworkers()]
    sub_Nzs[end] = ocn.Nz_bone - (nworkers()-1) * sub_Nz
    beg_z_idxs = [sub_Nz * (block_id - 1) + 1 for block_id = 1:nworkers()]

    @sync for (i, p) in enumerate(workers())
        
        # We have P processors, N workers, N blocks
        # Block ids are numbered from 1 to N
        @spawnat p let
            global subocn_hz = makeSubOcean(ocn, :hz, i, beg_z_idxs[i], sub_Nzs[i])
            global subocn_vt = makeSubOcean(ocn, :vt, i, beg_y_idxs[i], sub_Nys[i])
        end

    end
end

function run!(
    ocn :: Ocean;
    cfgs...
)

    (ocn.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(ocn.id)))

    @sync for (i, p) in enumerate(workers())

        @spawnat p let

            syncFromMaster!(subocn_hz)
            NKOM.stepOcean_hz!(subocn_hz.worker_ocn; cfgs...)
            syncToMaster(subocn_hz)

        end

    end

    @sync for (i, p) in enumerate(workers())

        @spawnat p let

            syncFromMaster!(subocn_vt)
            NKOM.stepOcean_vt!(subocn_vt.worker_ocn; cfgs...)
            syncToMaster(subocn_vt)

        end
    end

end

#=
function sync!(
    occ :: OceanColumnCollection;
)
    (occ.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(occ.id)))

    @sync for (i, p) in enumerate(workers())
        @spawnat p NKOM.syncToMaster(worker_occ)
    end


end
=#




