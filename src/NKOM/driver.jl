mutable struct SubOcean
    master_in_flds    :: InputFields
    master_ocn        :: Ocean
    worker_ocn        :: Ocean
    block_id          :: Integer

    pull_fr_rng2      :: AbstractArray
    pull_fr_rng3      :: AbstractArray

    push_fr_rng2      :: AbstractArray
    push_fr_rng3      :: AbstractArray
    push_to_rng2      :: AbstractArray
    push_to_rng3      :: AbstractArray

end

function makeSubOcean(
    master_ocn   :: Ocean,
    block_id     :: Integer,
    nblocks      :: Integer,
)

    println(format("{:03d} Entering makeSubOcean.", block_id))

    touch_southpole = block_id == 1
    touch_northpole = block_id == nblocks

    sub_Ny_wo_ghost = ceil(Integer, master_ocn.Ny / nblocks)

    if block_id != nblocks
        sub_Ny = sub_Ny_wo_ghost
    else
        sub_Ny = master_ocn.Ny - (nblocks-1) * sub_Ny_wo_ghost
    end

    if sub_Ny <= 0
        throw(ErrorException("sub_Ny <= 0. Please check your resolution and nblocks"))            
    end

    pull_fr_beg_y = push_to_beg_y = sub_Ny_wo_ghost * (block_id - 1) + 1
    pull_fr_end_y = push_to_end_y = pull_fr_beg_y + sub_Ny - 1
    
    push_fr_beg_y = 1
    push_fr_end_y = sub_Ny

    if ! touch_southpole
        # expand south boundary
        pull_fr_beg_y -= 1
        sub_Ny += 1
    
        # fix push from range.
        # We want to skip the expanded latitude
        push_fr_beg_y += 1
        push_fr_end_y += 1
    end

    if ! touch_northpole
        # expand north boundary
        pull_fr_end_y += 1
        sub_Ny += 1
        
        # No need to fix push from range.
        # It is not affected.
    end


    pull_fr_rng2 = [Colon(),          pull_fr_beg_y:pull_fr_end_y]
    pull_fr_rng3 = [Colon(), Colon(), pull_fr_beg_y:pull_fr_end_y]

    push_to_rng2 = [Colon(),          push_to_beg_y:push_to_end_y]
    push_to_rng3 = [Colon(), Colon(), push_to_beg_y:push_to_end_y]

    push_fr_rng2 = [Colon(),          push_fr_beg_y:push_fr_end_y]
    push_fr_rng3 = [Colon(), Colon(), push_fr_beg_y:push_fr_end_y]


    println("pull_fr_rng2: ", pull_fr_rng2)
    println("push_to_rng2: ", push_to_rng2)
    println("push_fr_rng2: ", push_fr_rng2)

    #=
    println("### rng3: ")
    println("pull_fr_rng3: ", pull_fr_rng3)
    println("push_to_rng3: ", push_to_rng3)
    println("push_fr_rng3: ", push_fr_rng3)
    =#


    if length(pull_fr_rng2[2]) != sub_Ny
        throw(ErrorException("Pull from dimension does not match sub_Ny"))
    end

    if length(push_fr_rng2[2]) != length(push_to_rng2[2])
        throw(ErrorException("Push from and push to dimensions do not match."))
    end


    return SubOcean(
        SubInputFields(master_ocn.in_flds, pull_fr_rng2...), 
        master_ocn,
        Ocean(
            id             = block_id,
            gridinfo_file  = master_ocn.gi_file,
            sub_yrng       = pull_fr_beg_y:pull_fr_end_y,
            Nx             = master_ocn.Nx,
            Ny             = sub_Ny,
            zs_bone        = master_ocn.zs_bone,
            Ts             = master_ocn.Ts[pull_fr_rng3...],
            Ss             = master_ocn.Ss[pull_fr_rng3...],
            K_T            = master_ocn.K_T,
            K_S            = master_ocn.K_S,
            T_ML           = master_ocn.T_ML[pull_fr_rng2...],
            S_ML           = master_ocn.S_ML[pull_fr_rng2...],
            h_ML           = master_ocn.h_ML[pull_fr_rng2...],
            h_ML_min       = master_ocn.h_ML_min[pull_fr_rng2...],
            h_ML_max       = master_ocn.h_ML_max[pull_fr_rng2...],
            we_max         = master_ocn.we_max,
            R              = master_ocn.R,
            ζ              = master_ocn.ζ,
            Ts_clim_relax_time = master_ocn.Ts_clim_relax_time,
            Ss_clim_relax_time = master_ocn.Ss_clim_relax_time,
            Ts_clim        = ( master_ocn.Ts_clim != nothing ) ? master_ocn.Ts_clim[pull_fr_rng3...] : nothing,
            Ss_clim        = ( master_ocn.Ss_clim != nothing ) ? master_ocn.Ss_clim[pull_fr_rng3...] : nothing,
            mask           = master_ocn.mask[pull_fr_rng2...],
            topo           = master_ocn.topo[pull_fr_rng2...],
            fs             = master_ocn.fs[pull_fr_rng2...],
            ϵs             = master_ocn.ϵs[pull_fr_rng2...],
            in_flds        = InputFields(:local, master_ocn.Nx, sub_Ny),
            arrange        = :zxy,
        ),
        block_id,
        pull_fr_rng2,
        pull_fr_rng3,
        push_fr_rng2,
        push_fr_rng3,
        push_to_rng2,
        push_to_rng3,
    )

end

function syncToMaster!(subocn::SubOcean;
        vars2 :: Tuple = (),
        vars3 :: Tuple = (),
)

    (subocn.worker_ocn.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
   
    w_rng3 = subocn.push_fr_rng3
    w_rng2 = subocn.push_fr_rng2
 
    m_rng3 = subocn.push_to_rng3
    m_rng2 = subocn.push_to_rng2
 
    for var in vars2
        getfield(master_ocn, var)[m_rng2...] = view(getfield(worker_ocn, var), w_rng2...)
    end

    for var in vars3
        getfield(master_ocn, var)[m_rng3...] = view(getfield(worker_ocn, var), w_rng3...)
    end

end

function syncFromMaster!(
    subocn::SubOcean
)

    (subocn.worker_ocn.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
 
    copyfrom!(worker_ocn.in_flds, subocn.master_in_flds)
end



function syncBoundaryToMaster!(subocn::SubOcean;
        vars2 :: Tuple = (),
        vars3 :: Tuple = (),
)

    (subocn.worker_ocn.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
  
    w_rng3 = subocn.push_fr_rng3
    w_rng2 = subocn.push_fr_rng2
 
    m_rng3 = subocn.push_to_rng3
    m_rng2 = subocn.push_to_rng2
 
    for var in vars2
        getfield(master_ocn, var)[:, m_rng2[2][  1]] = view(getfield(worker_ocn, var), :, w_rng2[2][  1])
        getfield(master_ocn, var)[:, m_rng2[2][end]] = view(getfield(worker_ocn, var), :, w_rng2[2][end])
    end

    for var in vars3
        getfield(master_ocn, var)[:, :, m_rng3[3][  1]] = view(getfield(worker_ocn, var), :, :, w_rng3[3][  1])
        getfield(master_ocn, var)[:, :, m_rng3[3][end]] = view(getfield(worker_ocn, var), :, :, w_rng3[3][end])
    end

end

function syncBoundaryFromMaster!(subocn::SubOcean;
        vars2 :: Tuple = (),
        vars3 :: Tuple = (),
)

    (subocn.worker_ocn.id == 0) && throw(ErrorException("`id` should not be 0 (master)."))

    master_ocn = subocn.master_ocn
    worker_ocn = subocn.worker_ocn
  
    m_rng3 = subocn.pull_fr_rng3
    m_rng2 = subocn.pull_fr_rng2
 
    for var in vars2
        getfield(worker_ocn, var)[:,   1] = view(getfield(master_ocn, var), :, m_rng2[2][  1])
        getfield(worker_ocn, var)[:, end] = view(getfield(master_ocn, var), :, m_rng2[2][end])
    end

    for var in vars3
        getfield(worker_ocn, var)[:, :,   1] = view(getfield(master_ocn, var), :, :, m_rng3[3][  1])
        getfield(worker_ocn, var)[:, :, end] = view(getfield(master_ocn, var), :, :, m_rng3[3][end])
    end

end





function init(ocn::Ocean)

    global  wkrs  =  workers()
    nwkrs = length(wkrs) 

    println("Number of all workers: ", length(wkrs))

    (ocn.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(ocn.id)))

    tmp_cols = ocn.cols
    tmp_lays = ocn.lays

    ocn.cols = nothing
    ocn.lays = nothing

    @sync for (i, p) in enumerate(wkrs)
            # We have P processors, N workers, N blocks
            # Block ids are numbered from 1 to N
            @spawnat p let
                global subocn = makeSubOcean(ocn, i, nwkrs)
            end

    end

    ocn.cols = ocn.cols
    ocn.lays = ocn.lays

end

function run!(
    ocn :: Ocean;
    Δt  :: Float64,
    substeps :: Integer,
    cfgs...
)

    (ocn.id == 0) || throw(ErrorException("`id` is not 0 (master). Id received: " * string(ocn.id)))

    dt = Δt / substeps

    @sync for (i, p) in enumerate(wkrs)
        @spawnat p let
            syncFromMaster!(subocn)
            cleanQflx2atm!(subocn.worker_ocn)
            stepOcean_prepare!(subocn.worker_ocn; cfgs...)
        end
    end


    sync_bnd_vars2 = (:T_ML, :S_ML, :h_ML, :FLDO)
    sync_bnd_vars3 = (:Ts,   :Ss)

    sync_to_master_vars2 = (:FLDO, :T_ML, :S_ML, :h_ML, :h_MO, :fric_u, :qflx2atm, :τx, :τy, :dTdt_ent, :dSdt_ent)
    sync_to_master_vars3 = (:Ts, :Ss, :bs, :u, :v, :w, :T_hadvs, :T_vadvs, :S_hadvs, :S_vadvs)

    #accumulative_vars2 = (:dTdt_ent, :dSdt_ent)
    #accumulative_vars3 = (:T_hadvs, :T_vadvs, :S_hadvs, :S_vadvs)

    for substep = 1:substeps

        @sync for (i, p) in enumerate(wkrs)
            @spawnat p let
                syncBoundaryFromMaster!(subocn; vars3 = sync_bnd_vars3, vars2 = sync_bnd_vars2)
                stepOcean_Flow!(subocn.worker_ocn; Δt = dt, cfgs...)
                stepOcean_MLDynamics!(subocn.worker_ocn; Δt = dt, cfgs...)
                syncBoundaryToMaster!(subocn; vars3 = sync_bnd_vars3, vars2 = sync_bnd_vars2)
                accumulate!(subocn.worker_ocn)
            end
        end

    end

    @sync for (i, p) in enumerate(wkrs)
        @spawnat p let
            stepOcean_slowprocesses!(subocn.worker_ocn; Δt = Δt, cfgs...)
            calQflx2atm!(subocn.worker_ocn; Δt=Δt)
            avg_accumulate!(subocn.worker_ocn; count=substeps)
            syncToMaster!(
                subocn;
                vars2 = sync_to_master_vars2,
                vars3 = sync_to_master_vars3,
            )
        end
    end

end

