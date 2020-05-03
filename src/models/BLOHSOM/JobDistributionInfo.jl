mutable struct YSplitInfo
    pull_fr_rng      :: UnitRange
    push_fr_rng      :: UnitRange
    push_to_rng      :: UnitRange
    
    pull_fr_rng_bnd  :: AbstractArray{Union{UnitRange, Nothing}}
    pull_to_rng_bnd  :: AbstractArray{Union{UnitRange, Nothing}}
    push_fr_rng_bnd  :: AbstractArray{Union{UnitRange, Nothing}}
    push_to_rng_bnd  :: AbstractArray{Union{UnitRange, Nothing}}
end

mutable struct JobDistributionInfo

    overlap :: Int64
    pids    :: AbstractArray{Integer, 1}
    
    dyn_slave_pid :: Int64
    
    # now assume tcr and mld share the same proc
    tmd_slave_pids   :: AbstractArray{Integer, 1}
    y_split_infos  :: AbstractArray{YSplitInfo, 1}

    function JobDistributionInfo(
        env     :: OcnEnv;
        overlap :: Int64 = 2,
    )
            
 
        pids  =  workers()
        dyn_slave_pid = pids[1]
        if length(pids) >= 2
            tmd_slave_pids = pids[2:end]
        elseif length(pids) == 1
            tmd_slave_pids = pids
        else
            throw(ErrorException("No available workers!"))
        end

        (

            pull_fr_rngs,
            push_fr_rngs,
            push_to_rngs,
            pull_fr_rngs_bnd,
            pull_to_rngs_bnd,
            push_fr_rngs_bnd,
            push_to_rngs_bnd 

        ) = calParallizationRange(N=env.Ny, P=length(tmd_slave_pids), L=overlap)
    
        y_split_infos = Array{YSplitInfo}(undef, length(tmd_slave_pids))
        for (i, p) in enumerate(tmd_slave_pids)

            y_split_infos[i] = YSplitInfo(
                pull_fr_rngs[i],
                push_fr_rngs[i],
                push_to_rngs[i],
                pull_fr_rngs_bnd[i, :],
                pull_to_rngs_bnd[i, :],
                push_fr_rngs_bnd[i, :],
                push_to_rngs_bnd[i, :],
            )

        end
      
        return new(
            overlap,
            pids,
            dyn_slave_pid,
            tmd_slave_pids,
            y_split_infos,
        ) 

    end

end

function calParallizationRange(;
    N = Integer,     # Total grids
    P = Integer,     # Number of procs
    L = Integer,     # Overlapping grids
)

    if ! (N >= max(1, L) * P)
        throw(ErrorException("Condition must be satisfied: N >= max(1, L) * P"))
    end

    n̄ = floor(Integer, N / P)
    R = N - n̄ * P


    pull_fr_rngs = Array{Union{UnitRange, Nothing}}(undef, P)
    push_fr_rngs = Array{Union{UnitRange, Nothing}}(undef, P)
    push_to_rngs = Array{Union{UnitRange, Nothing}}(undef, P)
    
    # 1: lower latitude side (south), 2: higher latitude side (north)
    pull_fr_rngs_bnd = Array{Union{UnitRange, Nothing}}(undef, P, 2)
    pull_to_rngs_bnd = Array{Union{UnitRange, Nothing}}(undef, P, 2)
    push_fr_rngs_bnd = Array{Union{UnitRange, Nothing}}(undef, P, 2)
    push_to_rngs_bnd = Array{Union{UnitRange, Nothing}}(undef, P, 2)



    cnt = 1
    for p = 1:P
        m = (p <= R) ? n̄ + 1 : n̄  # assigned grids


        pull_fr_rngs[p] = cnt-L:cnt+m-1+L
        push_fr_rngs[p] = L+1:L+m
        push_to_rngs[p] = cnt:cnt+m-1


        # Boundary
        pull_fr_rngs_bnd[p, 1] = cnt-L:cnt-1
        pull_fr_rngs_bnd[p, 2] = cnt+m:cnt+m+L-1

        pull_to_rngs_bnd[p, 1] = 1:L
        pull_to_rngs_bnd[p, 2] = L+m+1:L+m+L

        push_fr_rngs_bnd[p, 1] = L+1:L+L
        push_fr_rngs_bnd[p, 2] = L+m-L+1:L+m

        push_to_rngs_bnd[p, 1] = cnt:cnt+L-1
        push_to_rngs_bnd[p, 2] = cnt+m-L:cnt+m-1

        cnt += m
    end

    # South pole and north pole do not have boundaries
    pull_fr_rngs_bnd[1, 1] = nothing
    pull_to_rngs_bnd[1, 1] = nothing
    push_fr_rngs_bnd[1, 1] = nothing
    push_to_rngs_bnd[1, 1] = nothing

    pull_fr_rngs_bnd[end, 2] = nothing
    pull_to_rngs_bnd[end, 2] = nothing
    push_fr_rngs_bnd[end, 2] = nothing
    push_to_rngs_bnd[end, 2] = nothing

    # Adjust the first and last range (south pole and north pole)
    pull_fr_rngs[1] = (pull_fr_rngs[1][1]+L):pull_fr_rngs[1][end]
    pull_fr_rngs[end] = pull_fr_rngs[end][1]:(pull_fr_rngs[end][end]-L)
    push_fr_rngs[1] = 1:length(push_fr_rngs[1])


    # Change range because to southmost boundary is trimmed
    if P > 1
        pull_to_rngs_bnd[1, 2] = pull_to_rngs_bnd[1, 2] .- L
        push_fr_rngs_bnd[1, 2] = push_fr_rngs_bnd[1, 2] .- L
    end

    return pull_fr_rngs,
           push_fr_rngs,
           push_to_rngs,
           pull_fr_rngs_bnd,
           pull_to_rngs_bnd,
           push_fr_rngs_bnd,
           push_to_rngs_bnd

end

