module CyclicData

    using NCDatasets
    using CFTime
    using Formatting

    export CyclicDataManager, interpData!
 
    mutable struct CyclicDataManager

        filename  :: String

        beg_time  :: Float64
        cyc_time  :: Float64
        sub_yrng  :: Union{Colon, UnitRange}
        varnames  :: Array{String, 1}

        data_l   :: Union{Dict, Nothing}     # Data on the lower point. Used for interpolation
        data_r   :: Union{Dict, Nothing}     # Data on the upper point. Used for interpolation

        t_vec     :: AbstractArray{Float64, 1}

        # Time pointer. Indicating the last data read. If it is zero then no previous data is read
        t_ptr_l   :: Integer
        t_ptr_r   :: Integer

        function CyclicDataManager(;
            filename        :: String,
            varname_time    :: String,
            varnames        :: Array,
            beg_time        :: Float64,
            cyc_time        :: Float64,
            sub_yrng        :: Union{UnitRange, Colon} = Colon(),
        )

            data = Dict()
            t_vec = nothing

            Dataset(filename, "r") do ds
                    
                t_vec = nomissing(ds[varname_time][:])

                t_attrib = ds[varname_time].attrib
                t_vec = [ timeencode(t_vec[i], t_attrib["units"], t_attrib["calendar"]) for i = 1:length(t_vec) ]

            end

            println("t_vec: ", t_vec)
            #=

                for varname in varnames
                    var = ds[varname]
                    if length(size(var)) == 4  # 3D 
                        tmp = nomissing(ds[varname][spatial_rng..., :])
                        if xyz2zxy 
                            tmp = toXYZ(tmp, :zxy)
                        end
                    data[varname] = copy(reshape(tmp, :, length(t_vec)))
                end

            end
            =#

            #println(typeof(t_vec))

            if any( (t_vec[2:end] - t_vec[1:end-1]) .<= 0.0 )
                throw(ErrorException("Time dimension has to be monotonically increasing."))
            end

            if cyc_time <= 0.0
                throw(ErrorException("End time cannot be earlier than beg time"))
            end

            if t_vec[end] - t_vec[1] > cyc_time
               throw(ErrorException("Time has to be within one cycle time.")) 
            end

            return new(
                filename,
                beg_time,
                cyc_time,
                sub_yrng,
                varnames,
                nothing,
                nothing,
                t_vec,
                0,
                0,
            )
        end
        
    end


    function detectTimeBoundary(
        cdm :: CyclicDataManager,
        t   :: Float64,
    )

        # lcr = left, center, right
        # Determine interpolation position
        t_c = mod(t - cdm.beg_time, cdm.cyc_time) + cdm.beg_time
       
        if t_c >= cdm.t_vec[end]
            idx_l = length(cdm.t_vec)
            idx_r = 1
            t_l = cdm.t_vec[idx_l]
            t_r = cdm.t_vec[idx_r] + cdm.cyc_time
        elseif t_c < cdm.t_vec[1]
            idx_l = length(cdm.t_vec)
            idx_r = 1
            t_l = cdm.t_vec[idx_l]
            t_r = cdm.t_vec[idx_r] + cdm.cyc_time
            t_c += cdm.cyc_time
        else
            for i = 1:length(cdm.t_vec)-1
                if cdm.t_vec[i] <= t_c < cdm.t_vec[i+1]
                    idx_l = i
                    idx_r = i+1
                    t_l = cdm.t_vec[idx_l]
                    t_r = cdm.t_vec[idx_r]
                    break
                end
            end
        end

        return idx_l, idx_r, t_l, t_c, t_r 

    end

    function loadData(
        cdm      :: CyclicDataManager,
        t_idx    :: Int64,
    )

        local data = Dict()

        Dataset(cdm.filename, "r") do ds

            for varname in cdm.varnames
                var = ds[varname]
                s = size(var)

                if length(s) == 4  # 3D case
                    data[varname] = permutedims(nomissing(ds[varname][:, cdm.sub_yrng, :, t_idx], NaN), [3,1,2]) 
                elseif length(s) == 3  # 2D case
                    data[varname] = reshape(nomissing(ds[varname][cdm.sub_yrng, :, t_idx], NaN), 1, s[1:2]...)
                else
                    throw(ErrorException("Unknown dimension: " * string(s)))
                end

            end

        end
        
        return data
    end


    function interpData!(
        cdm      :: CyclicDataManager,
        t        :: Float64,
        data     :: Union{Dict, Nothing} = nothing;
        create   :: Bool = false,
    )

        if data == nothing
            data = Dict()
        end

        idx_l, idx_r, t_l, t_c, t_r = detectTimeBoundary(cdm, t)

        #println(format("{:d}, {:d}, {:.1f}", idx_l, idx_r, t_c))

        if cdm.t_ptr_l == 0 && cdm.t_ptr_r == 0 # Initialization

            #println("Initialization")
            # load idx_l into data_l
            # load idx_r into data_r

            cdm.t_ptr_l, cdm.t_ptr_r = idx_l, idx_r

            cdm.data_l = loadData(cdm, cdm.t_ptr_l)
            cdm.data_r = loadData(cdm, cdm.t_ptr_r)

        elseif cdm.t_ptr_r == idx_l

            #println("Move on!")
            # move data_r to data_l
            cdm.t_ptr_l, cdm.t_ptr_r = cdm.t_ptr_r, idx_r
            cdm.data_l = cdm.data_r
    
            # load idx_r into data_r
            cdm.data_r = loadData(cdm, cdm.t_ptr_r)

        elseif cdm.t_ptr_l == idx_l && cdm.t_ptr_r == idx_r

            # do nothing. 

        else
            throw(ErrorException(format("Unknown situation. cdm.t_ptr_l = {:d}, cdm.t_ptr_r = {:d}", cdm.t_ptr_l, cdm.t_ptr_r)))
        end


        Δt_lr = t_r - t_l
        Δt_lc = t_c - t_l 
        Δt_cr = t_r - t_c

#        println(Δt_lr, "; ", Δt_lc, "; ", Δt_cr)

        coe_r = Δt_lc / Δt_lr
        coe_l = Δt_cr / Δt_lr
 
        #println(format("{:.2f}, {:.2f}", coe_l, coe_r))

        if any([Δt_lr, Δt_lc, Δt_cr] .< 0)
            println("[Δt_lr, Δt_lc, Δt_cr] = ", [Δt_lr, Δt_lc, Δt_cr]) 
            throw(ErrorException("Δt sign error."))
        end

        for varname in cdm.varnames
    
            if ! haskey(data, varname)
                if create
                    data[varname] = zeros(Float64, size(cdm.data_l[varname])...)
                else
                    continue
                end
            end

            tmp = view(data[varname], :)

            # interpolation happens here
            data_l = view(cdm.data_l[varname], :)
            data_r = view(cdm.data_r[varname], :)
            #println(size(data_l), "; ", size(data_r))            
            @. tmp = data_l * coe_l + data_r * coe_r
        end

        return data
    end 
end
