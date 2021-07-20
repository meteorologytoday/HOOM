module CyclicData

    using NCDatasets
    using CFTime

    mutable struct CyclicDataManager

        filename :: String
        beg_time :: Float64
        cyc_time :: Float64
        data     :: Dict
        t_vec    :: AbstractArray{Float64, 1}

        function CyclicDataManager(;
            filename        :: String,
            varname_time    :: String,
            varnames        :: Array,
            beg_time        :: Float64,
            cyc_time        :: Float64,
        )

            data = Dict()
            t_vec = nothing

            Dataset(filename, "r") do ds
                    
                t_vec = nomissing(ds[varname_time][:])

                t_attrib = ds[varname_time].attrib
                t_vec = [ timeencode(t_vec[i], t_attrib["units"], t_attrib["calendar"]) for i = 1:length(t_vec) ]


                for varname in varnames
                    tmp = nomissing(ds[varname][:])
                    data[varname] = reshape(ds[varname][:], :, length(t_vec))
                end

            end

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
                data,
                t_vec,
            )
        end
        
    end

    function getData!(
        cdm      :: CyclicDataManager,
        t        :: Float64,
        varnames :: Array{String, 1},
        data     :: Dict,
    )
        # lcr = left, center, right

        # Determine interpolation position
        t_c = mod(t - cdm.beg_time, cdm.cyc_time)
        
        idx_l = nothing
        idx_r = nothing

        #println("t_c = ", t_c)

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
        #println("t_l, t_c, t_r = ", t_l, "; ", t_c, ", ", t_r)

        Δt_lr = t_r - t_l
        Δt_lc = t_c - t_l 
        Δt_cr = t_r - t_c

        coe_r = Δt_lc / Δt_lr
        coe_l = Δt_cr / Δt_lr

        if any([Δt_lr, Δt_lc, Δt_cr] .< 0)
            println("[Δt_lr, Δt_lc, Δt_cr] = ", [Δt_lr, Δt_lc, Δt_cr]) 
            throw(ErrorException("Δt sign error."))
        end

        for varname in varnames
            
            tmp = data[varname]
            # interpolation happens here
            data_l = view(cdm.data[varname], :, idx_l)
            data_r = view(cdm.data[varname], :, idx_r)
            
            @. tmp = data_l * coe_l + data_r * coe_r
        end
    end 
end
