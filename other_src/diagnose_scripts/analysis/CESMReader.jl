module CESMReader

    using NCDatasets
    using Formatting
    export FileHandler, getData  

    include("nanop.jl")

    mutable struct FileHandler
        filename_format :: String
        form   :: Symbol

        function FileHandler(;
            filename_format :: String,
            form        :: Symbol = :YEAR_MONTH,
        )

            if ! (form in [:YEAR, :YEAR_MONTH]  )
                throw(ErrorException("Error: only two forms allowed. Either :YEAR or :YEAR_MONTH. Now we got " * string(form)))
            end

            return new(
                filename_format,
                form,
            )
        end

    end

    function getData(
        handler       :: FileHandler,
        varname  :: String,
        year_rng :: Union{Tuple, Array},
        idx...;
        verbose=true,
    )
        return getData(handler, varname, (year_rng[1], 1), (year_rng[2], 12), idx...; verbose=verbose)
    end

    function getData(
        handler       :: FileHandler,
        varname  :: String,
        beg_time :: Union{Tuple, Array},
        end_time :: Union{Tuple, Array},
        idx...;
        verbose=true,
    )

        beg_y, beg_m = beg_time
        end_y, end_m = end_time

        beg_t = (beg_y - 1) * 12 + beg_m - 1
        end_t = (end_y - 1) * 12 + end_m - 1

        months = end_t - beg_t + 1
        
        if months <= 0
            throw(ErrorException("End time should be larger than begin time"))
        end

        if ! ( ( 1 <= beg_m <= 12 ) && ( 1 <= end_m <= 12 ) )
            throw(ErrorException("Invalid month"))
        end

        if length(idx) == 0
            throw(ErrorException("Spatial range not specified."))
        end

        local data, new_idx

        if handler.form == :YEAR_MONTH

            flag_1 = true

            for y in beg_y:end_y, m in beg_m:end_m

                current_t = (y-1) * 12 + (m-1) - beg_t + 1

                filename = format(handler.filename_format, y, m)
                verbose && println("Loading file: ", filename)
                ds = Dataset(filename, "r")

                partial_data = nomissing(ds[varname][idx..., 1])
                if flag_1
                    new_size = size(partial_data)
                    new_idx  = [Colon() for _ in 1:length(new_size)]
                    data = zeros(Float64, new_size..., months)
                    data .= NaN
                    flag_1 = false
                end

                data[new_idx..., current_t] = partial_data

                close(ds)

            end            
            
        elseif handler.form == :YEAR
     
            flag_1 = true

            for y in beg_y:end_y

                current_t = (y-1) * 12 - beg_t
                
                filename = format(handler.filename_format, y)
                verbose && println("Loading file: ", filename)
                ds = Dataset(filename, "r")

                if y == beg_y
                    rng = (idx..., beg_m:12)
                elseif y == end_y
                    rng = (idx..., 1:end_m)
                else
                    rng = (idx..., 1:12)
                end  

                partial_data = nomissing(ds[varname][rng...])

                if flag_1
                    new_size = size(partial_data)[1:end-1]
                    new_idx  = [Colon() for _ in 1:length(new_size)]
                    data = zeros(Float64, new_size..., months)
                    flag_1 = false
                end



                if y == beg_y
                    rng = beg_m:12
                elseif y == end_y
                    rng = 1:end_m
                else
                    rng = 1:12
                end
                
                offset = (y-beg_y) * 12
                data[new_idx..., rng .+ offset] = partial_data
                  

                close(ds)

            end

        end 

        return data
    end

end
