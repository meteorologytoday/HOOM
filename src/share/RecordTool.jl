module RecordTool
    using NCDatasets
    missing_value = 1e20
    
    mutable struct StatObj

        varname  :: AbstractString
        varref   :: AbstractArray{Float64}
        var       :: AbstractArray{Float64}

        dimnames :: Tuple

        weight    :: Float64

        function StatObj(varname, varref, dimnames)
            if length(size(varref)) != length(dimnames)
                ErrorException(format("Variable `{:s}` is {:d} dimensional while only {:d} dimension names are given.", varname, length(size(varref)), length(dimnames))) |> throw
            end

            var = zeros(Float64, size(varref)...)

            return new(varname, varref, var, dimnames)
        end
    
    end

    
    mutable struct Recorder

        filename :: Union{Nothing, AbstractString}
        time_ptr :: Integer  # The position of next record
        
        dims     :: Dict    # A dictionary of dimension name mapping to its length
        sobjs    :: Dict

        function Recorder(dims, vars, listfile)

            sobjs = Dict()

            if haskey(dims, "time")
                ErrorException("Dimension `time` is used for record dimension. It cannot be specified in dict `dims`.") |> throw
            end

            for (varname, varref, dimnames) in vars

                for dimname in dimnames
                    if !haskey(dims, dimname)
                        ErrorException(format("Variable `{:s}` contains a dimension `{:s}` not specified in `dims` dict.", varname, dimname)) |> throw
                    end
                end

                sobjs[varname] = StatObj(varname, varref, dimnames)
            end

            return new(nothing, 1, dims, sobjs)
        end

    end
   

    function record_wrap!(
        rec             :: Recorder;
        create_new_file :: Bool,
        avg_and_output  :: Bool,
        new_file_name   :: AbstractString,
    )

        if create_new_file
            setNewNCFile!(rec, new_file_name)
        end

        record!(rec; avg_and_output=avg_and_output)


    end

 
    function record!(
        rec::Recorder;
        avg_and_output::Bool
    )

        varnames = keys(rec.sobjs)

        for varname in varnames
            sobj = rec.sobjs[varname]
            sobj.var .+= sobj.varref
            sobj.weight += 1.0
        end

        if avg_and_output

            if rec.filename == nothing
                ErrorException("Undefined record filename") |> throw
            end
 
            # Do average
            for (varname, sobj) in rec.sobjs
                if sobj.weight == 0
                    ErrorException(format("StatObj for variable `{:s}` has weight 0 during normalization.", varname)) |> throw
                end
                sobj.var /= sobj.weight
            end
            
            # Output data
            Dataset(rec.filename, "a") do ds
                for (varname, sobj) in rec.sobjs
                    ds_var = ds[varname]
                    ds[varname][repeat([:,], length(sobj.dimnames))..., rec.time_ptr] = sobj.var
                end
            end
 

            # Reset StatObjs
            for (_, sobj) in rec.sobjs
                sobj.var .= 0.0
                sobj.weight = 0.0
            end
            
            # Increment of time
            rec.time_ptr += 1

        end
        
    end

    function setNewNCFile!(rec::Recorder, filename::AbstractString)
        rec.filename = filename
        rec.time_ptr = 1
        
        Dataset(filename, "c") do ds
            for (dimname, dim) in rec.dims 
                defDim(ds, dimname, dim)
            end

            defDim(ds, "time",   Inf)
            ds.attrib["_FillValue"] = missing_value

            for (varname, sobj) in rec.sobjs
                ds_var = defVar(ds, varname, Float64, (sobj.dimnames..., "time"))
                ds_var.attrib["_FillValue"] = missing_value
            end 
        end
        
    end

end
