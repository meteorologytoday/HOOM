"""
loadSnapshot
Usage: It loads a netcdf snapshot file and a config JSON file that is used in Env.
       If `timestamp` is provided, the function will check if it matches the 
       timestamp recorded in the netcdf file.
"""
function loadSnapshot(
    filename_nc   :: AbstractString,   # Field
    filename_cfg  :: AbstractString;   # Configuration
    timestamp     :: Union{AbstractCFDateTime, Nothing} = nothing,
)

    writeLog("Reading files: {:s} and {:s}", filename_nc, filename_cfg)
    
    if timestamp != nothing
        timestamp_str = Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS")
    end

    Dataset(filename_nc, "r") do ds

        if timestamp_str != nothing && timestamp_str != ds.attrib["timestamp"]
            throw(ErrorException(format
                "The provided timestamp is {:s}, but the netcdf file timestamp is {:s}",
                timestamp_str,
                ds.attrib["timestamp"], 
            ))
        end

        for (varname, data_unit) in data_table.data_units # HOOM.getDynamicVariableList(mb; varsets = [:ALL,])

            println("Reading ", varname, "... ")
            data_unit.sdata2[:, :, :] = nomissing(ds[varname][:, :, :, 1], NaN)

        end

    end

    open(filename_nc, "r") do ds
        cfg = JSON.json(filename_cfg)
    end

    return 
end

"""
takeSnapshot

Usage: It output netcdf file that is a copy of variables in the given DataTable
       It also converts a config Dict into JSON text file
"""
function takeSnapshot(
    timestamp     :: AbstractCFDateTime,
    mb            :: HOOM.ModelBlock,
    filename_nc   :: AbstractString,   # Field
    filename_cfg  :: AbstractString;   # Configuration
    missing_value :: Float64=1e20,
)

    timestamp_str = Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS")
    data_table = mb.dt

    Dataset(filename_nc, "c") do ds

        defDim(ds, "time",   Inf)
        for (dimname, dimvalue) in data_table.dims
            defDim(ds, dimname, dimvalue)
        end

        ds.attrib["_FillValue"] = missing_value
        ds.attrib["timestamp"] = timestamp_str

        varnames = keys(HOOM.getDynamicVariableList(mb; varsets = [:ALL,]))
        for (varname, data_unit) in data_table.data_units

            data_unit = data_table.data_units[varname]

            println("Writing ", varname, "... ")
            ds_var = defVar(ds, varname, eltype(data_unit.sdata2), (data_table.grid_dims2_str[data_unit.grid]..., "time"))
            ds_var.attrib["_FillValue"] = missing_value 
            ds_var = ds[varname]

            ds_var[:, :, :, 1] = data_unit.sdata2

        end

    end

    open(filename_cfg, "w") do io
        JSON.print(io, mb.ev.config, 2)
    end

    
end


