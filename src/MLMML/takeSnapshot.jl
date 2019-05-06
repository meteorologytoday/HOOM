function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds
        occ = makeBlankOceanColumnCollection(
            ds.dim["Nx"],
            ds.dim["Ny"],
            nomissing(ds["zs_bone"][:], NaN);
            mask = nomissing(ds["mask"][:], NaN),
            topo = nomissing(ds["topo"][:], NaN),
        )

        occ.K_T = ds.attrib["K_T"]
        occ.K_S = ds.attrib["K_S"]
        occ.h_ML_min = ds.attrib["h_ML_min"]
        occ.h_ML_max = ds.attrib["h_ML_max"]
        occ.we_max = ds.attrib["we_max"]

        occ.h_ML[:, :]  = nomissing( ds["h_ML"][:], NaN )

        occ.Ts[:, :, :] = nomissing( ds["Ts"][:], NaN )
        occ.Ss[:, :, :] = nomissing( ds["Ss"][:], NaN )

        occ.T_ML[:, :]  = nomissing( ds["T_ML"][:], NaN )
        occ.S_ML[:, :]  = nomissing( ds["S_ML"][:], NaN )


    end

    updateFLDO!(occ)
    updateB!(occ)

    return occ 
end


function takeSnapshot(
    occ::OceanColumnCollection,
    filename::AbstractString;
    missing_value::Float64=1e20
)
    _createNCFile(occ, filename, missing_value)

    Dataset(filename, "a") do ds
   
        _write2NCFile(ds, "zs", ("N_zs",), occ.zs, missing_value)
        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)

        _write2NCFile(ds, "Ts", ("Nx", "Ny", "Nz"), occ.Ts, missing_value)
        _write2NCFile(ds, "Ss", ("Nx", "Ny", "Nz"), occ.Ss, missing_value)
        _write2NCFile(ds, "T_ML", ("Nx", "Ny",), occ.T_ML, missing_value)
        _write2NCFile(ds, "S_ML", ("Nx", "Ny",), occ.S_ML, missing_value)
        _write2NCFile(ds, "h_ML", ("Nx", "Ny",), occ.h_ML, missing_value)
    end

end

function _createNCFile(
    occ::OceanColumnCollection,
    filename::AbstractString,
    missing_value::Float64,
)

    Dataset(filename, "c") do ds

        defDim(ds, "N_ocs", occ.N_ocs)
        defDim(ds, "Nx", occ.Nx)
        defDim(ds, "Ny", occ.Ny)
        defDim(ds, "Nz", occ.Nz)
        defDim(ds, "N_zs",   length(occ.zs))
        defDim(ds, "time",   Inf)
       
        ds.attrib["_FillValue"] = missing_value
        ds.attrib["K_T"] = occ.K_T
        ds.attrib["K_S"] = occ.K_S
        ds.attrib["h_ML_min"] = occ.h_ML_min
        ds.attrib["h_ML_max"] = occ.h_ML_max
        ds.attrib["we_max"] = occ.we_max
        
    end

end


function _write2NCFile(
    ds            :: Dataset,
    varname       :: AbstractString,
    dim           :: Tuple,
    var_data      :: AbstractArray{T},
    missing_value :: G) where T where G

    #println("Write : ", varname)

    ds_var = defVar(ds, varname, eltype(var_data), dim)
    ds_var.attrib["_FillValue"] = missing_value
    ds_var[:] = var_data
end

"""
    This function is meant to append field into an NC file
    along the time dimension
"""
function _write2NCFile_time(
    ds            :: Dataset,
    varname       :: String,
    dim           :: Tuple,
    time          :: Integer,
    var_data      :: AbstractArray{T};
    missing_value :: Union{T, Nothing} = nothing,
) where T where G

#    time        :: Union{Nothing, UnitRange, Integer} = nothing,
#    time_exists :: Bool = true,
#    missing_value :: Union{T, Nothing} = nothing,
#) where T <: float

    local ds_var

    # Create variable if it is not in the file yet
    if ! ( varname in keys(ds) )

        ds_var = defVar(ds, varname, T, (dim..., "time"))
        
        if missing_value != nothing
            ds_var.attrib["_FillValue"] = missing_value 
        end
    else
        ds_var = ds[varname]
    end

    
    ds_var[repeat([:,], length(dim))..., time] = var_data

end

