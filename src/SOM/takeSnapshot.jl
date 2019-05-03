function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds
        occ = makeBlankOceanColumnCollection(
            ds.dim["Nx"], ds.dim["Ny"];
            mask  = nomissing(ds["mask"][:], NaN),
        )


        occ.Kh_T = ds.attrib["Kh_T"]

        occ.T_ML[:, :]  = nomissing( ds["T_ML"][:], NaN )
        occ.h_ML[:, :]  = nomissing( ds["h_ML"][:], NaN )

    end

    return occ 
end


function takeSnapshot(
    occ::OceanColumnCollection,
    filename::AbstractString;
    missing_value::Float64=1e20
)
    _createNCFile(occ, filename, missing_value)

    Dataset(filename, "a") do ds
   
        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)
        _write2NCFile(ds, "T_ML", ("Nx", "Ny",), occ.T_ML, missing_value)
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
        defDim(ds, "time", Inf)
       
        ds.attrib["_FillValue"] = missing_value
        ds.attrib["Kh_T"] = occ.Kh_T
        
    end

end


function _write2NCFile(
    ds            :: Dataset,
    varname       :: AbstractString,
    dim           :: Tuple,
    var_data      :: AbstractArray{T},
    missing_value :: G ) where T where G

    ds_var = defVar(ds, varname, eltype(var_data), dim)
    ds_var.attrib["_FillValue"] = missing_value
    ds_var[:] = var_data
end


"""
    This function is meant to append 2D field into an NC file.
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

    ds_var[:, :, time] = var_data

end

