function loadSnapshot(
    filename :: AbstractString;
    gridinfo_file :: Union{Nothing, AbstractString} = nothing,
)

    local occ

    Dataset(filename, "r") do ds

        occ = OceanColumnCollection(;
            gridinfo_file = (gridinfo_file == nothing) ? ds.attrib["gridinfo_file"] : gridinfo_file,
            Nx       = ds.dim["Nx"],
            Ny       = ds.dim["Ny"],
            hs       = nomissing(ds["hs"][:], NaN),
            Ts       = nomissing(ds["Ts"][:], NaN),
            Ss       = nomissing(ds["Ss"][:], NaN),
            Kh_T     = ds.attrib["Kh_T"],
            Kh_S     = ds.attrib["Kh_S"],
            fs       = nomissing(ds["fs"][:], NaN),
            ϵs       = nomissing(ds["epsilons"][:], NaN),
            mask     = nomissing(ds["mask"][:], NaN),
            topo     = nomissing(ds["topo"][:], NaN),
        )

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
        ds.attrib["Kh_T"] = occ.Kh_T
        ds.attrib["Kh_S"] = occ.Kh_S
        ds.attrib["gridinfo_file"] = occ.gi_file
  
        _write2NCFile(ds, "hs", ("Nz",), occ.hs, missing_value)
        _write2NCFile(ds, "Ts", ("Nx", "Ny", "Nz"), occ.Ts, missing_value)
        _write2NCFile(ds, "Ss", ("Nx", "Ny", "Nz"), occ.Ss, missing_value)
        _write2NCFile(ds, "fs", ("Nx", "Ny"), occ.fs, missing_value)
        _write2NCFile(ds, "epsilons", ("Nx", "Ny"), occ.ϵs, missing_value)
        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)
        _write2NCFile(ds, "topo", ("Nx", "Ny",), occ.topo, missing_value)

    end

end

function _createNCFile(
    occ::OceanColumnCollection,
    filename::AbstractString,
    missing_value::Float64,
)

    Dataset(filename, "c") do ds

        defDim(ds, "Nx", occ.Nx)
        defDim(ds, "Ny", occ.Ny)
        defDim(ds, "Nz", 2)
        defDim(ds, "time",   Inf)
       
        ds.attrib["_FillValue"] = missing_value

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

