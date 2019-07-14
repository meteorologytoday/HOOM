function loadSnapshot(
    filename::AbstractString;
    gridinfo_file::Union{AbstractString, Nothing} = nothing,
)
    local occ

    Dataset(filename, "r") do ds

        Ts_clim_relax_time = nothing
        Ss_clim_relax_time = nothing

        Ts_clim = nothing
        Ss_clim = nothing

        if haskey(ds, "Ts_clim")
            Ts_clim_relax_time = ds.attrib["Ts_clim_relax_time"]
            Ts_clim = nomissing(ds["Ts_clim"][:], NaN)
        end

        if haskey(ds, "Ss_clim")
            Ss_clim_relax_time = ds.attrib["Ss_clim_relax_time"]
            Ss_clim = nomissing(ds["Ss_clim"][:], NaN)
        end

        occ = OceanColumnCollection(
            gridinfo_file = (gridinfo_file == nothing) ? ds.attrib["gridinfo_file"] : gridinfo_file,
            Nx       = ds.dim["Nx"],
            Ny       = ds.dim["Ny"],
            zs_bone  = nomissing(ds["zs_bone"][:], NaN);
            Ts       = nomissing(ds["Ts"][:], NaN),
            Ss       = nomissing(ds["Ss"][:], NaN),
            K_T      = ds.attrib["K_T"],
            K_S      = ds.attrib["K_S"],
            fs       = nomissing(ds["fs"][:], NaN),
            ϵs       = nomissing(ds["epsilons"][:], NaN),
            T_ML     = nomissing(ds["T_ML"][:], NaN),
            S_ML     = nomissing(ds["S_ML"][:], NaN),
            h_ML     = nomissing(ds["h_ML"][:], NaN),
            h_ML_min = nomissing(ds["h_ML_min"][:], NaN),
            h_ML_max = nomissing(ds["h_ML_max"][:], NaN),
            we_max   = ds.attrib["we_max"],
            Ts_clim_relax_time = Ts_clim_relax_time,
            Ss_clim_relax_time = Ss_clim_relax_time,
            Ts_clim  = Ts_clim,
            Ss_clim  = Ss_clim,
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
 
        ds.attrib["gridinfo_file"] = occ.gi_file
        ds.attrib["K_T"] = occ.K_T
        ds.attrib["K_S"] = occ.K_S

        ds.attrib["we_max"] = occ.we_max

        if occ.Ts_clim_relax_time != nothing
            ds.attrib["Ts_clim_relax_time"] = occ.Ts_clim_relax_time
        end
 
        if occ.Ss_clim_relax_time != nothing
            ds.attrib["Ss_clim_relax_time"] = occ.Ss_clim_relax_time
        end
       
        _write2NCFile(ds, "zs_bone", ("NP_zs_bone",), occ.zs_bone, missing_value)

        _write2NCFile(ds, "Ts", ("Nx", "Ny", "Nz_bone"), occ.Ts, missing_value)
        _write2NCFile(ds, "Ss", ("Nx", "Ny", "Nz_bone"), occ.Ss, missing_value)
        _write2NCFile(ds, "bs", ("Nx", "Ny", "Nz_bone"), occ.bs, missing_value)
        _write2NCFile(ds, "T_ML", ("Nx", "Ny",), occ.T_ML, missing_value)
        _write2NCFile(ds, "S_ML", ("Nx", "Ny",), occ.S_ML, missing_value)
        _write2NCFile(ds, "h_ML", ("Nx", "Ny",), occ.h_ML, missing_value)
        
        _write2NCFile(ds, "h_ML_min", ("Nx", "Ny",), occ.h_ML_min, missing_value)
        _write2NCFile(ds, "h_ML_max", ("Nx", "Ny",), occ.h_ML_max, missing_value)

        if occ.Ts_clim != nothing
            _write2NCFile(ds, "Ts_clim", ("Nx", "Ny", "Nz_bone"), occ.Ts_clim, missing_value)
        end

        if occ.Ss_clim != nothing
            _write2NCFile(ds, "Ss_clim", ("Nx", "Ny", "Nz_bone"), occ.Ss_clim, missing_value)
        end

        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)
        _write2NCFile(ds, "topo", ("Nx", "Ny",), occ.topo, missing_value)

        _write2NCFile(ds, "fs", ("Nx", "Ny"), occ.fs, missing_value)
        _write2NCFile(ds, "epsilons", ("Nx", "Ny"), occ.ϵs, missing_value)

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
        defDim(ds, "Nz_bone", occ.Nz_bone)
        defDim(ds, "NP_zs_bone",   length(occ.zs_bone))
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

