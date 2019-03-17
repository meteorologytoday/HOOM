function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds
        occ = makeBlankOceanColumnCollection(
            ds.dim["Nx"], ds.dim["Ny"], nomissing(ds["zs"][:], NaN);
            mask  = nomissing(ds["mask"][:], NaN),
        )


        occ.K_T = ds.attrib["K_T"]
        occ.K_S = ds.attrib["K_S"]

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

        _write2NCFile(ds, "Ts", ("Nx", "Ny", "N_lays"), occ.Ts, missing_value)
        _write2NCFile(ds, "Ss", ("Nx", "Ny", "N_lays"), occ.Ss, missing_value)
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
        defDim(ds, "N_lays", occ.Nz)
        defDim(ds, "N_zs",   length(occ.zs))
       
        ds.attrib["_FillValue"] = missing_value
        ds.attrib["K_T"] = occ.K_T
        ds.attrib["K_S"] = occ.K_S
        
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
