function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds
        occ = makeBlankOceanColumnCollection(
            ds.dim["Nx"], ds.dim["Ny"], nomissing(ds["zs"][:], NaN);
            mask  = nomissing(ds["mask"][:], NaN),
        )


        occ.K = ds.attrib["K"]
        occ.b_ML[:, :]  = nomissing( ds["b_ML"][:], NaN )
        occ.h_ML[:, :]  = nomissing( ds["h_ML"][:], NaN )
        occ.FLDO[:, :]  = nomissing( ds["FLDO"][:], NaN )
        occ.bs[:, :, :] = nomissing( ds["bs"][:], NaN )

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
   
        _write2NCFile(ds, "zs", ("N_zs",), occ.zs, missing_value)
        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)

        _write2NCFile(ds, "bs", ("Nx", "Ny", "N_lays"), occ.bs, missing_value)
        _write2NCFile(ds, "b_ML", ("Nx", "Ny",), occ.b_ML, missing_value)
        _write2NCFile(ds, "h_ML", ("Nx", "Ny",), occ.h_ML, missing_value)
        _write2NCFile(ds, "FLDO", ("Nx", "Ny",), convert(Array{Float64}, occ.FLDO), missing_value)
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
        ds.attrib["K"] = occ.K
        
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
