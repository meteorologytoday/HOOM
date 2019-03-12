using NCDatasets

function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds
        occ = OceanColumnCollection(
            Nx    = ds.dim["Nx"],
            Ny    = ds.dim["Ny"],
            N     = ds.dim["N_lays"],
            zs    = nomissing(ds["zs"][:], NaN),
            bs    = zeros(Float64, ds.dim["N_lays"]),
            K     = ds.attrib["K"],
            b_ML  = 0.0,
            h_ML  = 0.0,
            FLDO  = 0,
            mask  = nomissing(ds["mask"][:], NaN),
        )

        b_ML = nomissing( ds["b_ML"][:], NaN )
        h_ML = nomissing( ds["h_ML"][:], NaN )
        FLDO = nomissing( ds["FLDO"][:], NaN )

        for i = 1:occ.Nx, j = 1:occ.Ny

            if occ.mask_idx[i, j]
                continue
            end    

            occ.ocs[i, j].bs .= ds["bs"][i, j, :]
            occ.ocs[i, j].b_ML = b_ML[i, j]
            occ.ocs[i, j].h_ML = h_ML[i, j]
            occ.ocs[i, j].FLDO = FLDO[i, j]
        end

    end
    return occ 
end


function takeSnapshot(occ::OceanColumnCollection, filename::AbstractString; missing_value::Float64=1e20)

    bs   = zeros(Float64, occ.Nx, occ.Ny, occ.N)
    b_ML = zeros(Float64, occ.Nx, occ.Ny)
    h_ML = zeros(Float64, occ.Nx, occ.Ny)
    FLDO = zeros(Float64, occ.Nx, occ.Ny)

    for i=1:occ.Nx, j=1:occ.Ny

        if occ.mask_idx[i, j]
            bs[i, j, :] .=  missing_value
            b_ML[i, j] = missing_value
            h_ML[i, j] = missing_value
            FLDO[i, j] = Float64(missing_value)
        else
            bs[i, j, :] .=  occ.ocs[i, j].bs
            b_ML[i, j] = occ.ocs[i, j].b_ML
            h_ML[i, j] = occ.ocs[i, j].h_ML
            FLDO[i, j] = Float64(occ.ocs[i, j].FLDO)
        end
    end



    Dataset(filename, "c") do ds

        defDim(ds, "N_ocs", occ.N_ocs)
        defDim(ds, "Nx", occ.Nx)
        defDim(ds, "Ny", occ.Ny)
        defDim(ds, "N_lays", occ.N)
        defDim(ds, "N_zs",   length(occ.zs))
       
        ds.attrib["_FillValue"] = missing_value
        ds.attrib["K"] = occ.K
        
        _write2NCFile(ds, "zs", ("N_zs",), occ.zs, missing_value)
        _write2NCFile(ds, "mask", ("Nx", "Ny",), occ.mask, missing_value)

        _write2NCFile(ds, "bs", ("Nx", "Ny", "N_lays"), bs, missing_value)
        _write2NCFile(ds, "b_ML", ("Nx", "Ny",), b_ML, missing_value)
        _write2NCFile(ds, "h_ML", ("Nx", "Ny",), h_ML, missing_value)
        _write2NCFile(ds, "FLDO", ("Nx", "Ny",), FLDO, missing_value)
    end

end

function _write2NCFile(
    ds            :: Dataset,
    varname       :: AbstractString,
    dim           :: Tuple,
    var_data      :: AbstractArray{T},
    missing_value :: T) where T 

    println("Write : ", varname)

    ds_var = defVar(ds, varname, eltype(var_data), dim)
    ds_var.attrib["_FillValue"] = missing_value
    ds_var[:] = var_data
end
