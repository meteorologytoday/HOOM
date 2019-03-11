using NCDatasets

function loadSnapshot(filename::AbstractString)
    local occ

    Dataset(filename, "r") do ds

        occ = OceanColumnCollection(
            N_ocs = ds.attrib["N_cols"],
            N     = ds.attrib["N_lays"],
            zs    = ds["zs"][:],
            bs    = zeros(Float64, ds.attrib["N_lays"]),
            K     = ds.attrib["K"],
            b_ML  = 0.0,
            h_ML  = 0.0,
            FLDO  = 0,
            mask  = ds["mask"][:],
        )

    end

    b_ML = ds["b_ML"][:]
    h_ML = ds["h_ML"][:]
    FLDO = ds["FLDO"][:]

    for i = 1:occ.N_ocs

        if occ.mask_idx[i]
            continue
        end    

        occ.ocs[i].bs .= ds["bs"][i, :]
        occ.ocs[i].b_ML = b_ML[i]
        occ.ocs[i].h_ML = h_ML[i]
        occ.ocs[i].FLDO = FLDO[i]
    end

    return occ 
end


function takeSnapshot(occ::OceanColumnCollection, filename::AbstractString; missing_value::Float64=1e20)

    bs = zeros(Float64, occ.N_ocs, occ.N)
    b_ML = zeros(Float64, occ.N_ocs)
    h_ML = zeros(Float64, occ.N_ocs)
    FLDO = zeros(Float64, occ.N_ocs)

    for i=1:occ.N_ocs

        if occ.mask_idx[i]
            bs[i, :] .=  missing_value
            b_ML[i] = missing_value
            h_ML[i] = missing_value
            FLDO[i] = Float64(missing_value)
        else
            bs[i, :] .=  occ.ocs[i].bs
            b_ML[i] = occ.ocs[i].b_ML
            h_ML[i] = occ.ocs[i].h_ML
            FLDO[i] = Float64(occ.ocs[i].FLDO)
        end
    end



    Dataset(filename, "c") do ds

        defDim(ds, "N_ocs", occ.N_ocs)
        defDim(ds, "N_lays", occ.N)
        defDim(ds, "N_zs",   length(occ.zs))
       
        ds.attrib["_FillValue"] = missing_value
        ds.attrib["K"] = occ.K
        
        write2NCFile(ds, "zs", ("N_zs",), occ.zs, missing_value)
        write2NCFile(ds, "mask", ("N_ocs",), occ.mask, missing_value)

        write2NCFile(ds, "bs", ("N_ocs", "N_lays"), bs, missing_value)
        write2NCFile(ds, "b_ML", ("N_ocs",), b_ML, missing_value)
        write2NCFile(ds, "h_ML", ("N_ocs",), h_ML, missing_value)
        write2NCFile(ds, "FLDO", ("N_ocs",), FLDO, missing_value)
    end

end

function write2NCFile(
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
