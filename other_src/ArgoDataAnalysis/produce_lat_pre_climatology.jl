using NCDatasets

infos = [
    Dict(
        "file" => "RG_ArgoClim_Psal.nc",
        "var"  => "ARGO_SALINITY_MEAN",
    ), Dict(
        "file" => "RG_ArgoClim_Temp.nc",
        "var"  => "ARGO_TEMPERATURE_MEAN",
    )
]

out_lats = collect(Float64, range(-89.5, stop=89.5, step=1.0))
missing_value = 1e20
Dataset(infos[1]["file"], "r") do ds
    global out_press = ds["PRESSURE"][:]
end


mod360 = x -> mod(x, 360.0)
getLatIndex = in_lat -> findmin(abs.(out_lats .- in_lat))[2]

Dataset("output_count_all.nc", "c") do o_ds

    defDim(o_ds, "lat", length(out_lats))
    defDim(o_ds, "pre", length(out_press))
 
    v = defVar(o_ds, "lat", Float64, ("lat",))
    v[:] = out_lats
 
    v = defVar(o_ds, "pre", Float64, ("pre",))
    v[:] = out_press
  
    mean_profile = zeros(Float64, length(out_lats), length(out_press))
    mean_profile_cnt = zeros(Float64, length(out_lats), length(out_press))

    for info in infos
        
        Dataset(info["file"], "r") do i_ds
            index = getLatIndex.(i_ds["LATITUDE"][:])
            profile = nomissing(i_ds[info["var"]][:], NaN)

            mean_profile .= 0.0
            mean_profile_cnt .= 0.0

            for j = 1:i_ds.dim["LATITUDE"]
                lat_j = index[j]

                # Do as much data as possible
                for i = 1:i_ds.dim["LONGITUDE"], k = 1:i_ds.dim["PRESSURE"]
                    if isnan(profile[i, j, k])
                        continue
                    end

                    mean_profile[lat_j, k] += profile[i, j, k]
                    mean_profile_cnt[lat_j, k] += 1.0
                end

                #=
                # Count if all depth data are present
                for i = 1:i_ds.dim["LONGITUDE"]
                    if any(isnan.(profile[i, j, :])) # only need deep sea profile
                        continue
                    end

                    mean_profile[lat_j, :] += profile[i, j, :]
                    mean_profile_cnt[lat_j, :] .+= 1.0
                end
                =#

            end

            for j = 1:length(out_lats), k = 1:length(out_press)
                if mean_profile_cnt[j, k] == 0.0
                    mean_profile[j, k] = NaN
                else
                    mean_profile[j, k] /= mean_profile_cnt[j, k]
                end
            end


        end

        mean_profile[isnan.(mean_profile)] .= missing_value
        
        v = defVar(o_ds, info["var"], Float64, ("lat", "pre"))
        v.attrib["_FillValue"] = missing_value
        v.attrib["missing_value"] = missing_value
        v[:] = mean_profile 
        
    end

    
end



