using NCDatasets

infos = [
    Dict(
        "file" => "RG_ArgoClim_PSal.nc",
        "var"  => "ARGO_SALINITY_MEAN",
    ), Dict(
        "file" => "RG_ArgoClim_Temp.nc",
        "var"  => "ARGO_TEMPERATURE_MEAN",
    )
]

out_lats = collect(Float64, range(-89.5, stop=89.5, step=1.0))

Dataset(infos[1]["file"], "r") do ds
    global out_press = ds["PRESSURE"][:]
end


mod360 = x -> mod(x, 360.0)
getLatIndex = in_lat -> findmin(abs.(out_lats .- in_lat))[2]

Dataset("output.nc", "c") do o_ds

    defDim(o_ds, "lat", length(out_lats))
    defDim(o_ds, "pressure", length(out_press))
  
    mean_profile = zeros(Float64, length(out_lats), length(out_press))
    mean_profile_cnt = zeros(Float64, length(out_lats))

    for info in infos
        
        NCDataset(info["file"], "r") do ds
            index = getLatIndex.(ds[info["LATITUDE"][:]])
            profile = nomissing(ds[info["var"]][:], NaN)

            mean_profile .= 0.0
            mean_profile_cnt .= 0.0

            for j = 1:ds.dim["LATITUDE"]
                lat_j = index[j]
                for i = 1:ds.dim["LONGITUDE"]
                    if any(isnan.(profile[i, j, :])) # only need deep sea profile
                        continue
                    end

                    mean_profile[lat_j, :] += profile(i, j, :)
                    mean_profile_cnt[lat_j] += 1.0
                end

                if mean_profile_cnt[lat_j] == 0.0
                    mean_profile[lat_j, :] .= NaN
                else
                    mean_profile[lat_j, :] /= mean_profile_cnt[lat_j]
                end
            end 

        end

    end
end



