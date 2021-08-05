using NCDatasets
using Formatting


mreplace = (x,) -> replace(x, missing => NaN)

sic_file = "data/monthly_clim.nc" 

Dataset(sic_file, "r") do ds

    global area = ds["area"][:] |> mreplace
    global mask = ds["mask"][:] |> mreplace
    global lat  = ds["lat"][:] |> mreplace
    global sic  = ds["IFRAC_clim"][:] |> mreplace

end

area = 6371e3^2.0 * 4π * area / sum(area)
sic ./= 100.0

NH_idx = (lat .>= 0.0) .& (mask .== 1)
SH_idx = (lat .< 0.0) .& (mask .== 1)
ALL_idx = mask .== 1

total_ice_area = zeros(Float64, 12)
total_ice_area_NH = zeros(Float64, 12)
total_ice_area_SH = zeros(Float64, 12)

for m=1:12

    ice_area = area .* sic[:, :, m]

    total_ice_area[m]    = sum(ice_area[ALL_idx])
    total_ice_area_NH[m] = sum(ice_area[NH_idx])
    total_ice_area_SH[m] = sum(ice_area[SH_idx])

end

Dataset("total_ice_area.nc", "c") do ds

    defDim(ds, "time", Inf)

    for (varname, vardata, vardim, attrib) in [

        ("total_ice_area", total_ice_area, ("time",), Dict(
            "units"     => "m^2",
            "long_name" => "total seaice area",
        )),

        ("total_ice_area_NH", total_ice_area_NH, ("time",), Dict(
            "units"     => "m^2",
            "long_name" => "total seaice area in NH",
        )),


        ("total_ice_area_SH", total_ice_area_SH, ("time",), Dict(
            "units"     => "m^2",
            "long_name" => "total seaice area in SH",
        )),

    ]

        var = defVar(ds, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20
        
        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        println(var.attrib)

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end 
