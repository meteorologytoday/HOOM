using NCDatasets
using Printf
using Formatting

include("config.jl")


Buoy = (α * readModelVar(ifilename_Temp, "ARGO_TEMPERATURE_ANOMALY") - β * readModelVar(ifilename_Psal, "ARGO_SALINITY_ANOMALY")) * g


#, range_tuple=(lon_i, lat_i, :, :))

println("Start working...")
for k = 1:length(pres), j=1:length(lat), i=1:length(lon)

    if i==1 && j==1
        println("Doing pressure level: ", pres[k])
    end

    if any(isnan.(Buoy[i,j,k,:]))
        Buoy[i,j,k,:] .= NaN
        continue
    end

#    Buoy[i, j, k, :] = rmSeasonalCycle(x=Buoy[i, j, k, :], period=12)
end
println("done.")


nan2missing!(Buoy)

ofilename = "RG_ArgoClim_Buoy_rmSeasonalCycle.nc"
ofilename = "RG_ArgoClim_Buoy.nc"

outtype = Float32

ds = Dataset(ofilename,"c")
defDim(ds,"time", length(time))
defDim(ds,"lat", length(lat))
defDim(ds,"lon", length(lon))
defDim(ds,"pres", length(pres))

defVar(ds, "time", outtype, ("time",))[:] = time
defVar(ds, "lat",  outtype, ("lat",))[:]  = lat
defVar(ds, "lon",  outtype, ("lon",))[:]  = lon
defVar(ds, "pres", outtype, ("pres",))[:] = pres

for o in (
    [
        "BUOYA", Buoy, ("lon", "lat", "pres", "time"), Dict(
        #"long_name"=>"Argo buoyancy anomaly (Seasonal cycle removed)",
        "long_name"=>"Argo buoyancy anomaly",
        "units"=>"m / s^2",
        )
    ],
)
    varname, vardata, vardims, varatts = o
    println("Writing ", varname, " with size: ", size(vardata), " ; dim: ", vardims)

    ncvar = defVar(ds, varname, eltype(vardata), vardims)
    ncvar.attrib["_FillValue"] = convert(Float64, missing_value)
    for key in keys(varatts)
        ncvar.attrib[key] = varatts[key]
    end

    ncvar[:] = vardata
    println("done.")
end

close(ds)
println("Output file: ", ofilename)

