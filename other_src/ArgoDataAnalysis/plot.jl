include("config.jl")
using PyPlot

filename = "RG_ArgoClim_Temp_rmSeasonalCycle.nc"

lat_i = 96
lon_i = 171


time /= 12.0
z = - pres * 1e4 / (ρ * g)

TEMPA = readModelVar(filename, "TEMPA", ((lon_i-5):(lon_i+5), (lat_i-2:lat_i+2), :, :))

println(size(TEMPA))
TEMPA = mean(TEMPA, dims=(1,2))[1,1,:,:]

fig, ax = plt[:subplots](1,1,figsize=(16,4))

fig[:suptitle](format("Lat={:.2f}, Lon={:.2f}", lat[lat_i], lon[lon_i]))

clevs = (range(-1, stop=1, length=51) |> collect ) * 0.1
cmap = plt[:get_cmap]("jet")
cbmapping = ax[:contourf](time, z, TEMPA, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)

cb = plt[:colorbar](cbmapping, ax=ax)

cb[:set_label]("Temperature anomaly [\$\\mathrm{K} \$]")


ax[:set_xticks](1:15)
plt[:show]()


