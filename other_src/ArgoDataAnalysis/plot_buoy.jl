include("config.jl")
using PyPlot

filename = "RG_ArgoClim_Buoy_rmSeasonalCycle.nc"
filename = "RG_ArgoClim_Buoy.nc"

lat_i = 90
lon_i = 171

ρ    = 1027.0  # kg / m^3
g    = 9.8     # m / s^2
c_p  = 3985.0  # J / kg / K


time /= 12.0
z = - pres * 1e4 / (ρ * g)

data = readModelVar(filename, "BUOYA", ((lon_i-2):(lon_i+2), (lat_i-2:lat_i+2), :, :))




println(size(data))
data = mean(data, dims=(1,2))[1,1,:,:]
mean_buoy = mean(data, dims=1)[1, :]
for i = 1:size(data)[2]
    data[:, i] .-= mean(data[:, i])
end



fig, ax = plt[:subplots](2,1,figsize=(12, 6), sharex=true)

fig[:suptitle](format("Lat={:.2f}, Lon={:.2f}, [5 deg avg]", lat[lat_i], lon[lon_i]))

ax[1][:plot](time, mean_buoy)

clevs = (range(-1, stop=1, length=51) |> collect ) * 1.0 
cmap = plt[:get_cmap]("jet")
cbmapping = ax[2][:contourf](time, z, data * 1e3, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)

cb = plt[:colorbar](cbmapping, ax=ax)

cb[:set_label]("Buoyancy anomaly [\$\\times \\, 10^{-3} \\, \\mathrm{m} \\, \\mathrm{s}^{-2} \$]")


ax[2][:set_xticks](1:15)
plt[:show]()


