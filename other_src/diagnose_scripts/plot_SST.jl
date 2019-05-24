using NCDatasets
using Formatting

using Statistics

if length(ARGS) >= 1
    fn = ARGS[1]
end


function nanmean(arr, dim)
    

end



const Re = 6371000.0 
const cp = 1004.0

ds = Dataset(fn, "r")

#lat = ds["lat"][:] 
#lon = ds["lon"][:] 
#lev = ds["lev"][:]
#ilev = ds["ilev"][:]

#deg2rad = π/180.0
lat = collect(range(-90, stop=90, length=ds.dim["Ny"]))

T_ML = replace(ds["T_ML"][:], missing=>NaN) .- 273.15
T_ML = mean(T_ML, dims=(3,))[:, :, 1]

cnts = sum(isfinite.(T_ML), dims=(1,))[1, :]
T_ML[isnan.(T_ML)] .= 0.0
T_ML = sum(T_ML, dims=(1,3))[1, :] ./ cnts

close(ds)


using PyPlot

fig, ax = plt[:subplots](1, 1, figsize=(6,4))

ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc", dashes=(5,3))
ax[:plot](lat, T_ML, linewidth=2, color="black")
ax[:set_xlabel]("Latitude [deg]")
ax[:set_ylabel]("Sea Surface Temperature [deg C]")

ax[:set_title](fn)
ax[:set_xlim](-90, 90)
#ax[:set_ylim](-100, 2500)

plt[:show]()




