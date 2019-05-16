using NCDatasets
using Printf
using Formatting
using Statistics: mean

ρ    = 1027.0  # kg / m^3
g    = 9.8     # m / s^2

ifilename_Temp = "RG_ArgoClim_Temp.nc"
ifilename_Psal = "RG_ArgoClim_Psal.nc"

month = 7

Dataset(ifilename_Temp,"r") do ds
    
    global lon  = replace(ds["LONGITUDE"][:], missing=>NaN)
    global lat  = replace(ds["LATITUDE"][:], missing=>NaN)
    global pres = replace(ds["PRESSURE"][:], missing=>NaN)
    global time = collect(Float64, 1:length(ds["TIME"]))

    global T_mean =   replace(ds["ARGO_TEMPERATURE_MEAN"][:], missing=>NaN)
    T_mean += mean(replace(ds["ARGO_TEMPERATURE_ANOMALY"][:, :, :, month:12:end], missing=>NaN), dims=(4,))[:, :, :, 1]

end

Dataset(ifilename_Psal,"r") do ds
    
    global S_mean = replace(ds["ARGO_SALINITY_MEAN"][:], missing=>NaN)
    S_mean += mean(replace(ds["ARGO_SALINITY_ANOMALY"][:, :, :, month:12:end], missing=>NaN), dims=(4,))[:,:,:,1]

end

valid_data = isfinite.(T_mean)
replace!(T_mean, NaN=>0.0)
replace!(S_mean, NaN=>0.0)

time /= 12.0
z = pres * 1e4 / (ρ * g)


llon = repeat(reshape(lon, :, 1, 1) , outer=(1, length(lat), length(z)))
llat = repeat(reshape(lat, 1, :, 1), outer=(length(lon), 1, length(z)))

rngs = [ (l-2,l+2)  for l in [-60.0, -30.0, 0.0, 30.0, 60.0] ]
rng_lon = (0, 360)#169, 171)

using PyPlot

fig, axes = plt[:subplots](1,length(rngs),figsize=(12,6), sharey=true)

fig[:suptitle](format("Argo float climatology ({:02d} years) of temperature (solid) and salinity (dashed) in month {:d}. Longitude: [{:d}, {:d}] ", floor(Int64, length(time) / 12), month, rng_lon[1], rng_lon[2]))

for i in 1:length(rngs)

    rng = rngs[i]
    ax  = axes[i]

    blk = (llat .< rng[1]) .| (llat .> rng[2]) .| (llon .> rng_lon[2]) .| (llon .< rng_lon[1])

    tmpT = copy(T_mean)
    tmpT[blk] .= 0.0

    tmpS = copy(S_mean)
    tmpS[blk] .= 0.0

    tmpVD = copy(valid_data)
    tmpVD[blk] .= false
    
    cnts = sum(tmpVD, dims=(1,2))[1, 1, :]

    T_avg = sum(tmpT, dims=(1,2,))[1, 1, :] ./ cnts
    S_avg = sum(tmpS, dims=(1,2,))[1, 1, :] ./ cnts


    mid = (rng[1]+rng[2]) / 2

    
    suf = let
        if mid > 0
            "N"
        elseif mid == 0
            "E"
        else
            "S"
        end
    end
   
    ax[:invert_yaxis]()
    ax2 = ax[:twiny]()

    ax[:plot](T_avg, z, "k-", label="Temperature")
    ax2[:plot](S_avg, z, "k--", label="Salinity")
    ax[:text](25, 850, format("{:d}{}\$\\pm\$2", abs(mid), suf), ha="center")
    ax[:set_xlim]([0, 30])
    ax2[:set_xlim]([32, 36])
    ax[:set_ylim]([1000, -20])
    
    if i == 3
        ax[:set_xlabel]("Temperature [deg C]")
        ax2[:set_xlabel]("Salinity [g/kg]")
    end

end

axes[1][:set_ylabel]("Depth [m]")
plt[:show]()
fig[:savefig](format("Argofloat_climatology_m{:02d}_lon{:d}-{:d}.png", month, rng_lon[1], rng_lon[2]), dpi=200)

