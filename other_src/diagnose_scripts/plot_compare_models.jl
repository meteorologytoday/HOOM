using NCDatasets
using Formatting

using Statistics

deg2rad = π/180.0
Re  = 6371000.0 #  m
cp  = 1004.0    #  J / kg K
g   = 9.81      #  m / s^2
Lw  = 2.26e6    #  J / kg

casenames = [
    "lowres_STD_SOM",
    "lowres_SSM_SOM",
    "lowres_SSM_NK",
    "lowres_SSM_SOM_noQflux",
]
    
linestyle = ["k-", "r-", "g-", "b-"]

Dataset(casenames[1] * ".nc", "r") do ds
    global lat, lon, lev, ilev, rec

    lat = ds["lat"][:] 
    lon = ds["lon"][:] 
    lev = ds["lev"][:]
    ilev = ds["ilev"][:]
    rec = ds["time"][:]

end

# mass
Δp = (ilev[2:end] - ilev[1:end-1]) * 100.0


function calMeanTransport(v)

    # Get the zonal mean
    v = mean(v, dims=(1,))[1, :, :, :]

    # Weight with mass
    for j=1:length(lat), k=1:length(lev), l=1:length(rec)
        v[j, k, l] *= Δp[k] * cos(lat[j] * deg2rad)
    end

    v .*= 2 * π * Re / g

    # Get the vertical integration and time mean
    v = sum(v, dims=(2,))[:, 1, :]
    v = mean(v, dims=(2,))[:, 1]


    return v
end

print("Import Pyplot... ")
using PyPlot
println("Done. ")



# Meridional Static Energy Transport

fig, ax = plt[:subplots](1, 1, figsize=(12,8))

ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

for i = 1:length(casenames)

    casename = casenames[i]

    Dataset(casename * ".nc", "r") do ds
        global MHT, MGT, MLT
        MHT = cp * ds["VT"][:]  |> calMeanTransport
        MGT = ds["VZ"][:]       |> calMeanTransport
        MLT = Lw * ds["VQ"][:]  |> calMeanTransport

    end

    MMSET = (MHT + MGT + MLT) / 1e15


    ax[:plot](lat, MMSET, linestyle[i], linewidth=2,  label=casename)

end

ax[:legend]()
ax[:set_xlabel]("Latitude [deg]")
ax[:set_ylabel]("[PW]")

ax[:set_title]("10 years")
ax[:set_xlim](-90, 90)
ax[:set_ylim](-15, 20)

fig[:savefig]("compare_models_MSEFLX.png", dpi=200)



# Precipitation
fig, ax = plt[:subplots](1, 1, figsize=(12,8))

ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

for i = 1:length(casenames)

    casename = casenames[i]

    Dataset(casename * ".nc", "r") do ds
        global tot_prec = replace(ds["PRECC"][:] + ds["PRECL"][:], missing=>NaN)

        tot_prec = mean(tot_prec, dims=(1,3))[1, :, 1] # average over time and longitude
    end

    ax[:plot](lat, tot_prec * 365.0 * 86400.0 * 1000.0, linestyle[i], linewidth=2,  label=casename)
end

ax[:legend]()
ax[:set_xlabel]("Latitude [deg]")
ax[:set_ylabel]("Precipitation [mm / yr]")

ax[:set_title]("10 years")
ax[:set_xlim](-90, 90)
#ax[:set_ylim](, 20)

fig[:savefig]("compare_models_PRECIP.png", dpi=200)



