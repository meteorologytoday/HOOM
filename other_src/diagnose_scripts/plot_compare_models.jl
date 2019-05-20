using NCDatasets
using Formatting

using Statistics

Re  = 6371000.0 #  m
c_p  = 1004.0    #  J / kg K
g   = 9.81      #  m / s^2
Lw  = 2.26e6    #  J / kg

casenames = [
    "lowres_STD_SOM",
    "lowres_SSM_SOM",
    "lowres_SSM_NK",
    "lowres_SSM_SOM_noQflux",
]
    
linestyle = ["k-", "r-", "g-", "b-"]

Dataset(casenames[1] * ".h1.nc", "r") do ds
    global lat, lon, lev, ilev, rec

    lat = ds["lat"][:] 
    lon = ds["lon"][:] 
#    lev = ds["lev"][:]
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
        v[j, k, l] *= Δp[k] * cos(lat[j] |>  deg2rad)
    end

    v .*= 2 * π * Re / g

    # Get the vertical integration and time mean
    v = sum(v, dims=(2,))[:, 1, :]
    v = mean(v, dims=(2,))[:, 1]


    return v
end

function integrate(x, dydx)
    y = copy(dydx) * 0.0

    for i = 2:length(x)
        y[i] = y[i-1] + (dydx[i-1] + dydx[i]) * (x[i] - x[i-1]) / 2.0
    end

    return y
end


print("Import Pyplot... ")
using PyPlot
println("Done. ")

#=

### Meridional Static Energy Transport

fig, ax = plt[:subplots](1, 1, figsize=(12,8))

ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

for i = 1:length(casenames)

    casename = casenames[i]

    Dataset(casename * ".nc", "r") do ds
        global MHT, MGT, MLT
        MHT = c_p * ds["VT"][:]  |> calMeanTransport
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



### Precipitation
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

=#
### Energy convergence


fig_met, ax_met = plt[:subplots](1, 1, figsize=(12,8))
ax_met[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

for i = 1:length(casenames)

    fig, ax = plt[:subplots](1, 1, figsize=(12,8))
    ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

    casename = casenames[i]

    Dataset(casename * ".h0.nc", "r") do ds
        t_rng = (:, :,  240-5*12+1:240)
        FFLX_TOA = mean( - mean(ds["FSNT"][t_rng...], dims=(3,)) + mean(ds["FLNT"][t_rng...], dims=(3,)), dims=(1,))[1, :, 1]
        FFLX_SFC = mean( - mean(ds["FSNS"][t_rng...], dims=(3,)) + mean(ds["FLNS"][t_rng...], dims=(3,)), dims=(1,))[1, :, 1]
        HFLX_SFC = mean( mean(ds["SHFLX"][t_rng...], dims=(3,)) + mean(ds["LHFLX"][t_rng...], dims=(3,)), dims=(1,))[1, :, 1]
       

        FFLX_TOA .*= Re * cos.(deg2rad.(lat)) * 2π 
        FFLX_SFC .*= Re * cos.(deg2rad.(lat)) * 2π
        HFLX_SFC .*= Re * cos.(deg2rad.(lat)) * 2π
 
        EFLX_CONV = - ( FFLX_TOA - FFLX_SFC - HFLX_SFC )


        
        ax[:plot](lat, FFLX_TOA, "r--", linewidth=2,  label="FTOA")
        ax[:plot](lat, FFLX_SFC, "b--", linewidth=2,  label="FSFC")
        ax[:plot](lat, HFLX_SFC, "g--", linewidth=2,  label="HFLX")
        ax[:plot](lat, EFLX_CONV, "k-", linewidth=2,  label="CONV")
        
        MET = integrate(Re * deg2rad.(lat), EFLX_CONV)
        ax_met[:plot](lat, MET / 1e15, linestyle[i], linewidth=2, label=casename)

    end

    ax[:legend]()
    ax[:set_xlabel]("Latitude [deg]")
    ax[:set_ylabel]("[PW / m]")

    ax[:set_title](format("10 years in case {}", casename))
    ax[:set_xlim](-90, 90)

    fig[:savefig](format("compare_models_EBM_{}.png", casename), dpi=200)



end

ax_met[:legend]()
ax_met[:set_xlabel]("Latitude [deg]")
ax_met[:set_ylabel]("[PW]")

ax_met[:set_title](format("16-20 years annual implied northward heat transport"))
ax_met[:set_xlim](-90, 90)
ax_met[:set_ylim](-6, 6)

fig_met[:savefig](format("compare_models_MET.png"), dpi=200)


#ax[:set_ylim](-15, 20)



### Temperature trend

fig, ax = plt[:subplots](1, 1, figsize=(12,8))


dθ = deg2rad(lat[2] - lat[1])
for i = 1:length(casenames)

    casename = casenames[i]

    Dataset(casename * ".h1.nc", "r") do ds

        TREFHT = sum(mean(ds["TREFHT"][:], dims=(1,)), dims=(1,) )[1, :, :]
        Tavg = zeros(Float64, length(rec))
        for l = 1:length(rec)
            TREFHT[:, l] .*= cos.(deg2rad.(lat))
            Tavg[l] = sum(TREFHT[:, l]) * dθ / 2.0
        end

        Tavg .-= 273.15

        ax[:plot](collect(Float64, 0:length(rec)-1) / 365.0, Tavg, linestyle[i], linewidth=2,  label=casename)
    end

end

ax[:legend]()
ax[:set_xlabel]("Time [month]")
ax[:set_ylabel]("Temperature")

ax[:set_title]("Global Mean surface temperature of 20 years")
ax[:set_xlim](0, 20)
ax[:set_ylim](0, 20)

fig[:savefig]("compare_models_ECNVG.png", dpi=200)


