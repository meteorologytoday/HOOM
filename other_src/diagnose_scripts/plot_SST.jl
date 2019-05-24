include("nanop.jl")

using NCDatasets
using Formatting
using Statistics


Re  = 6371000.0 #  m
c_p  = 1004.0    #  J / kg K
g   = 9.81      #  m / s^2
Lw  = 2.26e6    #  J / kg

casenames = [
    "lowres_SSM_NK",
    "lowres_SSM_SOM",
    "lowres_SSM_SOM_noQflux",
]
    
linestyle = ["k-", "r-", "g-", "b-"]

Dataset("/glade/p/cesmdata/inputdata/share/domains/domain.lnd.fv4x5_gx3v7.091218.nc", "r") do ds
    global lat, lon, lev, ilev, rec

    lat = replace(ds["yc"][1,:], missing=>NaN)
    lon = replace(ds["xc"][:,1], missing=>NaN)


end


using PyPlot

fig, ax = plt[:subplots](1, 1, figsize=(12,8))
ax[:plot]([-90, 90], [0, 0], linewidth=2, color="#cccccc")

for i = 1:length(casenames)

    casename = casenames[i]

    Dataset("$casename/archive/ocn/$casename.ocn.h.ma.fv4x5.0001-0020.nc", "r") do ds

        t_rng = 240-5*12+1:240

        if match(r"SSM_NK", casename) != nothing
            T_ML = ds["T"][:, :, 1, t_rng]
        elseif match(r"SSM_SOM", casename) != nothing
            T_ML = ds["T_ML"][:, :, t_rng]
        else
            throw(ErrorException("ERROR"))
        end
    
        T_ML = replace(T_ML, missing=>NaN)
        T_ML = nanmean(T_ML, dims=(1,3,))[1, :, 1]

        ax[:plot](lat, T_ML, linestyle[i], linewidth=2,  label=casename)
        
    end

end

ax[:legend]()
ax[:set_xlabel]("Latitude [deg]")
ax[:set_ylabel]("Sea Surface Temperature [degC]")
ax[:set_xlim](-90, 90)

fig[:savefig](format("compare_models_SST.png", casename), dpi=200)
