using NCDatasets
using Formatting

using Statistics

Re  = 6371000.0 #  m
c_p  = 1004.0    #  J / kg K
g   = 9.81      #  m / s^2
Lw  = 2.26e6    #  J / kg

t_beg = 240-10*12+1
t_end = 240
t_rng = t_beg:t_end
years = (t_end - t_beg + 1) / 12.0
nc_file_dir = "extract_nc_files"
domain_nc_file_dir = "domain_nc_files"

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file"
            help = "Atm data file. Need variable: "
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
      
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

output_file = joinpath(dirname(parsed["data-file"]), "atm_anomalies.nc")


casenames = [
    "lowres_STD_SOM",
    "lowres_SSM_SOM",
    "lowres_SSM_NK",
    "lowres_SSM_SOM_noQflux",
]
    
mreplace = (x,) -> replace(x, missing=>NaN)

Dataset("$domain_nc_file_dir/domain.lnd.fv4x5_gx3v7.091218.nc", "r") do ds
    global mask
    mask = ds["mask"][:] |> mreplace
    
    mask[mask.!=0] .= 1
end



Dataset("$nc_file_dir/$(casenames[1]).h1.nc", "r") do ds
    global lat, lon, lev, ilev, rec

    lat = ds["lat"][:] |> mreplace
    lon = ds["lon"][:] |> mreplace

    ilev = (ds["ilev"][:] |> mreplace) * 100.0

    rec = ds["time"][:] |> mreplace

    lev = (ilev[2:end] + ilev[1:end-1]) / 2.0
end

# mass
Δp = (ilev[2:end] - ilev[1:end-1])


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

end

