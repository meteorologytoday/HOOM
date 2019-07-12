
using NCDatasets
using ArgParse
using JSON

src = joinpath("..", "..", "..", "src")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
        "--data-clim-T-file"
            help = "Climatological temperature data file and its varname separated by comma. Supposed to be 3D."
            arg_type = String
            required = true
 
        "--data-clim-S-file"
            help = "Climatological salinity data file and its varname separated by comma. Supposed to be 3D."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file, dimension names of lon, lat and varname of mask. Separated by comma. Ex: xxx.nc,ni,nj,mask"
            arg_type = String
            required = true
 
        "--domain-z-file"
            help = "Z coordinate file, varname of zlon. Separated by comma."
            arg_type = String
            required = true
 
        "--topo-file"
            help = "Ocean topology file (depth) and varname of depth. Separated by comma"
            arg_type = String
            required = true


        "--T-unit"
            help = "Unit of temperature. Valid string: C, K."
            arg_type = String
            required = true
 
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

if ! (parsed["T-unit"] in ["C", "K"])
    throw(ErrorException("Invalid --T-unit: ", parsed["T-unit"]))
end


println("Processing data...")

ps = split(parsed["domain-file"], ",")
println(ps)
domain_file = String(ps[1])
Dataset(ps |> popfirst!, "r") do ds
    global Nx = ds.dim[ps |> popfirst!]
    global Ny = ds.dim[ps |> popfirst!]
    global mask = convert(Array{Float64}, replace(ds[ps |> popfirst!][:], missing=>NaN))
end

ps = split(parsed["domain-z-file"], ",")
println(ps)
Dataset(ps |> popfirst!, "r") do ds
    global zs  = replace(ds[ps |> popfirst!][:], missing=>NaN)
    global Δzs = zs[1:end-1] - zs[2:end]
end

ps = split(parsed["data-clim-T-file"], ",")
println(ps)
Dataset(ps |> popfirst!, "r") do ds
    global Ts_clim = replace(ds[ps |> popfirst!][:, :, :, 1], missing=>NaN)
end

ps = split(parsed["data-clim-S-file"], ",")
println(ps)
Dataset(ps |> popfirst!, "r") do ds
    global Ss_clim = replace(ds[ps |> popfirst!][:, :, :, 1], missing=>NaN)
end

ps = split(parsed["topo-file"], ",")
println(ps)
topo_file = String(ps[1])
Dataset(ps |> popfirst!, "r") do ds
    global topo = zeros(Float64, 1, 1)
    topo = - replace(ds[ps |> popfirst!][:], missing => NaN)
end

if parsed["T-unit"] == "C"
    Ts_clim .+= 273.15
end

