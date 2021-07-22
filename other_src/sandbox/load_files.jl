
using NCDatasets
using ArgParse
using JSON

src = normpath(joinpath(@__DIR__, "..", "..", "..", "src"))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file, dimension names of lon, lat and varname of mask. Separated by comma. Ex: xxx.nc,ni,nj,mask"
            arg_type = String
            required = true
 
        "--zdomain-file"
            help = "Z coordinate file, varname of zlon. Separated by comma."
            arg_type = String
            required = true
 
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

println("Processing data...")

Dataset(parsed["domain-file"], "r") do ds
    global Nx = ds.dim["ni"]
    global Ny = ds.dim["nj"]
    global mask = convert(Array{Float64}, replace(ds["mask"][:], missing=>NaN))
end

Dataset(parsed["zdomain-file"], "r") do ds
    global zs  = replace(ds["zs"][:], missing=>NaN)
    global Δzs = zs[1:end-1] - zs[2:end]
end

