src = normpath(joinpath(@__DIR__, "SMARTSLAB-main", "src"))
include(joinpath(src, "models", "HOOM", "HOOM.jl"))

using NCDatasets
using ArgParse
using JSON

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
 
        "--topo-file"
            help = "Topography file."
            arg_type = String
            default  = ""
 
        "--zdomain-file"
            help = "Z coordinate file, varname of zlon. Separated by comma."
            arg_type = String
            required = true

        "--flip-mask"
            action = :store_true
 
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

mask .= 1.0
if parsed["flip-mask"]
    mask = 1.0 .- mask
end
mask[:, 1]   .= 0.0
mask[:, end] .= 0.0

Dataset(parsed["zdomain-file"], "r") do ds
    global zs  = replace(ds["zs"][:], missing=>NaN)
end


zs = collect(Float64, 0:-50:-1000)

if parsed["topo-file"] == ""
    global topo = nothing
else
    Dataset(parsed["topo-file"], "r") do ds
        global topo  = - replace(ds["depth"][:], missing=>NaN)
    end
end




#mask = 1.0 .- mask
#zs = collect(Float64, 0:-10:-400)

#mask .= 1.0

#mask[:,   1] .= 0.0
#mask[:, end] .= 0.0
#=
for i=1:Nx, j=1:Ny
    if (10 < j < 25) && ( i+j == 50 || i+j == 51 )
        mask[i, j] = 0.0
    end
end
=#
    
Ts_clim = collect(Float64, 0:length(zs)-2)/ (length(zs)-1) * (-0.0) .+ 10.0
Ss_clim = collect(Float64, 0:length(zs)-2)/ (length(zs)-1) * 5 .+ 35

ocn = HOOM.Ocean(
    gridinfo_file = parsed["domain-file"],
    Nx       = Nx,
    Ny       = Ny,
    zs_bone  = zs,
    Ts       = Ts_clim .+ 5.0,
    Ss       = Ss_clim .+ 3.0,
    T_ML     = 10.0,
    S_ML     = 35.0,
    h_ML     = 10.0, 
    h_ML_min = 10.0,
    h_ML_max = 1e5,             # make it unrestricted
    topo     = topo,
    Ts_clim_relax_time = 86400.0 * 10,
    Ts_clim            = nothing, #Ts_clim,
    Ss_clim_relax_time = 86400.0 * 10,
    Ss_clim            = nothing, #Ss_clim,
    arrange  = :xyz,
    do_convective_adjustment = true,
)

HOOM.takeSnapshot(ocn, parsed["output-file"])


