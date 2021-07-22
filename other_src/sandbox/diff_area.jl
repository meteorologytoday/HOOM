include("SMARTSLAB-main/src/share/DisplacedPoleCoordinate.jl")
include("SMARTSLAB-main/src/share/MapInfo.jl")
include("SMARTSLAB-main/src/share/constants.jl")

using .ModelMap
using .DisplacedPoleCoordinate

using NCDatasets
using ArgParse
using JSON

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--domain-file"
            help = "Domain file, dimension names of lon, lat and varname of mask. Separated by comma. Ex: xxx.nc,ni,nj,mask"
            arg_type = String
            required = true
 
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))


mi = ModelMap.MapInfo{Float64}(parsed["domain-file"])
gi = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv, mi.area; angle_unit=:deg)



area_CESM  = mi.area / sum(mi.area) * 4π * Re^2
area_SMART = gi.dσ   / sum(gi.dσ)   * 4π * Re^2
area_diff  = area_SMART - area_CESM
area_diff_ratio = (area_diff ./ area_CESM)

Dataset("area_difference.nc", "c") do ds


    defDim(ds, "Nx", mi.nx)
    defDim(ds, "Ny", mi.ny)

    for (varname, vardata, vardim, attrib) in [
        ("area_CESM",   area_CESM,  ("Nx", "Ny"), Dict()),
        ("area_SMART",  area_SMART,  ("Nx", "Ny"), Dict()),
        ("area_diff",   area_diff,  ("Nx", "Ny"), Dict()),
        ("area_diff_ratio",   area_diff_ratio,  ("Nx", "Ny"), Dict()),
        ("ds1_SMART",   gi.ds1,  ("Nx", "Ny"), Dict()),
        ("ds2_SMART",   gi.ds2,  ("Nx", "Ny"), Dict()),
    ]

        println("Doing var: ", varname)

        var = defVar(ds, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end
   

end


















