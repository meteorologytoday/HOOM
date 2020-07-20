
using NCDatasets
using Formatting
using ArgParse
using Statistics
using JSON

include("./lib/map_transform.jl")
include("constants.jl")
include("CESMReader.jl")

using .CESMReader
using .MapTransform

function mreplace(x)
    return replace(x, missing=>NaN)
end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin
 
        "--data-file-prefix"
            help = "Data filename prefix including folder and path until the timestamp. File extension `nc` is assumed."
            arg_type = String
            required = true
 
        "--data-file-timestamp-form"
            help = "Data filename timestamp form. Either `YEAR` or `YEAR_MONTH`."
            arg_type = String
            required = true

        "--output-file"
            help = "Output file."
            arg_type = String
            required = true

        "--domain-file"
            help = "Ocn domain file."
            arg_type = String
            required = true
 
        "--beg-year"
            help = "Year of begin."
            arg_type = Int64
            required = true

        "--end-year"
            help = "Year of end."
            arg_type = Int64
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

Dataset(parsed["domain-file"], "r") do ds
    global area = mreplace(ds["area"][:])
    global mask = mreplace(ds["mask"][:])
    global lat  = mreplace(ds["yc"][:])

    area .*= 4π * Re^2 / sum(area)
end

lat_bnd = collect(Float64, -90:1:90)
r = MapTransform.Relation(
    lat = lat,
    area = area,
    mask = mask,
    lat_bnd = lat_bnd,
)

_proxy = area * 0 .+ 1.0
sum_valid_area = MapTransform.∫∂a(r, _proxy)[end]

println("Sum of valid area: ", sum_valid_area, "; ratio: ", sum_valid_area / sum(area))


let
    global TFLUX_DIV_implied, OHT

    if parsed["data-file-timestamp-form"] == "YEAR"
        filename_format = format("{:s}{{:04d}}.nc", joinpath(parsed["data-file-prefix"]))
        form = :YEAR
    elseif parsed["data-file-timestamp-form"] == "YEAR_MONTH"
        filename_format = format("{:s}{{:04d}}-{{:02d}}.nc", joinpath(parsed["data-file-prefix"]))
        form = :YEAR_MONTH
    end
   
    fh = FileHandler(filename_format=filename_format, form=form)

    beg_t = (parsed["beg-year"] - 1) * 12 + 1
    end_t = (parsed["end-year"] - 1) * 12 + 12
 
    TFLUX_DIV, seaice_nudge_energy = getData(fh, ["TFLUX_DIV_implied", "seaice_nudge_energy"], (parsed["beg-year"], parsed["end-year"]), (:, :))

    _data = mean(TFLUX_DIV, dims=(3, ))[:, :, 1] * ρc + mean(seaice_nudge_energy, dims=(3, ))[:, :, 1] 


    TFLUX_DIV_implied = MapTransform.transform(r, _data) 
    OHT = MapTransform.∫∂a(r, _data)

end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "lat_bnd", length(r.lat_bnd))
    defDim(ds, "lat",     length(r.lat_bnd)-1)

    for (varname, vardata, vardim, attrib) in [
        ("TFLUX_DIV_implied", TFLUX_DIV_implied, ("lat",), Dict()),
        ("OHT",               OHT,               ("lat_bnd",), Dict()),
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

