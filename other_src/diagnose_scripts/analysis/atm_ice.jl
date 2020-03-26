using Statistics
using NCDatasets

using ArgParse
using JSON
using Formatting

include("LinearRegression.jl")
include("CESMReader.jl")
include("constants.jl")
using .CESMReader


correlation = (x1, x2) -> x1' * x2 / (sum(x1.^2)*sum(x2.^2)).^0.5

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
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--output-file"
            help = "Output atm temperature file."
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

output_file = parsed["output-file"]

let

    if parsed["data-file-timestamp-form"] == "YEAR"
        filename_format = format("{:s}{{:04d}}.nc", joinpath(parsed["data-file-prefix"]))
        form = :YEAR
    elseif parsed["data-file-timestamp-form"] == "YEAR_MONTH"
        filename_format = format("{:s}{{:04d}}-{{:02d}}.nc", joinpath(parsed["data-file-prefix"]))
        form = :YEAR_MONTH
    end
   
    fh = FileHandler(filename_format=filename_format, form=form)

    #global hi   = getData(fh, "hi",   (parsed["beg-year"], parsed["end-year"]), (:, :))
    global aice = getData(fh, "ICEFRAC", (parsed["beg-year"], parsed["end-year"]), (:, :))

    #aice ./= 100.0

    global (Nx, Ny, Nt) = size(aice)
    global nyears = Int64(Nt / 12)
    
    if mod(Nt, 12) != 0
        ErrorException("Time record is not multiple of 12") |> throw
    end

end

Dataset(parsed["domain-file"], "r") do ds
    global mask = replace(ds["mask"], missing=>NaN)
    global area = replace(ds["area"], missing=>NaN)
    global lat  = replace(ds["yc"], missing=>NaN)

    global mask_idx = (mask .== 0.0)
    global total_area = sum(area)
end

total_ice_volume = zeros(Float64, Nt)
total_ice_area   = zeros(Float64, Nt)
total_ice_extent_area = zeros(Float64, Nt)


x = collect(Float64, 1:Nt)

# total ice volume and area
for i = 1:Nt
    a = area .* view(aice, :, :, i)
    #v = a .* view(hi, :, :, i)  
    
    total_ice_area[i]   = sum( a[mask_idx]  )
    #total_ice_volume[i] = sum( tmp2[mask_idx] )
end

total_ice_area .*= 4π * Re^2.0 / total_area



# total ice extent area
for i = 1:Nt
    tmp = view(aice, :, :, i)

    total_ice_extent_area = sum(area[mask_idx .& (tmp .> .15)])
end

total_ice_extent_area .*= 4π * Re^2.0 / total_area



Dataset(output_file, "c") do ds

    defDim(ds, "months", 12)
    defDim(ds, "time", Inf)
    defDim(ds, "Nx", Nx)
    defDim(ds, "Ny", Ny)

    for (varname, vardata, vardim, attrib) in [
        ("total_ice_volume",  total_ice_volume, ("time",), Dict()),
        ("total_ice_area",    total_ice_area  , ("time",), Dict()),
        ("total_ice_extent_area",    total_ice_area  , ("time",), Dict()),
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
