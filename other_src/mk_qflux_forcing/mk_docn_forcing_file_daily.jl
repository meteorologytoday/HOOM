using NCDatasets
using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "File containing climatology of SST and Mixed-layer depth."
            arg_type = String
            required = true
        
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

        "--output-file"
            help = "Output file name."
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))


nomissing = (x) -> replace(x, missing=>NaN)

Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)
end

Dataset(parsed["input-file"], "r") do ds
    global Tclim = ds["Tclim"][:] |> nomissing
    global hblt  = ds["hblt"][:] |> nomissing
end 

dom = convert(Array{Float64}, [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
time = collect(Float64, 1:365) .- 0.5

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", length(time))

    empty = zeros(Float64, ni, nj, length(time))

    for (varname, vardata, vardim, attrib) in [

        ("time", time, ("time",), Dict(
            "calendar"  => "noleap",
            "long_name" => "observation time",
            "units"     => "days since 0001-01-01 00:00:00",
        ) ),

        ("area", area, ("ni", "nj"), Dict(
            "units"     => "radians squared",
            "long_name" => "area of grid cell",
        )),

        ("mask", mask, ("ni", "nj"), Dict(
            "units"     => "unitless",
            "long_name" => "domain mask",
        )),


        ("xc", xc, ("ni", "nj"), Dict(
            "units"     => "degrees east",
            "long_name" => "longitude of grid cell center",
        )),

        ("yc", yc, ("ni", "nj"), Dict(
            "units"     => "degrees north",
            "long_name" => "latitude of grid cell center",
        )),

        ("qdp", empty, ("ni", "nj", "time"), Dict(  
            "units"     => "W/m^2",
            "long_name" => "ocean heat flux convergence",
        )),

        ("hblt", hblt, ("ni", "nj", "time"), Dict(
            "units"     => "m",
            "long_name" => "Mixed-layer depth",
        )),

        ("Tclim", Tclim, ("ni", "nj", "time"), Dict(
            "units"     => "m",
            "long_name" => "Climatology of sea-surface-temperature",
        )),

        ("S",    empty , ("ni", "nj", "time"), Dict()),
        ("T",    empty , ("ni", "nj", "time"), Dict()),
        ("U",    empty , ("ni", "nj", "time"), Dict()),
        ("V",    empty , ("ni", "nj", "time"), Dict()),
        ("dhdx", empty , ("ni", "nj", "time"), Dict()),
        ("dhdy", empty , ("ni", "nj", "time"), Dict()),
    ]

        var = defVar(ds, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20
        
        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        println(var.attrib)

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end 
