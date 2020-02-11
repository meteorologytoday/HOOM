using NCDatasets
using ArgParse
using JSON
using Statistics

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Qflux file"
            arg_type = String
            required = true

        "--Qflux-blank"
            help = "If turned on then --Qflux-varname is ignored and output qflux = 0.0 field."
            action = :store_true

        "--Qflux-T-varname"
            help = "T Qflux variable name in file"
            arg_type = String
            default  = ""

        "--Qflux-S-varname"
            help = "S Qflux variable name in file"
            arg_type = String
            default  = ""


        "--MLD-varname"
            help = "MLD variable name in file"
            arg_type = String
            required = true

        "--MLD-average"
            help = "Make MLD constant as average."
            action = :store_true


        "--Tclim-varname"
            help = "Tclim variable name in file"
            arg_type = String
            default  = ""

        "--Sclim-varname"
            help = "Sclim variable name in file"
            arg_type = String
            default  = ""


        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

        "--output-file"
            help = "Output file name."
            arg_type = String
            required = true

        "--time"
            help = "Whether the 12 data points are the beginning of the month or the middle of the month."
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

nomissing = (x) -> replace(x, missing=>NaN)

dom = convert(Array{Float64}, [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
time = zeros(Float64, 12)

if parsed["time"] == "begin"
    time[1] = 0.0
    for t=2:12
        time[t] = time[t-1] + dom[t-1]
    end
elseif parsed["time"] == "mid"
    time[1] = dom[1] / 2.0
    for t=2:12
        time[t] = time[t-1] + (dom[t-1] + dom[t]) / 2.0
    end
else
    throw(ErrorException("Unknown keyword of `time`: " * parsed["time"]))
end
println("time: ", time)


Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)

end

empty = zeros(Float64, ni, nj, length(time))
for t=1:length(time)
    view(empty,:,:,t)[mask .== 0] .= NaN
end


Dataset(parsed["input-file"], "r") do ds
    global qflux_T, qflux_S, Tclim, Sclim

    if parsed["Qflux-blank"] 
        qflux_T = empty
        qflux_S = empty
    else
        qflux_T = ds[parsed["Qflux-T-varname"]][:] |> nomissing
        qflux_S = ds[parsed["Qflux-S-varname"]][:] |> nomissing
    end

    if parsed["Tclim-varname"] == "" 
        Tclim = empty
    else
        Tclim = ds[parsed["Tclim-varname"]][:] |> nomissing
    end

    if parsed["Sclim-varname"] == "" 
        Sclim = empty
    else
        Sclim = ds[parsed["Sclim-varname"]][:] |> nomissing
    end


    global h     = ds[parsed["MLD-varname"]][:] |> nomissing
    if parsed["MLD-average"] 
        for i=1:ni, j=1:nj
            h[i, j, :] .= mean(h[i, j, :])
        end
    end

end 

avg_qflux_T = 0.0
avg_qflux_S = 0.0
mask_idx = (mask .== 1)
for t=1:12
    if any( isnan.(qflux_T[:,:,t][mask_idx]) )
        throw(ErrorException("T Qflux contains NaN"))
    end

    if any( isnan.(qflux_S[:,:,t][mask_idx]) )
        throw(ErrorException("S Qflux contains NaN"))
    end

    global avg_qflux_T += sum( (area .* view(qflux_T, :, :, t))[mask_idx] )
    global avg_qflux_S += sum( (area .* view(qflux_S, :, :, t))[mask_idx] )
end
avg_qflux_T /= 12.0 * sum(area[mask_idx])
avg_qflux_S /= 12.0 * sum(area[mask_idx])

#qflux .-= avg_qflux

#println("# Pre-adjustment avg qflux : ", avg_qflux)

println("# Avg qflux_T : ", avg_qflux_T)
println("# Avg qflux_S : ", avg_qflux_S)



Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", length(time))

    ds.attrib["avg_qflux_T"] = avg_qflux_T
    ds.attrib["avg_qflux_S"] = avg_qflux_S

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

        ("qflux_T", qflux_T, ("ni", "nj", "time"), Dict(  
            "units"     => "W/m^2",
            "long_name" => "ocean heat flux correction",
        )),

        ("qflux_S", qflux_S, ("ni", "nj", "time"), Dict(  
            "units"     => "kg / s",
            "long_name" => "ocean salt flux correction",
        )),


        ("Tclim", Tclim, ("ni", "nj", "time"), Dict(  
            "units"     => "K",
            "long_name" => "Climatology of SST.",
        )),

        ("Sclim", Sclim, ("ni", "nj", "time"), Dict(  
            "units"     => "g/kg",
            "long_name" => "Climatology of SSS.",
        )),



        ("hblt", h, ("ni", "nj", "time"), Dict(
            "units"     => "m",
            "long_name" => "Mixed-layer depth",
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
