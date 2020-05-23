using Interpolations
using NCDatasets
using ArgParse
using JSON
using Statistics
using Formatting
function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help     = "Input filename."
            arg_type = String
            required = true

        "--varname"
            help     = "Variable name."
            arg_type = String
            required = true

        "--domain-file"
            help = "Domain file. Must contain variables: xc, yc, area, mask"
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

dom = convert(Array{Int64}, [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
doy = sum(dom)
time = zeros(Float64, 12)
extended_time = zeros(Float64, 14)
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

extended_time[2:13] .= time
extended_time[1]     = time[end] - doy
extended_time[14]    = time[1]   + doy


println("time: ", time)
println("extended_time: ", extended_time)

Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, Nx, Ny
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    Nx, Ny = size(area)
end

Dataset(parsed["input-file"], "r") do ds
    global data = ds[parsed["varname"]][:] |> nomissing
end

data_interp = zeros(Float64, Nx, Ny, doy)

data_slice = zeros(14)
for i=1:Nx, j=1:Ny

    if mask[i, j] == 0
        data_interp[i, j, :] .= -999
        continue
    end

    data_slice[2:13] .= data[i, j, :]
    data_slice[1]     = data_slice[12 + 1]
    data_slice[14]    = data_slice[1  + 1]

    itp = LinearInterpolation( (extended_time, ), data_slice )
    
    for t=1:doy
        data_interp[i, j, t] = itp(t-0.5)
    end
end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", Nx)
    defDim(ds, "nj", Ny)
    defDim(ds, "time", doy)

#    ds.attrib["avg_Qflx_T"] = avg_Qflx_T
#    ds.attrib["avg_Qflx_S"] = avg_Qflx_S

    for (varname, vardata, vardim, attrib) in [
        ("time", collect(1:doy).-0.5, ("time",), Dict(
            "calendar"  => "noleap",
            "long_name" => "observation time",
            "units"     => "days since 0001-01-01 00:00:00",
        ) ),

        (parsed["varname"], data_interp, ("ni", "nj", "time"), Dict(  
            "units"     => "",
            "long_name" => "",
        )),
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
