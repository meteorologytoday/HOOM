using NCDatasets
using ArgParse
using JSON
using Statistics
using Formatting
function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-Qflx-file"
            help     = "If not specified, then output Qflx_S = Qflx_T = 0.0"
            arg_type = String
            default  = ""

        "--input-MLD-file"
            help     = "If not specified, then output MLD = 0.0"
            arg_type = String
            default  = ""

        "--input-clim-file"
            help     = "If not specified, then output T_clim = S_clim = 0.0"
            arg_type = String
            default  = ""

        "--domain-file"
            help = "Domain file. Must contain variables: xc, yc, area, mask"
            arg_type = String
            required = true

        "--Qflx_T-varname"
            help = "T Qflux variable name in Qflx file."
            arg_type = String
            default  = "Qflx_T"

        "--Qflx_S-varname"
            help = "S Qflux variable name in Qflx file."
            arg_type = String
            default  = "Qflx_S"

        "--T_clim-varname"
            help = "T_clim variable name in clim file."
            arg_type = String
            default  = "T_clim"

        "--S_clim-varname"
            help = "S_clim variable name in clim file."
            arg_type = String
            default  = "S_clim"

        "--MLD-varname"
            help = "MLD variable name in MLD file."
            arg_type = String
            default = "MLD"

        "--MLD-average"
            help = "Make MLD constant as average."
            action = :store_true

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

global Qflx_T = empty
global Qflx_S = empty


if parsed["input-Qflx-file"] != ""


    Dataset(parsed["input-Qflx-file"], "r") do ds
        if haskey(ds, parsed["Qflx_T-varname"])
            global Qflx_T = ds[parsed["Qflx_T-varname"]][:] |> nomissing
        else
            println(format("Variable {:s} does not exist. Skip it.", parsed["Qflx_T-varname"]))
        end

        if haskey(ds, parsed["Qflx_S-varname"])
            global Qflx_S = ds[parsed["Qflx_S-varname"]][:] |> nomissing
        else
            println(format("Variable {:s} does not exist. Skip it.", parsed["Qflx_S-varname"]))
        end

    end
end


if parsed["input-MLD-file"] == ""
    
    MLD = empty
    
else

    Dataset(parsed["input-MLD-file"], "r") do ds

        global MLD = ds[parsed["MLD-varname"]][:] |> nomissing
        if parsed["MLD-average"] 
            for i=1:ni, j=1:nj
                MLD[i, j, :] .= mean(MLD[i, j, :])
            end
        end

    end
 
end


if parsed["input-clim-file"] == ""

    T_clim = empty
    S_clim = empty

else

    Dataset(parsed["input-clim-file"], "r") do ds

        global T_clim, S_clim

        T_clim = ds[parsed["T_clim-varname"]][:] |> nomissing
        S_clim = ds[parsed["S_clim-varname"]][:] |> nomissing

    end
end

avg_Qflx_T = 0.0
avg_Qflx_S = 0.0
mask_idx = (mask .== 1)
for t=1:12
    if any( isnan.(Qflx_T[:,:,t][mask_idx]) )
        throw(ErrorException("T Qflux contains NaN"))
    end

    if any( isnan.(Qflx_S[:,:,t][mask_idx]) )
        throw(ErrorException("S Qflux contains NaN"))
    end

    global avg_Qflx_T += sum( (area .* view(Qflx_T, :, :, t))[mask_idx] )
    global avg_Qflx_S += sum( (area .* view(Qflx_S, :, :, t))[mask_idx] )
end
avg_Qflx_T /= 12.0 * sum(area[mask_idx])
avg_Qflx_S /= 12.0 * sum(area[mask_idx])

println("# Avg Qflx_T : ", avg_Qflx_T)
println("# Avg Qflx_S : ", avg_Qflx_S)

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", length(time))

    ds.attrib["avg_Qflx_T"] = avg_Qflx_T
    ds.attrib["avg_Qflx_S"] = avg_Qflx_S

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
            "long_name" => "Longitude of grid cell center",
        )),

        ("yc", yc, ("ni", "nj"), Dict(
            "units"     => "degrees north",
            "long_name" => "Latitude of grid cell center",
        )),

        ("Qflx_T", Qflx_T, ("ni", "nj", "time"), Dict(  
            "units"     => "W / m^2",
            "long_name" => "Ocean heat flux correction",
        )),

        ("Qflx_S", Qflx_S, ("ni", "nj", "time"), Dict(  
            "units"     => "kg / s / m^2",
            "long_name" => "Ocean salt flux correction",
        )),

        ("T_clim", T_clim, ("ni", "nj", "time"), Dict(  
            "units"     => "K",
            "long_name" => "Climatology of SST.",
        )),

        ("S_clim", S_clim, ("ni", "nj", "time"), Dict(  
            "units"     => "g / kg",
            "long_name" => "Climatology of SSS.",
        )),

        ("MLD", MLD, ("ni", "nj", "time"), Dict(
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
