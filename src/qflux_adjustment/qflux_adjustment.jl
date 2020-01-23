using NCDatasets
using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Qflux file"
            arg_type = String
            required = true

        "--output-file"
            help = "Output file name."
            arg_type = String
            required = true

        "--target-SST-file"
            help = "Target SST file name."
            arg_type = String
            required = true

        "--current-SST-file"
            help = "Target SST file name."
            arg_type = String
            required = true

        "--sensitivity-file"
            help = "SST sensitivity file name."
            arg_type = String
            default = "sensitivity.nc"

        "--default-sensitivity"
            help     = "Default climate sensitivity in unit of K / (W/m^2)"
            arg_type = Float64
            default  = - 0.67

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

nomissing = (x) -> replace(x, missing=>NaN)

Dataset(parsed["input-file"], "r") do ds

    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:]  |> nomissing
    xc   = ds["xc"][:]    |> nomissing
    yc   = ds["yc"][:]    |> nomissing
    area = ds["area"][:]  |> nomissing
    ni, nj = size(area)
    
    global qflux, h
    qflux = ds["qdp"][:]  |> nomissing
    h     = ds["hblt"][:] |> nomissing

    #global qflux_time
    #qflux_time = timeencode(ds["time"][:] |> nomissing, "days since 0001-01-01 00:00:00", calendar="noleap")

    #if length(qflux_time) != 12
    #    throw(ErrorException("qflux needs to be exactly 12 months"))
    #end

end 

Dataset(parsed["target-SST-file"], "r") do ds

    global target_SST
    target_SST = ds["SST"][:,:,1,:]  |> nomissing

end 

Dataset(parsed["current-SST-file"], "r") do ds

    global current_SST
#    current_SST = ds["SST"][:]  |> nomissing
    current_SST = ds["T_ML"][:]  |> nomissing

end 

println(size(current_SST))
println(size(target_SST))

SST_offset = current_SST - target_SST


dom = convert(Array{Float64}, [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
qflux_time = zeros(Float64, 12)

qflux_time[1] = dom[1] / 2.0
for t=2:12
    qflux_time[t] = qflux_time[t-1] + (dom[t-1] + dom[t]) / 2.0
end

# Reading or creating sensitivity file
if !isfile(parsed["sensitivity-file"])

    Dataset(parsed["sensitivity-file"], "c") do ds
        defDim(ds, "ni", ni)
        defDim(ds, "nj", nj)
        defDim(ds, "time", 12)
        defDim(ds, "trial", Inf)

        var = defVar(ds, "qflux_record", Float64, ("ni", "nj", "time", "trial"))
        var.attrib["_FillValue"] = 1e20
#        var[:, :, :, 1] = qflux
 
        var = defVar(ds, "SST_offset", Float64, ("ni", "nj", "time", "trial"))
        var.attrib["_FillValue"] = 1e20
#        var[:, :, :, 1] = SST_offset

        var = defVar(ds, "crude_sensitivity", Float64, ("ni", "nj", "time", "trial"))
        var.attrib["_FillValue"] = 1e20
 
    end

end

# Append latest result
Dataset(parsed["sensitivity-file"], "a") do ds

    local trial = ds.dim["trial"] + 1

    ds["qflux_record"][:, :, :, trial] = qflux
    ds["SST_offset"][:, :, :, trial]   = SST_offset 

end

# Calculate updated qflux
Dataset(parsed["sensitivity-file"], "a") do ds

    _, _, _, trial = size(ds["SST_offset"])
    S = zeros(Float64, ni, nj, 12)    

    println("trial: ", trial)

    if trial == 1
        S .= parsed["default-sensitivity"]
    else

        qfluxs      = ds["qflux_record"][:, :, :, trial-1:trial]      |> nomissing
        SST_offsets = ds["SST_offset"][:, :, :, trial-1:trial]  |> nomissing
        
        for i=1:ni, j=1:nj, t=1:12
            S[i, j, t] = ( SST_offsets[i, j, t, 2] - SST_offsets[i, j, t, 1] ) / ( qfluxs[i, j, t, 2] - qfluxs[i, j, t, 1] )
        end
            
    end

    if any(S .== 0.0)
        throw(ErrorException("Sensitivity 0 derived !!"))
    end
    
    ds["crude_sensitivity"][:, :, :, trial] = S

    global new_qflux = qflux .- SST_offset ./ S

    valid_idx = mask .== 1

    for t=1:12
        new_qflux[:, :, t] .-= sum(new_qflux[:, :, t][valid_idx] .* area[valid_idx]) / sum(area[valid_idx])
    end

end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", length(qflux_time))

    empty = zeros(Float64, ni, nj, length(qflux_time))

    for (varname, vardata, vardim, attrib) in [

        ("time", qflux_time, ("time",), Dict(
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

        ("qdp", new_qflux, ("ni", "nj", "time"), Dict(  
            "units"     => "W/m^2",
            "long_name" => "ocean heat flux convergence",
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
