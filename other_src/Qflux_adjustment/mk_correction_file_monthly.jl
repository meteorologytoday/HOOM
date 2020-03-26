using NCDatasets
using ArgParse
using JSON
using Statistics
using Formatting

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file-prefix"
            help = "Data filename prefix including folder and path until the timestamp. File extension `nc` is assumed."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            default = ""

        "--input-file"
            help = "Input Qflx file name."
            arg_type = String
            required = true

        "--output-file"
            help = "Output file name."
            arg_type = String
            required = true

        "--beg-year"
            help = "Begin year of the data."
            arg_type = Int64
            required = true

        "--end-year"
            help = "End year of the data."
            arg_type = Int64
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

nomissing = (x) -> replace(x, missing=>NaN)

if parsed["domain-file"] == ""
    parsed["domain-file"] = parsed["input-file"]
end

Dataset(parsed["domain-file"], "r") do ds

    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)



end

idx       = mask .== 1.0
_area     = area[idx]
_area_sum = sum(_area)

calAreaAvg = ( field , ) -> sum(_area .* (field[idx])) / _area_sum
function calMonthlyAreaAvg(field)
    return convert(Array{Float64}, [calAreaAvg(view(field, :, :, m)) for m=1:12])
end

Dataset(parsed["input-file"], "r") do ds

    global Qflx_T_old, Qflx_S_old
    
    Qflx_T_old = ds["Qflx_T"][:] |> nomissing
    Qflx_S_old = ds["Qflx_S"][:] |> nomissing

end


let

    global Qflx_T_correction = zeros(Float64, ni, nj, 12)
    global Qflx_S_correction = zeros(Float64, ni, nj, 12)

    global TSAS = zeros(Float64, ni, nj, 12)
    global SSAS = zeros(Float64, ni, nj, 12)


    for m=1:12

        println(format("Doing month : {:02d}", m))

        for y=parsed["beg-year"]:parsed["end-year"]
            Dataset(format("{:s}{:04d}-{:02d}.nc", parsed["data-file-prefix"], y, m)) do ds
                Qflx_T_correction[:, :, m] +=  nomissing( ds["qflx_T_correction"][:, :, 1] )
                Qflx_S_correction[:, :, m] +=  nomissing( ds["qflx_S_correction"][:, :, 1] )
                TSAS[:, :, m]              +=  nomissing( ds["TSAS_clim"][:, :, 1] )
                SSAS[:, :, m]              +=  nomissing( ds["SSAS_clim"][:, :, 1] )
            end
        end


    end
   
    nyears = parsed["end-year"] - parsed["beg-year"] + 1 
    Qflx_T_correction /= nyears
    Qflx_S_correction /= nyears
    TSAS /= nyears
    SSAS /= nyears

    # Transform into energy flux and salt flux
    TSAS *= 3996.0 * 1026.0

    #=
    # Apply mask
    idx = mask .== 0.0
    for m=1:12
        view(Qflx_T_correction, :, :, m)[idx] .= NaN
        view(Qflx_S_correction, :, :, m)[idx] .= NaN
    end

    println("Size of TSAS", size(TSAS))

    TSAS[idx] .= NaN
    SSAS[idx] .= NaN
    =#
end

# calculate average information
avg_Qflx_T_old = Qflx_T_old |> calMonthlyAreaAvg
avg_Qflx_S_old = Qflx_S_old |> calMonthlyAreaAvg

avg_Qflx_T_correction = Qflx_T_correction |> calMonthlyAreaAvg
avg_Qflx_S_correction = Qflx_S_correction |> calMonthlyAreaAvg

avg_TSAS = TSAS |> calMonthlyAreaAvg
avg_SSAS = SSAS |> calMonthlyAreaAvg

avg_Qflx_T_new = avg_Qflx_T_old + avg_Qflx_T_correction
avg_Qflx_S_new = avg_Qflx_S_old + avg_Qflx_S_correction

avg_Qflx_T_imbalance = avg_Qflx_T_new - (-avg_TSAS)
avg_Qflx_S_imbalance = avg_Qflx_S_new - (-avg_SSAS)

# new_2 is adjusted such that new_2 + SAS = 0
# avg_Qflx_T_new -= avg_Qflx_T_imbalance
# avg_Qflx_S_new -= avg_Qflx_S_imbalance

Qflx_T_new_before = Qflx_T_old + Qflx_T_correction
Qflx_S_new_before = Qflx_S_old + Qflx_S_correction

Qflx_T_new = copy(Qflx_T_new_before)
Qflx_S_new = copy(Qflx_S_new_before)
for m=1:12
    Qflx_T_new[:, :, m] .-= avg_Qflx_T_imbalance[m]
    Qflx_S_new[:, :, m] .-= avg_Qflx_S_imbalance[m]
end

avg_Qflx_T_new_before = Qflx_T_new_before |> calMonthlyAreaAvg
avg_Qflx_S_new_before = Qflx_S_new_before |> calMonthlyAreaAvg

avg_Qflx_T_new = Qflx_T_new |> calMonthlyAreaAvg
avg_Qflx_S_new = Qflx_S_new |> calMonthlyAreaAvg

avg_Qflx_T_imbalance_new = avg_Qflx_T_new + avg_TSAS
avg_Qflx_S_imbalance_new = avg_Qflx_S_new + avg_SSAS

for m = 1:12
    println("Month: ", m)
    println("avg_Qflx_T_old: ", avg_Qflx_T_old[m] )
    println("avg_Qflx_S_old: ", avg_Qflx_S_old[m] )
    println()
    println("avg_Qflx_T_correction: ", avg_Qflx_T_correction[m] )
    println("avg_Qflx_S_correction: ", avg_Qflx_S_correction[m] )
    println()
    println("avg_TSAS       [ W/m^2]: ", avg_TSAS[m] )
    println("avg_SSAS       [kg/m^2]: ", avg_SSAS[m] )
    println()
    println("avg_Qflx_T_imbalance: ", avg_Qflx_T_imbalance[m] )
    println("avg_Qflx_S_imbalance: ", avg_Qflx_S_imbalance[m] )
    println()
    println("avg_Qflx_T_new: ", avg_Qflx_T_new[m] )
    println("avg_Qflx_S_new: ", avg_Qflx_S_new[m] )
    println()
    println("avg_Qflx_T_imbalance_new: ", avg_Qflx_T_imbalance_new[m])
    println("avg_Qflx_S_imbalance_new: ", avg_Qflx_S_imbalance_new[m])
end

println("========================")
println("avg_Qflx_T_imbalance: ", mean(avg_Qflx_T_imbalance) )
println("avg_Qflx_S_imbalance: ", mean(avg_Qflx_S_imbalance) )
println("========================")


Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", 12)

    ds.attrib["AVG_Qflx_T_old"] = mean(avg_Qflx_T_old)
    ds.attrib["AVG_Qflx_S_old"] = mean(avg_Qflx_S_old)
    ds.attrib["AVG_TSAS"] = mean(avg_TSAS)
    ds.attrib["AVG_SSAS"] = mean(avg_SSAS)
    ds.attrib["AVG_Qflx_T_new_before"] = mean(avg_Qflx_T_new_before)
    ds.attrib["AVG_Qflx_S_new_before"] = mean(avg_Qflx_S_new_before)
    ds.attrib["AVG_Qflx_T_new"] = mean(avg_Qflx_T_new)
    ds.attrib["AVG_Qflx_S_new"] = mean(avg_Qflx_S_new)
    ds.attrib["AVG_Qflx_T_imbalance"] = mean(avg_Qflx_T_imbalance)
    ds.attrib["AVG_Qflx_S_imbalance"] = mean(avg_Qflx_S_imbalance)
    ds.attrib["AVG_Qflx_T_imbalance_new"] = mean(avg_Qflx_T_imbalance_new)
    ds.attrib["AVG_Qflx_S_imbalance_new"] = mean(avg_Qflx_S_imbalance_new)


    for (varname, vardata, vardim, attrib) in [
        ("Qflx_T_old", Qflx_T_old, ("ni", "nj", "time"), Dict()),
        ("Qflx_S_old", Qflx_S_old, ("ni", "nj", "time"), Dict()),
        ("TSAS", TSAS, ("ni", "nj", "time"), Dict()),
        ("SSAS", SSAS, ("ni", "nj", "time"), Dict()),
        ("Qflx_T_correction", Qflx_T_correction, ("ni", "nj", "time"), Dict()),
        ("Qflx_S_correction", Qflx_S_correction, ("ni", "nj", "time"), Dict()),
        ("Qflx_T_new", Qflx_T_new, ("ni", "nj", "time"), Dict()),
        ("Qflx_S_new", Qflx_S_new, ("ni", "nj", "time"), Dict()),
        ("area", area, ("ni", "nj"), Dict()),
        ("mask", mask, ("ni", "nj"), Dict()),
        ("avg_Qflx_T_old", avg_Qflx_T_old, ("time",), Dict()),
        ("avg_Qflx_S_old", avg_Qflx_S_old, ("time",), Dict()),
        ("avg_TSAS", avg_TSAS, ("time",), Dict()),
        ("avg_SSAS", avg_SSAS, ("time",), Dict()),
        ("avg_Qflx_T_new", avg_Qflx_T_new, ("time",), Dict()),
        ("avg_Qflx_S_new", avg_Qflx_S_new, ("time",), Dict()),
        ("avg_Qflx_T_imbalance_new", avg_Qflx_T_imbalance_new, ("time",), Dict()),
        ("avg_Qflx_S_imbalance_new", avg_Qflx_T_imbalance_new, ("time",), Dict()),
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



