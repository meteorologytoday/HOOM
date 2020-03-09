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

Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)i

end


let

    global Qflx_T_correction = zeros(Float64, ni, nj, 12)
    global Qflx_S_correction = zeros(Float64, ni, nj, 12)

    global TSAS = zeros(Float64, ni, nj)
    global SSAS = zeros(Float64, ni, nj)


    for m=1:12

        println(format("Doing month : {:02d}", m))

        for y=parsed["beg-year"]:parsed["end-year"]
            Dataset(format("{:s}{:04d}-{:02d}.nc", parsed["data-file-prefix"], y, m)) do ds
                Qflx_T_correction[:, :, m] +=  nomissing( ds["qflx_T_correction"][:, :, 1] )
                Qflx_S_correction[:, :, m] +=  nomissing( ds["qflx_S_correction"][:, :, 1] )
                TSAS +=  nomissing( ds["TSAS_clim"][:, :, 1] )
                SSAS +=  nomissing( ds["SSAS_clim"][:, :, 1] )
            end
        end


    end
   
    nyears = parsed["end-year"] - parsed["beg-year"] + 1 
    Qflx_T_correction /= nyears
    Qflx_S_correction /= nyears
    TSAS /= nyears * 12
    SSAS /= nyears * 12

    TSAS *= 3996.0 * 1026.0

    # Apply mask
    idx = mask .== 0.0
    for m=1:12
        view(Qflx_T_correction, :, :, m)[idx] .= NaN
        view(Qflx_S_correction, :, :, m)[idx] .= NaN
    end

    println("Size of TSAS", size(TSAS))

    TSAS[idx] .= NaN
    SSAS[idx] .= NaN

    idx = mask .== 1.0
    _area = area[idx]
    _area_sum = sum(_area)
    global avg_Qflx_T_correction = sum(_area .* (mean(Qflx_T_correction, dims=3)[:, :, 1])[idx]) / _area_sum
    global avg_Qflx_S_correction = sum(_area .* (mean(Qflx_S_correction, dims=3)[:, :, 1])[idx]) / _area_sum
    global avg_TSAS = sum(_area .* TSAS[idx]) / _area_sum
    global avg_SSAS = sum(_area .* SSAS[idx]) / _area_sum
end

println("avg_Qflx_T_correction: ", avg_Qflx_T_correction )
println("avg_Qflx_S_correction: ", avg_Qflx_S_correction )
println("avg_TSAS: ", avg_TSAS )
println("avg_SSAS: ", avg_SSAS )

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", Inf)

    ds.attrib["avg_Qflx_T_correction"] = avg_Qflx_T_correction
    ds.attrib["avg_Qflx_S_correction"] = avg_Qflx_S_correction
    ds.attrib["avg_TSAS"] = avg_TSAS
    ds.attrib["avg_SSAS"] = avg_SSAS

    for (varname, vardata, vardim, attrib) in [
        ("Qflx_T_correction", Qflx_T_correction, ("ni", "nj", "time"), Dict()),
        ("Qflx_S_correction", Qflx_S_correction, ("ni", "nj", "time"), Dict()),
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



