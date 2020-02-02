using NCDatasets
using Formatting
using ArgParse
using Statistics
using JSON

function mreplace(x)
    return replace(x, missing=>NaN)
end

function integrate(x, dydx)
    y = copy(dydx) * 0.0

    for i = 2:length(x)
        y[i] = y[i-1] + (dydx[i-1] + dydx[i]) * (x[i] - x[i-1]) / 2.0
    end

    return y
end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file"
            help = "Ocean data file."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--t-idx"
            help = "Time slot"
            arg_type = Int64
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

Dataset(parsed["domain-file"], "r") do ds

    global mask = ds["mask"][:] |> mreplace
    global area = ds["area"][:] |> mreplace

    global lnd_idx = (mask .!= 1.0)

    area[lnd_idx] .= 0.0

    global area_sum = sum(area)
end

Dataset(parsed["data-file"], "r") do ds
 
    rng = (:,:,parsed["t-idx"]) 


    function getWeightedData(varname)
        d = replace(ds[varname][rng...], missing=>0.0, NaN=>0.0)
        println(rng)
        println(sum(d))
        return sum(area .* d) / area_sum
    end

    global neb  = getWeightedData("neb")
    global dHdt = getWeightedData("dHdt")

    global total_energy_transport = neb - dHdt
end

using Formatting
println(format("# Total Area : {:e}", area_sum))
println(format("# Total Energy Transport : {:e}", total_energy_transport))
println(format("# Neb  : {:e}", neb))
println(format("# dHdt : {:e}", dHdt))
