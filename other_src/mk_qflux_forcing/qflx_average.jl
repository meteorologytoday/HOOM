using NCDatasets
using ArgParse
using JSON
using Formatting

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Qflux file"
            arg_type = String
            required = true

        "--domain-file"
            help = "Domain file."
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
    global qflux = ds["qdp"][:] |> nomissing
end 

idx = mask .== 1.0
q_avg = zeros(Float64, 12)

area_sum = sum(area[idx])
for t = 1:12
    q_avg[t] = sum((area .* qflux[:, :, t])[idx]) / area_sum
    println(format("Average Qflux for month {:02d} : {:.2e}", t, q_avg[t]))
end

println(format("Average Qflux all : {:.2e}", sum(q_avg)/12.0))
