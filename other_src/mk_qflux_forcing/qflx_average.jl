using NCDatasets
using ArgParse
using JSON
using Formatting
using Statistics

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Qflux file"
            arg_type = String
            required = true
        
    end

    return parse_args(ARGS, s)
end


parsed = parse_commandline()
print(json(parsed, 4))

nomissing = (x) -> replace(x, missing=>NaN)

Dataset(parsed["input-file"], "r") do ds

    global mask, xc, yc, area, ni, nj, qflux

    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing

    ni, nj = size(area)
    qflx_T = ds["Qflx_T"][:] |> nomissing
    qflx_S = ds["Qflx_S"][:] |> nomissing
    
end 

idx = mask .== 1.0
qflx_T_avg = zeros(Float64, 12)
qflx_S_avg = zeros(Float64, 12)

area_sum = sum(area[idx])
for t = 1:12
    qflx_T_avg[t] = sum((area .* qflx_T[:, :, t])[idx]) / area_sum
    qflx_S_avg[t] = sum((area .* qflx_S[:, :, t])[idx]) / area_sum
    println(format("Average Qflx_T for month {:02d} : {:.2e}", t, qflx_T_avg[t]))
    println(format("Average Qflx_S for month {:02d} : {:.2e}", t, qflx_S_avg[t]))
end

println(format("Average Qflx_T all : {:.2e}", mean(qflx_T_avg)))
println(format("Average Qflx_S all : {:.2e}", mean(qflx_S_avg)))

