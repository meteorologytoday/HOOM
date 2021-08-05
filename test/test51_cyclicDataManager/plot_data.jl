include("share/CyclicData.jl")

using Plots
using .CyclicData
using CFTime

cdm =  CyclicDataManager(;
    timetype     = DateTimeNoLeap,
    filename     = "paper2021_CTL_POP2.100years.nc",
    varnames     = ["HBLT", "TEMP", "SALT", "TAUX", "TAUY"],
    beg_time     = DateTimeNoLeap(0, 1, 1, 0, 0, 0),
    end_time     = DateTimeNoLeap(1, 1, 1, 0, 0, 0),
    align_time   = DateTimeNoLeap(0, 1, 1, 0, 0, 0),
)

t_vec = timedecode(collect(Float64, 0:2:365), "days since 0001-01-01 00:00:00", DateTimeNoLeap)

t_num_vec = timeencode(t_vec, "days since 0001-01-01 00:00:00", DateTimeNoLeap)


d_vec = Dict()
for v in ["HBLT", "TEMP", "SALT", "TAUX", "TAUY"]
    d_vec[v] = zeros(length(t_vec))
end

data = makeDataContainer(cdm)

for (i, t) in enumerate(t_vec)
    
    interpData!(cdm, t, data)

    for (k, v) in d_vec
            v[i] = data[k][1, 210, 250]
    end

end

p = []


for (k, v) in d_vec
    push!(p, plot(t_num_vec, d_vec[k], markershape=:circle, label=k))
end

plot(p...)
