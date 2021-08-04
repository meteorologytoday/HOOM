include("share/CyclicData.jl")

using Plots
using .CyclicData

cdm =  CyclicDataManager(;
    filename     = "paper2021_CTL_POP2.100years.nc",
    varname_time = "time",
    varnames     = ["HBLT", "TEMP", "SALT", "TAUX", "TAUY"],
    beg_time     = 0.0,
    cyc_time     = 365.0,
)

#t_vec = cdm.t_vec#collect(Float64, 0:5:365)
t_vec = collect(Float64, 0:2:365)


d_vec = Dict()

for v in ["HBLT", "TEMP", "SALT", "TAUX", "TAUY"]
    d_vec[v] = t_vec * 0.0
end

data = interpData!(cdm, 0.0 ; create = true)

for (i, t) in enumerate(t_vec)
    
    interpData!(cdm, t, data)

    for (k, v) in d_vec
            v[i] = data[k][1, 210, 250]
    end

end

p = []


for (k, v) in d_vec
    push!(p, plot(t_vec, d_vec[k], markershape=:circle, label=k))
end

plot(p...)
