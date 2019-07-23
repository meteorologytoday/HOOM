include("../MLM2L.jl")



period = 12.0
t_old = collect(Float64, 0:12)[1:end-1]
t_new = collect(Float64, range(0, stop=period, length=30))[1:end-1]

d_old = zeros(Float64, 3, length(t_old))

d_old[1,:] = sin.(t_old / period * 2π) 
d_old[2,:] = sin.(t_old / period * 2π) + cos.(t_old / period * 4π) 
d_old[3,:] = rand(length(t_old))


d_new = zeros(Float64, 3, length(t_new))


MLM2L.interpolatePeriodic!(
    time1 = t_old,
    time2 = t_new,
    data1 = d_old,
    data2 = d_new,
    period = period,
)

using PyPlot

fig, ax = plt[:subplots](3, 1, sharex=true)

for i=1:3
    ax[i][:scatter](t_old, d_old[i, :], s=40, marker="o", color="r")
    ax[i][:scatter](t_new, d_new[i, :], s=20, marker="x", color="k")
end

plt[:show]()
