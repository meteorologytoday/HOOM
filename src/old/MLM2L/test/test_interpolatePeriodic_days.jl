include("../MLM2L.jl")



days_of_month = convert(Array{Float64}, [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
period = sum(days_of_month)

t = zeros(Float64, length(days_of_month)+1)
for i=2:length(t)
    t[i] = t[i-1] + days_of_month[i-1]
end

t_old = (t[2:end] + t[1:end-1]) / 2.0
t_new = collect(Float64, 1:Int(period)) .- 1

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
