
using Printf
using Statistics: mean
using Formatting
using NCDatasets

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

filename = joinpath(dirname(@__FILE__), "..", "..", "data", "Gaussian_SSM_output_002.nc")
println(filename)
vars = ["mld", "sst", "sumflx", "fric_u"]

avg_vars = Dict()

ds = Dataset(filename, "r")

lat = ds["lat"]
time = 1:ds.dim["time"] |> collect
for var in vars
    println("Doing average of var: ", var)
    v = nomissing(ds[var][:], NaN)
    println(size(v))
    cnt = sum(isfinite.(v), dims=1)

    v[isnan.(v)] .= 0
    avg_vars[var] = (sum(v, dims=1) ./ cnt)[1, :, :]
    
end



gs0 = GS.GridSpec(4, 2, width_ratios=[100,5])
#gs_l = GS.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], height_ratios=[1, 4])

fig = plt[:figure](figsize=(8, 12))
ax1 = plt[:subplot](gs0[1,1])
ax2 = plt[:subplot](gs0[2,1])
ax3 = plt[:subplot](gs0[3,1])
ax4 = plt[:subplot](gs0[4,1])


cax1 = plt[:subplot](gs0[1,2])
cax2 = plt[:subplot](gs0[2,2])
cax3 = plt[:subplot](gs0[3,2])
cax4 = plt[:subplot](gs0[4,2])


cmap = plt[:get_cmap]("coolwarm")
clevs = range(0, stop=20, step=2) |> collect
cbmapping = ax1[:contourf](time, lat, avg_vars["sst"] .- 273.15, clevs, cmap=cmap, extend="both")
cb = plt[:colorbar](cbmapping, cax=cax1)
cb[:set_label]("SST [\$\\mathrm{K}\$]")
ax1[:set_ylabel]("Lat [deg]")

cmap = plt[:get_cmap]("binary")
clevs = range(0, stop=500, step=50) |> collect
cbmapping = ax2[:contourf](time, lat, avg_vars["mld"], clevs, cmap=cmap, extend="max")
cb = plt[:colorbar](cbmapping, cax=cax2)
cb[:set_label]("MLD [\$\\mathrm{m}\$]")
ax2[:set_ylabel]("Lat [deg]")

cmap = plt[:get_cmap]("binary")
clevs = range(0, stop=20, step=4) |> collect
cbmapping = ax3[:contourf](time, lat, avg_vars["fric_u"] * 1e2, 10, #=clevs,=# cmap=cmap, extend="max")
cb = plt[:colorbar](cbmapping, cax=cax3)
cb[:set_label]("Friction velocity [\$\\mathrm{cm} \\, \\mathrm{s}^{-1}\$]")
ax3[:set_ylabel]("Lat [deg]")

cmap = plt[:get_cmap]("bwr")
clevs = range(-200, stop=200, step=20) |> collect
cbmapping = ax4[:contourf](time, lat, avg_vars["sumflx"], clevs, cmap=cmap, extend="both")
cb = plt[:colorbar](cbmapping, cax=cax4)
cb[:set_label]("Total heat flux [\$\\mathrm{W} \\, \\mathrm{m}^{-1}\$]")
ax4[:set_ylabel]("Lat [deg]")
ax4[:set_xlabel]("Time [day]")

fig[:suptitle]("Second year run")

fig[:savefig]("diag.png", dpi=100)
#=
cmap = plt[:get_cmap]("jet")
clevs = (range(-1, stop=1, length=51) |> collect ) * 20
cbmapping = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec * 1e4, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, cax=cax)

cb[:set_label]("Buoyancy anomaly [\$\\times\\,10^{-3}\\,\\mathrm{m} \\, \\mathrm{s}^{-2}\$]")

ax2[:plot](t_day, - h_rec , "r--", linewidth=2, zorder=10)
ax2[:set_ylim]([-D, 0])

tlim = [t_day[1], t_day[end]]
ax1[:set_xlim](tlim)
ax2[:set_xlim](tlim)


ax1[:set_ylabel]("Insolation Flux [\$\\mathrm{m}^{2}\\, \\mathrm{s}^{-3}\$]")
ax2[:set_ylabel]("Z [m]")
ax2[:set_xlabel]("Time [year]")

xticks      = collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0
xticklabels = [format("{:d}", i-1) for i=1:length(xticks)]

ax1[:set_xticks](collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0)  
ax2[:set_xticks](collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0)

ax1[:set_xticklabels](xticklabels)
ax2[:set_xticklabels](xticklabels)

using Formatting
fig[:suptitle](format("Buoyancy anomaly (annual cycle removed) with Δt = 1 day, Δz = {:.1f}", zs[1]-zs[2]))
p
=#
