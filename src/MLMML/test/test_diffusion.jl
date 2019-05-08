
include("test_profile_config.jl")

if do_calculation 

Δt = 86400.0
sim_steps = 360 * 20

t = collect(Float64, 0:sim_steps)

_fric_u  = zeros(Float64, 1, 1)
_swflx  = zeros(Float64, 1, 1)
_nswflx = zeros(Float64, 1, 1)
_frwflx = zeros(Float64, 1, 1)

rec = Dict(
    "T"  => zeros(Float64, occ.Nz[1, 1], length(t)),
    "S"  => zeros(Float64, occ.Nz[1, 1], length(t)),
)

init_profile = Dict(
    "T" => occ.Ts[1, 1, :],
    "S" => occ.Ss[1, 1, :],
)


for k = 1:length(t)
    if (k-1)%30 == 0
        print(format("\rIteration = {:d}/{:d} ({:.1f}%)", k, length(t), k / length(t) * 100.0))
    end

    if k != 0
        MLMML.stepOceanColumnCollection!(
            occ;
            fric_u = _fric_u,
            swflx  = _swflx,
            nswflx = _nswflx,
            frwflx = _frwflx,
            Δt     = Δt,
        )
    end

    rec["T"][:, k]   = occ.Ts[1, 1, :]
    rec["S"][:, k]   = occ.Ss[1, 1, :]

end

println()
println("===== Simulation is done =====")

#rec["T"] .-= repeat(rec["T"][:, end], outer=(1, length(t)))
#rec["S"] .-= repeat(rec["S"][:, end], outer=(1, length(t)))



end


using PyPlot


# Fig 1: Single lines

plot_interval = Int(sim_steps / 10)

mid_zs = (zs[1:end-1] + zs[2:end]) / 2.0

fig, ax = plt[:subplots](1, 2, figsize=(6, 4))

for i = 1:plot_interval:length(t)
    ax[1][:plot](rec["T"][:,i]      , mid_zs)
    ax[2][:plot](rec["S"][:,i] * 1e3, mid_zs)
end

ax[1][:set_title]("Init T [K]")
ax[2][:set_title]("Init S [g/Kg]")


# Fig 2: Hovmoeller Diagram

fig, (ax1, ax2) = plt[:subplots](2, 1, figsize=(12, 8))

# Temperature profile
cmap = plt[:get_cmap]("jet")
cbmapping = ax1[:contourf](t, mid_zs, rec["T"], cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, ax=ax1)
cb[:set_label]("Temperature [\$ \\mathrm{K} \$]")

# Salinity profile
cmap = plt[:get_cmap]("jet")
cbmapping = ax2[:contourf](t, mid_zs, rec["S"] * 1e3, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, ax=ax2)
cb[:set_label]("Salinity [\$ g/Kg \$]")


#=
#convadj_rec[convadj_rec .== 0.0] .= NaN
#ax2[:scatter](t_day, -600.0 * convadj_rec, marker="^")
#ax2[:plot](t_day, - h_rec , "r--", linewidth=2, zorder=10)

ax1[:set_ylabel]("SST anomaly\n[\$\\mathrm{K}\$]")
ax2[:set_ylabel]("Z [m]")
ax2[:set_xlabel]("Time [year]")

ticks      = collect(0:5:YEARS_WANTED)
ticklabels = [format("{:d}", ticks[i]) for i=1:length(ticks)]

ax1[:set_xticks](ticks)
ax2[:set_xticks](ticks)

ax1[:set_xticklabels](ticklabels)
ax2[:set_xticklabels](ticklabels)


ax2[:set_ylim]([-1500, 0])
ax2[:set_ylim]([-1100, 0])

tlim = [t_yr[1], t_yr[end]]


ax1[:set_ylim]([-0.5, 0.5])
ax1[:set_xlim](tlim)
ax2[:set_xlim](tlim)


fig, ax = plt[:subplots](1, 2, figsize=(6, 4))

for i = 1:plot_interval:length(t)
    ax[1][:plot](rec["T"][:,i]      , mid_zs)
    ax[2][:plot](rec["S"][:,i] * 1e3, mid_zs)
end

ax[1][:set_title]("Init T [K]")
ax[2][:set_title]("Init S [g/Kg]")

=#
