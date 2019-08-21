println("The test currently includes: TOPO, DIFFUSION, CLIM, ALL ")

if length(ARGS) > 1
    test_type = ARGS[1]
end

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
    "T"  => zeros(Float64, occ.Nz_bone, length(t)),
    "S"  => zeros(Float64, occ.Nz_bone, length(t)),
)

init_profile = Dict(
    "T" => occ.Ts[1, 1, :],
    "S" => occ.Ss[1, 1, :],
)

#println("init_T")
#println(init_profile["T"])

for k = 1:length(t)
    if (k-1)%30 == 0
        print(format("\rIteration = {:d}/{:d} ({:.1f}%)", k, length(t), k / length(t) * 100.0))
    end

    if k != 1
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


