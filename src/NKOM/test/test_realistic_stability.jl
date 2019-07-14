include("LinearRegression.jl")

using Printf
using Statistics: mean
using Formatting
using StatsBase

@printf("Importing PyPlot... ")
using PyCall
plt = pyimport("matplotlib.pyplot")
GS  = pyimport("matplotlib.gridspec")
@printf("done.\n")

test_type = "ALL"

include("test_profile_config.jl")

if(do_calculation)



SECS_PER_DAY = 86400.0
DAYS_PER_MON = 30
MONS_PER_YEAR= 12
DAYS_PER_YEAR = DAYS_PER_MON * MONS_PER_YEAR
SECS_PER_YEAR = DAYS_PER_YEAR * SECS_PER_DAY

SPINUP_YEARS = 20
YEARS_WANTED = 100
TOTAL_YEARS = SPINUP_YEARS + YEARS_WANTED

TOTAL_DAYS = TOTAL_YEARS * DAYS_PER_YEAR
TOTAL_SECS = TOTAL_DAYS * SECS_PER_DAY

SPINUP_DAYS = SPINUP_YEARS * DAYS_PER_YEAR
DAYS_WANTED = YEARS_WANTED * DAYS_PER_YEAR
MONS_WANTED = YEARS_WANTED * MONS_PER_YEAR

ω = 2π/360.0/86400.0
t_sim = collect(Float64, range(0.0, step=SECS_PER_DAY, stop=TOTAL_SECS))[1:end-1]
Δt = t_sim[2] - t_sim[1]

t = t_sim[SPINUP_DAYS+1:end]
t .-= t[1]


swflx  = 125.0 * cos.(ω*t_sim)
nswflx = swflx * 0.0

frwflx = 0*(40.0 / 1000.0 / 86400.0) * (- cos.(ω*t_sim))



U10 = zeros(Float64, length(t_sim))
U10 .= 8.0 .+ 2.0 * cos.(ω*t_sim)


using PyPlot
mid_zs = (zs[1:end-1] + zs[2:end]) / 2.0
fig, ax = plt[:subplots](1, 3, sharey=true)
ax[1][:plot](occ.bs[1,1,:], mid_zs)
ax[2][:plot](occ.Ts[1,1,:], mid_zs)
ax[3][:plot](occ.Ss[1,1,:] * 1e3, mid_zs)

ax[1][:set_title]("Init b [m/s^2]")
ax[2][:set_title]("Init T [K]")
ax[3][:set_title]("Init S [g/Kg]")



rec = Dict(
    "b"  => zeros(Float64, Nz, DAYS_WANTED),
    "T"  => zeros(Float64, Nz, DAYS_WANTED),
    "S"  => zeros(Float64, Nz, DAYS_WANTED),
    "h"  => zeros(Float64, DAYS_WANTED),
    "hb" => zeros(Float64, DAYS_WANTED),
    "swflx" => zeros(Float64, DAYS_WANTED),
    "nswflx" => zeros(Float64, DAYS_WANTED),
    "frwflx" => zeros(Float64, DAYS_WANTED),
)

fric_u  = zeros(Float64, 1, 1)
_swflx  = zeros(Float64, 1, 1)
_nswflx = zeros(Float64, 1, 1)
_frwflx = zeros(Float64, 1, 1)


for k = 1:length(t_sim)
    if true || (k-1)%DAYS_PER_YEAR == 0
        print(format("\rIteration = {:d}/{:d} ({:.1f}%)", k, length(t_sim), k / length(t_sim) * 100.0))
    end

    if all(isfinite.(occ.Ts))
        println("At simulation")
        throw(ErrorException(occ.Ts[1, 1, 65:68] |> string))
    end


    fric_u[1,1]  = NKOM.getFricU(ua=U10[k])
    _swflx[1,1]  = swflx[k]
    _nswflx[1,1] = nswflx[k]
    _frwflx[1,1] = frwflx[k]

    if k != 0
        NKOM.stepOceanColumnCollection!(
            occ;
            fric_u = fric_u,
            swflx  = _swflx,
            nswflx = _nswflx,
            frwflx = _frwflx,
            Δt     = Δt,
        )
    end

    if k > SPINUP_DAYS

        i = k - SPINUP_DAYS
        rec["swflx"][i]   = swflx[k]
        rec["nswflx"][i]  = nswflx[k]
        rec["frwflx"][i]  = frwflx[k]
        rec["h"][i]       = occ.h_ML[1, 1]
        #rec["b"][i]       = occ.b_ML[1, 1]
        rec["hb"][i]      = NKOM.OC_getIntegratedBuoyancy(occ, 1, 1)
        rec["b"][:, i]    = occ.bs[1, 1, :]
        rec["T"][:, i]    = occ.Ts[1, 1, :]
        rec["S"][:, i]    = occ.Ss[1, 1, :]

        #println(rec["b"][65:68])
    end
end

println()
using Statistics


# Doing Dayily RMSC
day = rec
day_rmsc = Dict()
day_means  = Dict()     # means  of every layer of every day
day_len = length(t)

for k in ("b", "T", "S")
    
    day_means[k]    = zeros(Nz, DAYS_PER_YEAR)
    day_rmsc[k]     = zeros(Nz, day_len)
    for i=1:Nz
        d_timeseries = day[k][i, :]

        if any(isnan.(d_timeseries))
            break
        end

        # Because detrending is not uniform across different layers,
        # here I just simply average them to get avg profile
        day_means[k][i, :] = mean( reshape( d_timeseries, DAYS_PER_YEAR, :), dims=2)[:,1]
        β = LinearRegression(t, d_timeseries)
        d_timeseries -= β[1] .+ β[2] * t
       
        seasonality = repeat(
            mean( reshape( d_timeseries, DAYS_PER_YEAR, :), dims=2)[:,1]
            , outer=(YEARS_WANTED,)
        )
        day_rmsc[k][i, :] = d_timeseries - seasonality
    end
end




# Doing Monthly Average
mon_len = Int(length(t)/DAYS_PER_MON)
mon = Dict(
    "t" => zeros(mon_len),
    "b" => zeros(Nz, mon_len),
    "T" => zeros(Nz, mon_len),
    "S" => zeros(Nz, mon_len),
    "h" => zeros(mon_len),
)

mon_rmsc = Dict()
total_trends = Dict()     # trends of every layer
mon_means  = Dict()     # means  of every layer and every month


for i = 1:mon_len
    mon_rng = (1+DAYS_PER_MON*(i-1)):DAYS_PER_MON*i
    mon["t"][i] = mean(t[mon_rng])
    mon["h"][i] = mean(rec["h"][mon_rng])
    mon["b"][:, i] = mean(rec["b"][:, mon_rng], dims=2)
    mon["T"][:, i] = mean(rec["T"][:, mon_rng], dims=2)
    mon["S"][:, i] = mean(rec["S"][:, mon_rng], dims=2)
end



mon["t"] .-= mon["t"][1]


for k in ("b", "T", "S")
    total_trends[k] = zeros(Nz)
    mon_means[k]    = zeros(Nz, MONS_PER_YEAR)
    mon_rmsc[k]     = zeros(Nz, mon_len)

    for i=1:Nz
        d_timeseries = mon[k][i, :]

        # Because detrending is not uniform across different layers,
        # here I just simply average them to get avg profile
        mon_means[k][i, :] = mean( reshape( d_timeseries, MONS_PER_YEAR, :), dims=2)[:,1]

        β = LinearRegression(mon["t"], d_timeseries)
        d_timeseries -= β[1] .+ β[2] * mon["t"]
        total_trends[k][i] = β[2]
       
        seasonality = repeat(
            mean( reshape( d_timeseries, MONS_PER_YEAR, :), dims=2)[:,1]
            , outer=(YEARS_WANTED,)
        )
        mon_rmsc[k][i, :] = d_timeseries - seasonality
    end
end


end

println(format("Average of total trends: {:.2e} K / yr", mean(total_trends["T"][1:Nz]) * SECS_PER_YEAR))
println(format("Average of total trends: {:.2e} g/Kg / yr", mean(total_trends["S"][1:Nz]) * SECS_PER_YEAR))


t_yr = mon["t"] / SECS_PER_YEAR
mid_zs = (zs[1:end-1] + zs[2:end]) / 2.0

#=
auto_SST = autocor(SST_rec, collect(0:120); demean=true) 
auto_SST_yr = collect(Float64, 0:length(auto_SST)-1) / 12.0

plt[:figure]()
plt[:plot](auto_SST_yr, auto_SST)
=#

# Hovmoller diagram
gs0 = GS.GridSpec(1, 2, width_ratios=[100,5])
gs_l = GS.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], height_ratios=[1, 4])

fig = plt[:figure](figsize=(12, 6))
ax1 = plt[:subplot](gs_l[1])
ax2 = plt[:subplot](gs_l[2])
cax = plt[:subplot](gs0[2])


# SST record
ax1[:plot]([t_yr[1], t_yr[end]], [0, 0], "k--")
ax1[:plot](t_yr, mon_rmsc["T"][1, :], label="SST")
#ax1[:plot](t / SECS_PER_YEAR, day_rmsc["T"][1, :], label="SST")


cmap = plt[:get_cmap]("jet")

clevs = - [0.5, 0.2]
append!(clevs, collect(range(-0.1, stop=0.0, length=6)))
append!(clevs, -clevs[end-1:-1:1])

clevs = (range(-0.5, stop=0.5, length=21) |> collect )  
#cbmapping = ax2[:contourf](t_yr, mid_zs, mon_rmsc["T"] * 10.0, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cbmapping = ax2[:contourf](t/SECS_PER_YEAR, mid_zs, day_rmsc["T"] * 10.0, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, cax=cax)

cb[:set_label]("Temperature anomaly [\$ \\times 10^{-1} \\, \\mathrm{K} \$]")

#ax2[:plot](t_yr, - mon["h"], "r--")
ax2[:plot](t / SECS_PER_YEAR, - rec["h"], "r--")


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


using Formatting
fig[:suptitle](
    format(
        "Temperature anomaly (annual cycle removed)\nΔt = 1 day, spin up time = {:d} years", SPINUP_YEARS)
    )
plt[:show]()


# Monthly structure
fig, ax = plt[:subplots](1, 3, figsize=(15,6), sharey=true)

fig[:suptitle]("Vertical profile of monthly mean buoyancy")

for i = 1:size(mon_means["T"])[2]

    offset =  i
    ax[1][:plot](mon_means["T"][:, i], mid_zs, label="$i")
    ax[1][:text](mon_means["T"][1, i], zs[1]+offset, "$i", va="bottom", ha="center")

end

ax[1][:set_title]("Original")

ax[1][:set_ylim]([-200, 20])
ax[1][:legend]()

ax[1][:set_xlabel]("Temperature anomaly [\$\\mathrm{K}\$]")


ax[2][:plot]([0, 0], [-200, 20], "--", color="#888888")
ax[2][:plot](total_trends["T"] * SECS_PER_YEAR, mid_zs, "r-", label="trend")
ax[2][:set_title]("Trend for each layer")

ax[3][:plot](mean(mon_means["T"], dims=2)[:, 1], mid_zs, "r-", label="mean")
ax[3][:set_title]("Mean for each layer")

ax[2][:set_xlabel]("Trend of Temperature [\$\\mathrm{K}\\,\\mathrm{yr}^{-1} \$]")
ax[3][:set_xlabel]("Mean Temperature [\$\\mathrm{K}\$]")
plt[:show]()


# ===== MONTHLY STRUCTURE =====

fig, ax = plt[:subplots](1, 3, figsize=(15,6), sharey=true)

fig[:suptitle]("Monthly mean structure")

for i = 1:size(mon_means["S"])[2]
    offset =  i 
    ax[1][:plot](mon_means["b"][:, i], mid_zs, label="$i")
    ax[1][:text](mon_means["b"][1, i], zs[1]+offset, "$i", va="bottom", ha="center")

    ax[2][:plot](mon_means["T"][:, i] .- NKOM.T_ref, mid_zs, label="$i")
    ax[2][:text](mon_means["T"][1, i] .- NKOM.T_ref, zs[1]+offset, "$i", va="bottom", ha="center")

    ax[3][:plot]((mon_means["S"][:, i].- NKOM.S_ref) * 1e3, mid_zs, label="$i")
    ax[3][:text]((mon_means["S"][1, i].- NKOM.S_ref) * 1e3, zs[1]+offset, "$i", va="bottom", ha="center")
end

ax[1][:set_title]("Mean b")
ax[2][:set_title]("Mean T")
ax[3][:set_title]("Mean S")

ax[1][:set_xlabel]("b [m / s^2]")
ax[2][:set_xlabel]("T - T_ref [K]")
ax[3][:set_xlabel]("S - S_ref [g/Kg]")


ax[1][:set_ylim]([-200, 20])
plt[:show]()
