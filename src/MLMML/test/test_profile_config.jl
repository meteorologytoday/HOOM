include("../MLMML.jl")
using .MLMML
using NCDatasets
using Formatting

"""
zs = collect(Float64, range(0.0, stop=100.0 , length=51))
push!(zs, collect(Float64, range(100, step=10, stop=2000))[2:end]...)
push!(zs, collect(Float64, range(2000, step=50, stop=4000))[2:end]...)
zs *= -1.0
"""
# ===== [BEGIN] topo and climatology =====

pick = (i=65, j=60)


Dataset("clim_LENS_B1850C5CN_005_gx3v7_zMLMML_TEMP.nc", "r") do ds
    global Ts_clim = replace(ds["TEMP"][pick.i, pick.j, :, 1], missing=>NaN)
    println(ds["lat"][pick.i, pick.j])
    println(ds["lon"][pick.i, pick.j])

    global zs = replace(ds["zs"][:], missing=>NaN)
end

Dataset("clim_LENS_B1850C5CN_005_gx3v7_zMLMML_SALT.nc", "r") do ds
    global Ss_clim = replace(ds["SALT"][pick.i, pick.j, :, 1], missing=>NaN)
end


println(zs)

Dataset("ocean_topog_gx3v7.nc", "r") do ds
    global topo = zeros(Float64, 1, 1)
    println(ds["depth"][pick.i, pick.j])
    topo[1, 1] = - ds["depth"][pick.i, pick.j]
end

# ===== [END] topo and climatology =====



N = length(zs) - 1

ΔT_0 = 5.0
ΔS_0 = 0.0
T_ML_0 = 300.0
S_ML_0 = 35e-3
h_ML_0 = 10.0
h_ML_min = 10.0
h_ML_max = 1000.0
we_max = 1.0

T_slope = 2.0 / 4000.0
S_slope = 0.0



if test_type == "ALL"
    occ = MLMML.makeBasicOceanColumnCollection(
        Nx       = 1,
        Ny       = 1,
        zs_bone  = zs,
        T_slope  = T_slope,
        S_slope  = S_slope,
        T_ML     = T_ML_0,
        S_ML     = S_ML_0,
        h_ML     = h_ML_0,
        ΔT       = ΔT_0,
        ΔS       = ΔS_0,
        K_T      = 1e-5,
        K_S      = 1e-5,
        h_ML_min = h_ML_min,
        h_ML_max = h_ML_max,
        we_max   = we_max,
        topo     = topo,
        Ts_clim_relax_time = 86400.0 * 360 * 10,
        Ts_clim  = Ts_clim,
        Ss_clim_relax_time = 86400.0 * 360 * 10,
        Ss_clim  = Ss_clim,
    )

elseif test_type == "TOPO"

    occ = MLMML.makeBasicOceanColumnCollection(
        Nx       = 1,
        Ny       = 1,
        zs_bone  = zs,
        T_slope  = T_slope,
        S_slope  = S_slope,
        T_ML     = T_ML_0,
        S_ML     = S_ML_0,
        h_ML     = h_ML_0,
        ΔT       = ΔT_0,
        ΔS       = ΔS_0,
        K_T      = 0.0,
        K_S      = 0.0,
        h_ML_min = h_ML_min,
        h_ML_max = h_ML_max,
        we_max   = we_max,
        topo     = topo,
    )

elseif test_type == "DIFFUSION"
    occ = MLMML.makeBasicOceanColumnCollection(
        Nx       = 1,
        Ny       = 1,
        zs_bone  = zs,
        T_slope  = T_slope,
        S_slope  = S_slope,
        T_ML     = T_ML_0,
        S_ML     = S_ML_0,
        h_ML     = h_ML_0,
        ΔT       = ΔT_0,
        ΔS       = ΔS_0,
        K_T      = 1e-5,
        K_S      = 1e-5,
        h_ML_min = h_ML_min,
        h_ML_max = h_ML_max,
        we_max   = we_max,
    )

elseif test_type == "CLIM"
    occ = MLMML.makeBasicOceanColumnCollection(
        Nx       = 1,
        Ny       = 1,
        zs_bone  = zs,
        T_slope  = T_slope,
        S_slope  = S_slope,
        T_ML     = T_ML_0,
        S_ML     = S_ML_0,
        h_ML     = h_ML_0,
        ΔT       = ΔT_0,
        ΔS       = ΔS_0,
        K_T      = 0.0,
        K_S      = 0.0,
        h_ML_min = h_ML_min,
        h_ML_max = h_ML_max,
        we_max   = we_max,
        Ts_clim_relax_time = 86400.0 * 360 * 10,
        Ts_clim  = Ts_clim,
        Ss_clim_relax_time = 86400.0 * 360 * 10,
        Ss_clim  = Ss_clim,
    )
end
Nz = occ.Nz_bone
