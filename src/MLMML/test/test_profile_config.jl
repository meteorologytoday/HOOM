include("../MLMML.jl")

zs = collect(Float64, range(0.0, stop=100.0 , length=51))
push!(zs, collect(Float64, range(100, step=10, stop=2000))[2:end]...)
push!(zs, collect(Float64, range(2000, step=50, stop=4000))[2:end]...)
zs *= -1.0

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
    K_T      = 1e-3,
    K_S      = 1e-3,
    h_ML_min = h_ML_min,
    h_ML_max = h_ML_max,
    we_max   = we_max,
)










