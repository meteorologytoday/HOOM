include("include_packages.jl")

include("OcnEnv.jl")
include("SharedData.jl")


ocn_env = OcnEnv(
    "grid.nc",
    "topo.nc",
    "bg_TS.nc",
    3600.0,
    24,
    8,
    8,
    360,
    180,
    20,
    3,
    collect(Float64, range(0.0, -1000.0, length=21)),
    [1, 5, 14],
    2,
    0,
    1000.0,
    1e-3,
    1e-3,
    [1e-3, 1e-3], 
    [1e-3, 1e-3], 
    10.0,
    1000.0,
    0.48,
    23.0,
)

sd = SharedData(ocn_env)

regVariable!(
    sd,
    :T,
    :fT,
    :xyz,
    Float64,
)
