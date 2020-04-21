include("include_packages.jl")

include("OcnEnv.jl")
include("SharedData.jl")
include("JobDistributionInfo.jl")

Nx = 360
Ny = 180
Nz_f = 20
Nz_c = 3

ocn_env = OcnEnv(
    "grid.nc",
    "topo.nc",
    "bg_TS.nc",
    3600.0,
    24,
    8,
    8,
    Nx,
    Ny,
    Nz_f,
    Nz_c,
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

# create new variable
regVariable!(
    sd,
    :T,
    :fT,
    :xyz,
    Float64,
)

# add existing variable
S_cT = SharedArray{Float64}(Nx, Ny, Nz_c)
regVariable!(
    sd,
    :S_cT,
    :cT,
    :xyz,
    Float64;
    data=S_cT
)


jdi = JobDistributionInfo(env=ocn_env)
