include("../MLMML.jl")
include("../../CESM_driver/julia_lib/NetCDFIO.jl")

using .NetCDFIO
using .MLMML

domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
map = NetCDFIO.MapInfo{Float64}(domain_file)
zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5



init_b_ML     = 280.0 * MLMML.g * MLMML.α
init_h_ML     = MLMML.h_ML_min
init_b_slope  = 30.0 / 5000.0 * MLMML.g * MLMML.α
init_Δb       = 1.0 * MLMML.g * MLMML.α


occ = MLMML.makeBasicOceanColumnCollection(
    map.nx, map.ny, zs;
    b_slope = init_b_slope,
    b_ML  = init_b_ML,
    h_ML  = init_h_ML,
    Δb    = init_Δb,
    K     = K,
    mask  = map.mask,
)

MLMML.takeSnapshot(occ, "snapshot01.nc")
occ2 = MLMML.loadSnapshot("snapshot01.nc")

MLMML.takeSnapshot(occ, "snapshot02.nc")
