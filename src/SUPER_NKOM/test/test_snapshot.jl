include("../MLMML.jl")
include("../../CESM_driver/julia_lib/NetCDFIO.jl")
include("../../share/constants.jl")

using .NetCDFIO
using .MLMML

domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
map = NetCDFIO.MapInfo{Float64}(domain_file)
zs_bone = collect(Float64, range(0, -500, step=-5))

T_ML     = T_ref
S_ML     = S_ref
ΔT       = 1.0
ΔS       = 0.0
T_slope  = 30.0 / 5000.0
S_slope  = 0.0
K_T      = 1e-5
K_S      = 1e-5
h_ML_min = 10.0
h_ML_max = 500.0
we_max   = 1.0

mask = nothing
topo = nothing
    
occ = MLMML.makeBasicOceanColumnCollection(
    Nx      = map.nx,
    Ny      = map.ny,
    zs_bone = zs_bone,
    T_slope = T_slope,
    S_slope = S_slope,
    T_ML    = T_ML,
    S_ML    = S_ML,
    h_ML    = h_ML_min,
    ΔT      = ΔT,
    ΔS      = ΔS,
    K_T     = K_T,
    K_S     = K_S,
    h_ML_min= h_ML_min,
    h_ML_max= h_ML_max,
    we_max  = we_max,
    mask    = mask,
    topo    = topo,
)

MLMML.takeSnapshot(occ, "snapshot01.nc")
occ2 = MLMML.loadSnapshot("snapshot01.nc")

MLMML.takeSnapshot(occ, "snapshot02.nc")
