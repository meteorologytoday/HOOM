include("../MLMML.jl")
include("../../CESM_driver/julia_lib/NetCDFIO.jl")
include("../../share/constants.jl")

using .NetCDFIO
using .MLMML

domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
map = NetCDFIO.MapInfo{Float64}(domain_file)
zs = collect(Float64, range(0, -500, step=-5))

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
    
    Nx      :: Integer,
    Ny      :: Integer,
    zs_bone :: AbstractArray{Float64, 1};
    T_slope :: Float64,
    S_slope :: Float64,
    T_ML    :: Float64,
    S_ML    :: Float64,
    h_ML    :: Float64,
    ΔT      :: Float64,
    ΔS      :: Float64,
    K_T     :: Float64,
    K_S     :: Float64,
    h_ML_min:: Float64,
    h_ML_max:: Float64,
    we_max  :: Float64,
    mask    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    topo    :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,


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
