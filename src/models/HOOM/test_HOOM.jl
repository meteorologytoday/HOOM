include("HOOM.jl")

Nz_f = 20
Nz_c = 3

hrgrid_file = "/seley/tienyiah/CESM_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/ocean_topog_gx1v6.nc"

ocn_env = HOOM.OcnEnv(
    hrgrid_file,
    topo_file,
    "bg_TS.nc",
    3600.0,
    24,
    8,
    Nz_f,
    Nz_c,
    collect(Float64, range(0.0, -1000.0, length=21)),
    [1, 5, 14],
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


model = HOOM.Model()

HOOM.init!(model; ocn_env=ocn_env)
