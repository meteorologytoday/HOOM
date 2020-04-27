include("PUHSOM.jl")

Nz_f = 20
Nz_c = 3

hrgrid_file = "/seley/tienyiah/CESM_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/ocean_topog_gx1v6.nc"

ocn_env = PUHSOM.OcnEnv(
    hrgrid_file,
    topo_file,
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
    0.48,
    23.0,
    1e-2,
    [10.0, 1000.0],
    ["T.nc", "S.nc"],
    ["T", "S"],
    [NaN, NaN],
    "exponential_decay",
    true,
)


PUHSOM.init!(ocn_env)
