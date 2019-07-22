include(joinpath("..", "load_files.jl"))
include(joinpath(src, "NKOM", "NKOM.jl"))
using .NKOM

occ = NKOM.OceanColumnCollection(
    gridinfo_file = parsed["domain-file"],
    Nx       = Nx,
    Ny       = Ny,
    zs_bone  = zs,
    Ts       = copy(Ts_clim),
    Ss       = copy(Ss_clim),
    K_T      = 1e-5,
    K_S      = 1e-5,
    T_ML     = Ts_clim[:, :, 1],
    S_ML     = Ss_clim[:, :, 1],
    h_ML     = 10.0, 
    h_ML_min = 10.0,
    h_ML_max = 1e5,
    we_max   = 1e-2,
    mask     = mask,
    topo     = topo,
    Ts_clim_relax_time = 86400.0 * 365 * 10, # 10 years
    Ts_clim            = copy(Ts_clim),
    Ss_clim_relax_time = 86400.0 * 365 * 10, # 10 years
    Ss_clim            = copy(Ss_clim),
    arrange  = "xyz",
)

NKOM.takeSnapshot(occ, parsed["output-file"])


