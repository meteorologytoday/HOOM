include(joinpath("..", "load_files.jl"))
include(joinpath(src, "NKOM", "NKOM.jl"))
using .NKOM


zs = zeros(Float64, 2)

zs[1] = 0.0
zs[2] = mininum(topo)

Ts_clim = Ts_clim[:, :, 1]
Ss_clim = Ss_clim[:, :, 1]

occ = NKOM.OceanColumnCollection(
    gridinfo_file = parsed["domain-file"],
    Nx       = Nx,
    Ny       = Ny,
    zs_bone  = zs,
    Ts       = Ts_clim,
    Ss       = Ss_clim,
    K_T      = 1e-5,
    K_S      = 1e-5,
    T_ML     = Ts_clim,
    S_ML     = Ss_clim,
    h_ML     = h_ML[:, :, 1], 
    h_ML_min = 0.0,
    h_ML_max = -zs[end],             # make it unrestricted
    we_max   = 1e5,
    mask     = mask,
    topo     = topo,
    Ts_clim_relax_time = 0.0,
    Ts_clim            = nothing,
    Ss_clim_relax_time = 0.0,
    Ss_clim            = nothing,
)

NKOM.takeSnapshot(occ, parsed["output-file"])


