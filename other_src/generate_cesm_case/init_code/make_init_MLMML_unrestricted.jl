include("load_files.jl")
include(joinpath(src, "MLMML", "MLMML.jl"))
using .MLMML

occ = MLMML.OceanColumnCollection(
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
    h_ML_max = 1e5,             # make it unrestricted
    we_max   = 1e5,             # make it unrestricted
    mask     = mask,
    topo     = topo,
    Ts_clim_relax_time = nothing, # Make it unrestricted
    Ts_clim            = nothing,
    Ss_clim_relax_time = nothing,
    Ss_clim            = nothing,
)

MLMML.takeSnapshot(occ, parsed["output-file"])


