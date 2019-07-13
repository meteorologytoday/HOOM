include(joinpath("..", "load_files.jl"))
include(joinpath(src, "SOM", "SOM.jl"))
using .SOM


occ = SOM.OceanColumnCollection(
    Nx       = Nx,
    Ny       = Ny,
    Kh_T     = 0.0,
    T_ML     = Ts_clim[:, :, 1],
    h_ML     = 30.0,
    mask     = mask,
)

SOM.takeSnapshot(occ, parsed["output-file"])


