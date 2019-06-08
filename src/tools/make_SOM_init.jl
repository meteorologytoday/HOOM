include("../SOM/SOM.jl")

using .SOM
using NCDatasets
using Formatting

data_path = joinpath("../../../data")

# ===== [BEGIN] topo and climatology =====

Dataset("$data_path/domain.ocn.gx3v7.120323.nc", "r") do ds
    global Nx = ds.dim["ni"]
    global Ny = ds.dim["nj"]
    global mask = convert(Array{Float64}, replace(ds["mask"][:], missing=>NaN))
end

# ===== [END] topo and climatology =====

occ = SOM.OceanColumnCollection(
    Nx       = Nx,
    Ny       = Ny,
    Kh_T     = 0.0,
    T_ML     = 15.0 + 273.15,
    h_ML     = 30.0,
    mask = mask,
)

SOM.takeSnapshot(occ, "$data_path/SSM_SOM_init.nc")


