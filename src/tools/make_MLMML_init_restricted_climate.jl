include("../MLMML/MLMML.jl")

using .MLMML
using NCDatasets
using Formatting

data_path = joinpath("../../../data")

# ===== [BEGIN] topo and climatology =====

Dataset("$data_path/domain.ocn.gx3v7.120323.nc", "r") do ds
    global Nx = ds.dim["ni"]
    global Ny = ds.dim["nj"]
    global mask = convert(Array{Float64}, replace(ds["mask"][:], missing=>NaN))
end

Dataset("$data_path/clim_LENS_B1850C5CN_005_gx3v7_zMLMML_TEMP.nc", "r") do ds
    global Ts_clim = replace(ds["TEMP"][:, :, :, 1], missing=>NaN)
    global zs = replace(ds["zs"][:], missing=>NaN)
end

Dataset("$data_path/clim_LENS_B1850C5CN_005_gx3v7_zMLMML_SALT.nc", "r") do ds
    global Ss_clim = replace(ds["SALT"][:, :, :, 1], missing=>NaN)
end

Dataset("$data_path/ocean_topog_gx3v7.nc", "r") do ds
    global topo = zeros(Float64, 1, 1)
    topo = - replace(ds["depth"][:], missing => NaN)
end


Ts_clim .+= 273.15

# ===== [END] topo and climatology =====

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
    h_ML_max = 1e5,
    we_max   = 1e5,
    mask     = mask,
    topo     = topo,
    Ts_clim_relax_time = 86400.0 * 365 * 10,
    Ts_clim  = Ts_clim,
    Ss_clim_relax_time = 86400.0 * 365 * 10,
    Ss_clim  = Ss_clim,
)

MLMML.takeSnapshot(occ, "SSM_NK_init_clim_LENS_B1850C5CN_005_restricted_climate.nc")


