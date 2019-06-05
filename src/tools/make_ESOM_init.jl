include("../ESOM/ESOM.jl")

using .ESOM
using NCDatasets
using Formatting


data_path = joinpath("../../../data")

# ===== [BEGIN] topo and climatology =====

Dataset("$data_path/domain.ocn.gx3v7.120323.nc", "r") do ds
    global Nx = ds.dim["ni"]
    global Ny = ds.dim["nj"]
    global mask = convert(Array{Float64}, replace(ds["mask"][:], missing=>NaN))
end

Dataset("$data_path/clim_LENS_B1850C5CN_005_gx3v7_TEMP.nc", "r") do ds
    global _Ts_clim = replace(ds["TEMP"][:, :, :, 1], missing=>NaN)
    z_w_top = replace(ds["z_w_top"][:], missing=>NaN) / 100.0
    z_w_bot = replace(ds["z_w_bot"][:], missing=>NaN) / 100.0

    global zs = zeros(Float64, length(z_w_top) + 1)
    zs[1:end-1] = z_w_top
    zs[end] = z_w_bot[end]

    zs .*= -1

end

Dataset("$data_path/clim_LENS_B1850C5CN_005_gx3v7_SALT.nc", "r") do ds
    global _Ss_clim = replace(ds["SALT"][:, :, :, 1], missing=>NaN)
end

Dataset("$data_path/ocean_topog_gx3v7.nc", "r") do ds
    global topo = zeros(Float64, 1, 1)
    topo = - replace(ds["depth"][:], missing => NaN)
end

dims = size(_Ts_clim)

bnds_1 = 1:4
bnds_2 = 4:27

layers_1 = bnds_1[1] : bnds_1[2] - 1
layers_2 = bnds_2[1] : bnds_2[2] - 1


zs_1 = zs[bnds_1]
zs_2 = zs[bnds_2]

hs_1 = zs_1[1:end-1] - zs_1[2:end]
hs_2 = zs_2[1:end-1] - zs_2[2:end]

Ts_clim = zeros(Float64, dims[1:2]..., 2)
Ss_clim = zeros(Float64, dims[1:2]..., 2)

Ts_clim .= NaN
Ss_clim .= NaN


sum_hs_1 = sum(hs_1)
sum_hs_2 = sum(hs_2)

for i=1:dims[1], j=1:dims[2]

    if mask[i, j] == 0
        continue
    end
    
    Ts_clim[i, j, 1] = sum(hs_1 .* _Ts_clim[i, j, layers_1]) / sum_hs_1
    Ts_clim[i, j, 2] = sum(hs_2 .* _Ts_clim[i, j, layers_2]) / sum_hs_2
 
    Ss_clim[i, j, 1] = sum(hs_1 .* _Ss_clim[i, j, layers_1]) / sum_hs_1
    Ss_clim[i, j, 2] = sum(hs_2 .* _Ss_clim[i, j, layers_2]) / sum_hs_2
   
    
    if isnan(Ts_clim[i, j, 1])
        Ts_clim[i, j, 1] = _Ts_clim[i, j, 1]
        Ts_clim[i, j, 2] = _Ts_clim[i, j, 1]
    elseif isnan(Ts_clim[i, j, 2])
        Ts_clim[i, j, 2] = Ts_clim[i, j, 1]
    end
 
    if isnan(Ss_clim[i, j, 1])
        Ss_clim[i, j, 1] = _Ss_clim[i, j, 1]
        Ss_clim[i, j, 2] = _Ss_clim[i, j, 1]
    elseif isnan(Ss_clim[i, j, 2])
        Ss_clim[i, j, 2] = Ss_clim[i, j, 1]
    end
     
end


Ts_clim .+= 273.15

# ===== [END] topo and climatology =====

occ = ESOM.OceanColumnCollection(
    gridinfo_file = "$data_path/domain.ocn.gx3v7.120323.nc",
    Nx       = Nx,
    Ny       = Ny,
    hs       = [50.0, 250.0],
    Ts       = Ts_clim,
    Ss       = Ss_clim,
    Kh_T     = 25000.0,
    Kh_S     = 25000.0,
    fs       = nothing,
    ϵs       = 1.0 / (86400.0 * 10), # 10 days
    mask     = mask,
    topo     = topo,
)

ESOM.takeSnapshot(occ, "ESOM_init_clim_LENS_B1850C5CN_005.nc")


