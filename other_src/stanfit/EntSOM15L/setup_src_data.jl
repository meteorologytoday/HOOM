using NCDatasets

src_data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data", "NCAR_LENS")

F_filename = joinpath(src_data_path,   "transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SHF.100001-109912.nc")
SST_filename = joinpath(src_data_path, "transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SST.100001-109912.nc")
MLD_filename = joinpath(src_data_path, "transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.HMXL.100001-109912.nc")
T_filename = joinpath(src_data_path, "transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.HMXL.100001-109912.nc")

Dataset(F_filename, "r") do ds
    global nlon = ds.dim["nlon"]
    global nlat = ds.dim["nlat"]
end



