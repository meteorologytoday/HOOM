include("WeightGeneration.jl")

using .WeightGeneration

using Formatting

in_filename  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc"
out_filename = "test_SHF.nc"
wgt_filename = "wgt_gx1v6_to_gx3v7.nc"

WeightGeneration.convertFile(in_filename, out_filename, wgt_filename, varnames2D=("SHF",))
