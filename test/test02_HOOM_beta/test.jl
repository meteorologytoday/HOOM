using NCDatasets
include("HOOM/src/models/HOOM_beta/HOOM.jl")

domain_file = "domain.ocn_aqua.fv4x5_gx3v7.091218.nc"
Dataset(domain_file, "r") do ds
    global Ny = ds.dim["nj"]
end
ev = HOOM.Env(;
        id = 1,
        gf_filename = "domain.ocn_aqua.fv4x5_gx3v7.091218.nc",
        sub_yrng = 1:Ny,
        z_w      = collect(0.0:-10.0:-400.0),
        τ_TEMP   = nothing,
        τ_SALT   = nothing,
)
 
