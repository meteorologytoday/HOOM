include(joinpath("..", "load_files.jl"))
include(joinpath(src, "NKOM", "NKOM.jl"))
using .NKOM


Δz_1 =  50.0 
Δz_2 = 250.0


Dataset(parsed["forcing-file"], "r") do ds

    h_ML = replace(ds["hblt"][:], missing=>NaN)

    global h_ML_max = maximum(h_ML[isfinite.(h_ML)])
    global replace    



end


for i = 1:length(zs)
    if 0.0 - zs[i] >= Δz_1
        global bnd_1 = i
        break
    end
end

for i = bnd_1:length(zs)
    if zs[bnd_1]- zs[i] >= Δz_2
        global bnd_2 = i
        break
    end
end

bnds_1 =     1:bnd_1
bnds_2 = bnd_1:bnd_2

layers_1 = bnds_1[1] : bnds_1[2] - 1
layers_2 = bnds_2[1] : bnds_2[2] - 1


zs_1 = zs[bnds_1]
zs_2 = zs[bnds_2]

hs_1 = zs_1[1:end-1] - zs_1[2:end]
hs_2 = zs_2[1:end-1] - zs_2[2:end]

_Ts_clim = zeros(Float64, Nx, Ny, 2)
_Ss_clim = zeros(Float64, Nx, Ny, 2)

_Ts_clim .= NaN
_Ss_clim .= NaN


sum_hs_1 = sum(hs_1)
sum_hs_2 = sum(hs_2)

for i=1:Nx, j=1:Ny

    if mask[i, j] == 0
        continue
    end
    
    _Ts_clim[i, j, 1] = sum(hs_1 .* Ts_clim[i, j, layers_1]) / sum_hs_1
    _Ts_clim[i, j, 2] = sum(hs_2 .* Ts_clim[i, j, layers_2]) / sum_hs_2
 
    _Ss_clim[i, j, 1] = sum(hs_1 .* Ss_clim[i, j, layers_1]) / sum_hs_1
    _Ss_clim[i, j, 2] = sum(hs_2 .* Ss_clim[i, j, layers_2]) / sum_hs_2
   
    # Might touch the bottom of ocean
    if isnan(_Ts_clim[i, j, 1])
        _Ts_clim[i, j, 1] = Ts_clim[i, j, 1]
        _Ts_clim[i, j, 2] = Ts_clim[i, j, 1]
    elseif isnan(Ts_clim[i, j, 2])
        _Ts_clim[i, j, 2] = Ts_clim[i, j, 1]
    end
 
    if isnan(_Ss_clim[i, j, 1])
        _Ss_clim[i, j, 1] = Ss_clim[i, j, 1]
        _Ss_clim[i, j, 2] = Ss_clim[i, j, 1]
    elseif isnan(Ss_clim[i, j, 2])
        _Ss_clim[i, j, 2] = Ss_clim[i, j, 1]
    end
     
end

# Check if weird
for i=1:Nx, j=1:Ny
    if mask[i, j] != 0
        if any(isnan.(_Ts_clim[i, j, :])) || any(isnan.(_Ss_clim[i, j, :]))
            throw(ErrorException("Some data are missing at ", i, ", ", j))
        end
    end
end
occ = ESOM.OceanColumnCollection(
    gridinfo_file = parsed["domain-file"],
    Nx       = Nx,
    Ny       = Ny,
    hs       = [Δz_1, Δz_2],
    Ts       = _Ts_clim,
    Ss       = _Ss_clim,
    Kh_T     = 25000.0,
    Kh_S     = 25000.0,
    fs       = nothing,
    ϵs       = 1e-5,    # 1 day
    mask     = mask,
    topo     = topo,
)

ESOM.takeSnapshot(occ, parsed["output-file"])




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
    h_ML_max = 1e5,             # make it unrestricted
    we_max   = 1e-2,             # make it unrestricted
    mask     = mask,
    topo     = topo,
    Ts_clim_relax_time = 86400.0 * 365 * 10, # 10 years
    Ts_clim            = copy(Ts_clim),
    Ss_clim_relax_time = 86400.0 * 365 * 10, # 10 years
    Ss_clim            = copy(Ss_clim),
    arrange  = "xyz",
)

NKOM.takeSnapshot(occ, parsed["output-file"])


