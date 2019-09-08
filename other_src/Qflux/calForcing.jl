include("../../src/share/DisplacedPoleCoordinate.jl")
include("../../src/share/MapInfo.jl")

using Formatting
using NCDatasets
using Statistics
ρ    = 1026.0  # kg / m^3
c_p  = 3996.0  # J / kg / K

dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
sum(dom) == 365 || throw(ErrorException("Sum of dom is not 365"))

Δts = (dom + circshift(dom, -1))/2.0 * 86400.0

println("Δts = ", Δts)

in_SST  = "lowres_SST.nc"
in_SHF  = "lowres_SHF.nc" 
in_TAUX = "lowres_TAUX.nc"
in_TAUY = "lowres_TAUY.nc"
in_HMXL = "lowres_HMXL.nc"
in_TEMP = "lowres_TEMP.nc"

in_zcoord = "b.e11.B1850C5CN.f45_g37.005.pop.h.TEMP.100001-109912.nc"

out_file = "forcing.gx3v7.nc"

missing_value = 1e20

function getTEMPProfile(t::Integer)

    local TEMP, BOT_TEMP
    Dataset(in_TEMP, "r") do ds
        TEMP = replace(ds["TEMP"][:, :, :, t], missing=>NaN)

        # Construct the last value of temperature
        BOT_TEMP = TEMP[:, :, 1]

        for i=1:Nx, j=1:Ny
            isnan(BOT_TEMP[i, j]) && continue
            for k=Nz:-1:1
                if isfinite(TEMP[i, j, k])
                    BOT_TEMP[i, j] = TEMP[i, j, k]
                    break
                end
            end
        end


    end
    return TEMP, BOT_TEMP
end


Dataset(in_SST, "r") do ds
    global SST = replace(ds["SST"][:, :, 1, :], missing=>NaN)
    global Nx, Ny, Nt = size(SST)

    (Nt%12 == 0) || throw(ErrorException("Time is not multiple of 12"))

    global nyears = Int(Nt / 12)

end

Dataset(in_SHF, "r") do ds
    global SHF = replace(ds["SHF"][:], missing=>NaN)
end


Dataset(in_HMXL, "r") do ds
    global HMXL = replace(ds["HMXL"][:], missing=>NaN) / 100.0
end

Dataset(in_TAUX, "r") do ds
    global TAUX = replace(ds["TAUX"][:], missing=>NaN)
end

Dataset(in_TAUY, "r") do ds
    global TAUY = replace(ds["TAUY"][:], missing=>NaN)
end


Dataset(in_zcoord, "r") do ds
    z_top = - replace(ds["z_w_top"][:], missing=>NaN) / 100.0
    z_bot = - replace(ds["z_w_bot"][:], missing=>NaN) / 100.0
    
    global zs = [z_top..., z_bot[end]]
    global hs = z_top - z_bot
    global Nz = length(zs) - 1

end

function findZInd(z)
    (z == 0.0) && (return 1)
    for k=1:Nz
        (zs[k] > z >= zs[k+1]) && (return k)
    end
end

function getTd(z, i, j, temp, bot_temp)
    T = temp[i, j, findZInd(z)]
    return ( isfinite(T) ) ? T : bot_temp[i, j]
end

function integrate(zb, zt, zs, Ts)
    ( zb > zt ) && throw(ErrorException("zb > zt error. zb = " * string(zb) * ", zt = " * string(zt)))

    kb = findZInd(zb)
    kt = findZInd(zt)

    if kb == kt
        return (zt - zb) * Ts[kb]
    else

        s = (zt - zs[kt+1]) * Ts[kt]

        for k = kt+1:kb-1
            s += hs[k] * Ts[k]
        end

        s += (zs[kb] - zb) * Ts[kb]

    end

end

mi = ModelMap.MapInfo{Float64}(gridinfo_file)
gi = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv; angle_unit=:deg)

Qs = zeros(Float64, Nx, Ny, 12)


dTdts = zeros(Float64, Nx, Ny, 12)
Fs    = zeros(Float64, Nx, Ny, 12)
Ents  = zeros(Float64, Nx, Ny, 12)
Qs    = zeros(Float64, Nx, Ny, 12)
hs    = zeros(Float64, Nx, Ny, 12) 
EKs   = zeros(Float64, Nx, Ny, 12) 

τx = zeros(Float64, Nx, Ny)
τy = zeros(Float64, Nx, Ny)

u      = zeros(Float64, Nx, Ny)
v      = zeros(Float64, Nx, Ny)
DIV    = zeros(Float64, Nx, Ny)

uT     = zeros(Float64, Nx, Ny)
vT     = zeros(Float64, Nx, Ny)
DIV_UT = zeros(Float64, Nx, Ny)

wT_bot = zeros(Float64, Nx, Ny)


for t = 13:length(T) - 12

    global TEMP, BOT_TEMP

    # month
    m = (t-1)%12+1

    # Loading Temperature profile
    TEMP_0, BOT_TEMP_0 = getTEMPProfile(t  )
    TEMP_1, BOT_TEMP_1 = getTEMPProfile(t+1)

    # Transform input wind stress vector first
    
    DisplacedPoleCoordinate.project!(
        gi,
        (view(TAUX, :, :, t) + view(TAUX, :, :, t+1))/2.0,
        (view(TAUY, :, :, t) + view(TAUY, :, :, t+1))/2.0,
        τx, τy,
        direction=:Forward,
    )

    # Calculate Ekman flow
    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue

        h = (HMXL[i, j, t] + HMXL[i, j, t+1] ) / 2.0
        s̃ = ϵs[i, j] + fs[i, j] * im
        H̃ = √(ocn.K_v / s̃)
        H = abs(H̃)
        
        M̃ = (τx[i, j] + τy[i, j] * im) / (ρ * s̃)
        ṽ_ek =   M̃ / H_ek
        u_ek, v_ek = real(ṽ_ek), imag(ṽ_ek)

        SST_mean = (SST[i, j, t] + SST[i, j, t+1])/2.0

        u[i, j] = u_ek
        v[i, j] = v_ek

        uT[i, j] = u_ek * SST_mean
        vT[i, j] = v_ek * SST_mean

    end
   
    DisplacedPoleCoordinate.DIV!(gi, u, v,  DIV, gi.mask)
    DisplacedPoleCoordinate.DIV!(gi, uT, vT,  DIV_UT, gi.mask)

    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue
        
        h_mean   = (HMXL[i, j, t] + HMXL[i, j, t+1] ) / 2.0
        Td_mean = ( getTd(-h[i, j, t+1], i, j, TEMP_1, BOT_TEMP_1) +  getTd(-h[i, j, t], i, j, TEMP_0, BOT_TEMP_0) ) / 2.0
        wT_bot[i, j] =  h_mean * DIV[i, j] * Td_mean
    end

    # Calculate each term
    for i=1:Nx, j=1:Ny

        isnan(SST[i, j, 1]) && continue

        # Temperature change
        dTdts[i, j, m] += (T_next[t] - T[t]) * (h[t] + h_next[t]) / 2.0 / Δts[m]

        # Sensible/Solar heat flux
        Fs[i, j, m] += (F[i, j, t] + F[i, j, t+1]) / 2.0
    
        # Entrainment
        Δh = h[i, j, t+1] - h[i, j, t]

        if Δh > 0.0
            Ents[i, j, m] += - ρ * c_p * (
                (SST[i, j, t+1] - getTd(-h[i, j, t+1], i, j, TEMP_1, BOT_TEMP_1))
              + (SST[i, j, t  ] - getTd(-h[i, j, t  ], i, j, TEMP_0, BOT_TEMP_0))
            ) * Δh  / 2.0 / Δts[m]

        end

        # Ekman advection
        EKs[i, j, m] += - DIV_UT[i, j] - wT_bot[i, j]

    end

end

for var in [dTds, Fs, Ents, EKs]

    # remember we discarded the first and last year
    var /= (nyears - 2.0)

end
    
Qs = ρ * c_p dTdts - Fs - Ents - EKs


for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
        
        for var in (Qs, dTds, Fs, Ents, EKs)
            # remember we discarded the first and last year
            var[i, j, :] .= NaN
        end
 
    else 

        for var in (Qs, dTds, Fs, Ents, EKs)
            var = view(var, i, j, :)
            var[:] = (var + circshift(var, 1) ) / 2.0
        
            if any(isnan.(var))
                println(format("Position (i, j) = ({:d}, {:d}) has NaN", i, j))
            end
        end


    end

end

# MLD
for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
        hs[i, j, :] .= NaN
        continue
    end 
    
    hs[i, j, :] = mean(reshape(h, 12, :), dims=2)[:, 1]

end

println("Output file...")

Dataset(out_file, "c") do ds

    defDim(ds,"time", 12)
    defDim(ds,"Nx", Nx)
    defDim(ds,"Ny", Ny)

    for o in (
        [
            "h", h, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "qflux", Qs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            )
        ],

    )
        varname, vardata, vardims, varatts = o
        println("Writing ", varname, " with size: ", size(vardata), " ; dim: ", vardims)

        ncvar = defVar(ds, varname, eltype(vardata), vardims)
        ncvar.attrib["_FillValue"] = missing_value
        for key in keys(varatts)
            ncvar.attrib[key] = varatts[key]
        end

        ncvar[:] = vardata
        println("done.")
    end
end

println("Output file: ", out_file)
