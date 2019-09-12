include("../../src/share/DisplacedPoleCoordinate.jl")
include("../../src/share/MapInfo.jl")
include("../../src/share/constants.jl")

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
in_zcoord = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
gridinfo_file = "domain.ocn.gx3v7.120323.nc"

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

    Nt=60

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

function getdTdz(z, i, j, temp)

    zidx = findZInd(z)

    local dTdz
    if zidx == 1
        dTdz = (temp[i, j, zidx] - temp[i, j, zidx+1]) / ((hs[zidx] + hs[zidx+1]) / 2.0)
    else
        dTdz = (temp[i, j, zidx-1] - temp[i, j, zidx+1]) / (0.5*hs[zidx-1] + hs[zidx] + 0.5*hs[zidx+1])
    end
    
    return ( isfinite(dTdz) ) ? dTdz : 0.0
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

ϵs = mi.xc * 0.0 .+ 1e-5
fs = 2 * Ωe * cos.(mi.yc * π/180.0)

K_v = 1e-2

Qs = zeros(Float64, Nx, Ny, 12)


hdTdts = zeros(Float64, Nx, Ny, 12)
Fs    = zeros(Float64, Nx, Ny, 12)
Ents  = zeros(Float64, Nx, Ny, 12)
Qs    = zeros(Float64, Nx, Ny, 12)
hs    = zeros(Float64, Nx, Ny, 12) 
EKs   = zeros(Float64, Nx, Ny, 12) 

τx = zeros(Float64, Nx, Ny)
τy = zeros(Float64, Nx, Ny)

u       = zeros(Float64, Nx, Ny)
v       = zeros(Float64, Nx, Ny)
T       = zeros(Float64, Nx, Ny)
T_hadvs = zeros(Float64, Nx, Ny)
T_vadvs = zeros(Float64, Nx, Ny)
DIV     = zeros(Float64, Nx, Ny)

uT     = zeros(Float64, Nx, Ny)
vT     = zeros(Float64, Nx, Ny)
DIV_UT = zeros(Float64, Nx, Ny)

wT_bot = zeros(Float64, Nx, Ny)


for t = 13:Nt - 12

    print(format("\rProgress: {:d}", t))

    global TEMP, BOT_TEMP

    # month
    m = (t-1)%12+1

    # Loading Temperature profile
    TEMP_0, BOT_TEMP_0 = getTEMPProfile(t  )
    TEMP_1, BOT_TEMP_1 = getTEMPProfile(t+1)
    TEMP_half = (TEMP_0 + TEMP_1) / 2.0

    # Transform input wind stress vector first
    
    DisplacedPoleCoordinate.project!(
        gi,
        (view(TAUX, :, :, t) + view(TAUX, :, :, t+1))/2.0,
        (view(TAUY, :, :, t) + view(TAUY, :, :, t+1))/2.0,
        τx, τy,
        direction=:Forward,
    )

#    τx = (view(TAUX, :, :, t) + view(TAUX, :, :, t+1))/2.0
#    τy = (view(TAUY, :, :, t) + view(TAUY, :, :, t+1))/2.0


    # Calculate Ekman flow
    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue

        h = (HMXL[i, j, t] + HMXL[i, j, t+1] ) / 2.0
        s̃ = ϵs[i, j] + fs[i, j] * im
        H̃ = √(K_v / s̃)
        H = abs(H̃)
        H_ek = max(h, 2*H)
        
        M̃ = (τx[i, j] + τy[i, j] * im) / (ρ * s̃)
        ṽ_ek =   M̃ / H_ek
        u_ek, v_ek = real(ṽ_ek), imag(ṽ_ek)



        u[i, j] = u_ek
        v[i, j] = v_ek
        T[i, j] = (SST[i, j, t] + SST[i, j, t+1])/2.0

        #uT[i, j] = u_ek * SST_mean
        #vT[i, j] = v_ek * SST_mean

    end
   
    DisplacedPoleCoordinate.hadv_upwind!(
        gi,
        T_hadvs,
        u,
        v,
        T,
        mi.mask,
    )

    DisplacedPoleCoordinate.DIV!(gi, u,  v,   DIV,    mi.mask)
    #DisplacedPoleCoordinate.DIV!(gi, uT, vT,  DIV_UT, mi.mask)
    #=
    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue
        
        h_mean   = (HMXL[i, j, t] + HMXL[i, j, t+1] ) / 2.0
        Td_mean = ( getTd(-HMXL[i, j, t+1], i, j, TEMP_1, BOT_TEMP_1) +  getTd(-HMXL[i, j, t], i, j, TEMP_0, BOT_TEMP_0) ) / 2.0
        wT_bot[i, j] =  h_mean * DIV[i, j] * Td_mean

        if (i,j) == (70,55)
            println("h_mean: ", h_mean, "; ", DIV[i, j], "; Td:", Td_mean)
            println("DIV UT: ", DIV_UT[i, j])
        end
    end
    =#

    # Calculate each term
    for i=1:Nx, j=1:Ny

        isnan(SST[i, j, 1]) && continue

        # Temperature change
        hdTdts[i, j, m] += (SST[i, j, t+1] - SST[i, j, t]) * (HMXL[i, j, t] + HMXL[i, j, t+1]) / 2.0 / Δts[m]

        # Sensible/Solar heat flux
        Fs[i, j, m] += (SHF[i, j, t] + SHF[i, j, t+1]) / 2.0
    
        # Entrainment
        Δh = HMXL[i, j, t+1] - HMXL[i, j, t]

        if Δh > 0.0
            Ents[i, j, m] += - ρ * c_p * (
                (SST[i, j, t+1] - getTd(-HMXL[i, j, t+1], i, j, TEMP_1, BOT_TEMP_1))
              + (SST[i, j, t  ] - getTd(-HMXL[i, j, t  ], i, j, TEMP_0, BOT_TEMP_0))
            ) * Δh  / 2.0 / Δts[m]

        end

        # Ekman advection
        #EKs[i, j, m] += - ρ * c_p * ( DIV_UT[i, j] )# + wT_bot[i, j] )
        h_mean = (HMXL[i, j, t] + HMXL[i, j, t+1]) / 2.0
        T_vadvs[i, j] = - DIV[i, j] * h_mean * getdTdz(-h_mean, i, j, TEMP_half)
        EKs[i, j, m] += ρ * c_p * (T_hadvs[i, j] + T_vadvs[i, j])

    end

end
    
#for i=1:Nx, j=1:Ny, m=1:12
#    Fs[i, j, m] /= (nyears - 2.0)
#end

for var in [hdTdts, Ents, EKs, Fs]

    for i=1:Nx, j=1:Ny, m=1:12
        # remember we discarded the first and last year
        var[i, j, m] /= (nyears - 2.0)
    end

end
    
Qs = ρ * c_p * hdTdts - Fs - Ents - EKs


for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
        
        for var in (Qs, hdTdts, Fs, Ents, EKs)
            # remember we discarded the first and last year
            var[i, j, :] .= NaN
        end
 
    else 

        for var in (Qs, hdTdts, Fs, Ents, EKs)
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
    
    hs[i, j, :] = mean(reshape(view(HMXL, i, j, :), 12, :), dims=2)[:, 1]

end

println("Output file...")

Dataset(out_file, "c") do ds

    defDim(ds,"time", 12)
    defDim(ds,"Nx", Nx)
    defDim(ds,"Ny", Ny)

    for o in (
        [
            "h", hs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "qflux", Qs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            )
        ], [
            "Fs", Fs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"SHF flux",
            "units"=>"W / m^2",
            )
        ], [
            "Ents", Ents, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Entrainment flux",
            "units"=>"W / m^2",
            )

        ], [
            "EKs", EKs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"W / m^2",
            )

        ], [
            "hdTdts", hdTdts * ρ * c_p, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
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
