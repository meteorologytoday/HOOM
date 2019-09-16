include("../../src/share/DisplacedPoleCoordinate.jl")
include("../../src/share/MapInfo.jl")
include("../../src/share/constants.jl")

using Formatting
using NCDatasets
using Statistics
using SharedArrays

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




in_SST  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SST.100001-109912.nc"
in_SHF  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc" 
in_TAUX = "b.e11.B1850C5CN.f09_g16.005.pop.h.TAUX.100001-109912.nc"
in_TAUY = "b.e11.B1850C5CN.f09_g16.005.pop.h.TAUY.100001-109912.nc"
in_HMXL = "b.e11.B1850C5CN.f09_g16.005.pop.h.HMXL.100001-109912.nc"
in_TEMP = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
in_QFLUX = "b.e11.B1850C5CN.f09_g16.005.pop.h.QFLUX.100001-109912.nc"
in_MELTH_F = "b.e11.B1850C5CN.f09_g16.005.pop.h.MELTH_F.100001-109912.nc"





in_zcoord = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
gridinfo_file = "domain.ocn.gx3v7.120323.nc"
gridinfo_file = "domain.ocn.gx1v6.090206.nc"

out_file = "forcing.gx3v7.nc"
out_file = "forcing.gx1v6.nc"

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
    #Nt=36

    global nyears = Int(Nt / 12)

end

Dataset(in_SHF, "r") do ds
    global SHF = replace(ds["SHF"][:], missing=>NaN)
end


Dataset(in_HMXL, "r") do ds
    global HMXL = replace(ds["HMXL"][:], missing=>NaN) / 100.0
end

Dataset(in_TAUX, "r") do ds
    global TAUX = replace(ds["TAUX"][:], missing=>NaN) / 10.0
end

Dataset(in_TAUY, "r") do ds
    global TAUY = replace(ds["TAUY"][:], missing=>NaN) / 10.0
end

Dataset(in_QFLUX, "r") do ds
    global QFLUX = replace(ds["QFLUX"][:], missing=>NaN)
end

Dataset(in_MELTH_F, "r") do ds
    global MELTH_F = replace(ds["MELTH_F"][:], missing=>NaN)
end


SEAICE_HFLX = QFLUX + MELTH_F



Dataset(in_zcoord, "r") do ds
    z_top = - replace(ds["z_w_top"][:], missing=>NaN) / 100.0
    z_bot = - replace(ds["z_w_bot"][:], missing=>NaN) / 100.0

    println("z_top: ", z_top)
    println("z_bot: ", z_bot)
    
    global zs = [z_top..., z_bot[end]]
    global Δhs = z_top - z_bot
    global Nz = length(zs) - 1

end


println("zs: ", zs)
println("Δhs: ", Δhs)

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
        dTdz = (temp[i, j, zidx] - temp[i, j, zidx+1]) / ((Δhs[zidx] + Δhs[zidx+1]) / 2.0)
    elseif isnan(temp[i, j, zidx+1])
        dTdz = (temp[i, j, zidx-1] - temp[i, j, zidx]) / (0.5*Δhs[zidx-1]  + 0.5*Δhs[zidx])
    else
        dTdz = (temp[i, j, zidx-1] - temp[i, j, zidx+1]) / (0.5*Δhs[zidx-1] + Δhs[zidx] + 0.5*Δhs[zidx+1])
    end
    if isnan(dTdz)

        println("z: ", z)
        println("zidx: ", zidx)
        println("Δhs[zidx+1]: ", Δhs[zidx+1])
        println("Δhs[zidx-1]: ", Δhs[zidx-1])
        println("temp: ", temp[i, j, :])
        throw(ErrorException("dTdz nan"))
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
            s += Δhs[k] * Ts[k]
        end

        s += (zs[kb] - zb) * Ts[kb]

    end

end

mi = ModelMap.MapInfo{Float64}(gridinfo_file)
gi = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv; angle_unit=:deg)


ent_thickness = 10.0

ϵs = mi.xc * 0.0 .+ 1e-5
fs = 2 * Ωe * sin.(mi.yc * π/180.0)

K_v = 1e-2

Qs_energy      = SharedArray{Float64}(Nx, Ny, 12)
var_Qs_energy  = SharedArray{Float64}(Nx, Ny, 12)

Qs_energy_OLD      = SharedArray{Float64}(Nx, Ny, 12)
var_Qs_energy_OLD  = SharedArray{Float64}(Nx, Ny, 12)

Qs_OLD      = SharedArray{Float64}(Nx, Ny, 12)
var_Qs_OLD  = SharedArray{Float64}(Nx, Ny, 12)



dTdts = SharedArray{Float64}(Nx, Ny, 12)
Fs    = SharedArray{Float64}(Nx, Ny, 12)
Ents  = SharedArray{Float64}(Nx, Ny, 12)
Qs    = SharedArray{Float64}(Nx, Ny, 12)
hs    = SharedArray{Float64}(Nx, Ny, 12) 
#EKs   = SharedArray{Float64}(Nx, Ny, 12) 
sea_ices = SharedArray{Float64}(Nx, Ny, 12) 

T_hadvs = SharedArray{Float64}(Nx, Ny, 12)
T_vadvs = SharedArray{Float64}(Nx, Ny, 12)
DIVs = SharedArray{Float64}(Nx, Ny, 12)
us = SharedArray{Float64}(Nx, Ny, 12)
vs = SharedArray{Float64}(Nx, Ny, 12)

var_dTdts = SharedArray{Float64}(Nx, Ny, 12)
var_Fs    = SharedArray{Float64}(Nx, Ny, 12)
var_Ents  = SharedArray{Float64}(Nx, Ny, 12)
var_Qs    = SharedArray{Float64}(Nx, Ny, 12)
var_hs    = SharedArray{Float64}(Nx, Ny, 12) 
#var_EKs   = SharedArray{Float64}(Nx, Ny, 12) 
var_sea_ices = SharedArray{Float64}(Nx, Ny, 12) 

var_T_hadvs = SharedArray{Float64}(Nx, Ny, 12)
var_T_vadvs = SharedArray{Float64}(Nx, Ny, 12)




τx = SharedArray{Float64}(Nx, Ny)
τy = SharedArray{Float64}(Nx, Ny)

u           = SharedArray{Float64}(Nx, Ny)
v           = SharedArray{Float64}(Nx, Ny)
T           = SharedArray{Float64}(Nx, Ny)
tmp_DIV     = SharedArray{Float64}(Nx, Ny)
tmp_T_hadvs = SharedArray{Float64}(Nx, Ny)

wT_bot = SharedArray{Float64}(Nx, Ny)

# MLD
for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
        hs[i, j, :] .= NaN
        continue
    end 
    
    hs[i, j, :] = mean(reshape(view(HMXL, i, j, :), 12, :), dims=2)[:, 1]

end


for t = 13:Nt - 12

    print(format("\rCalculate New Qflux Progress: {:d}", t))

    global TEMP, BOT_TEMP

    # month
    m = (t-1)%12+1

    # Loading Temperature profile
    TEMP_0, BOT_TEMP_0 = getTEMPProfile(t  )
    TEMP_1, BOT_TEMP_1 = getTEMPProfile(t+1)
    TEMP_half = (TEMP_0 + TEMP_1) / 2.0

    # Transform input wind stress vector first

    #=    
    DisplacedPoleCoordinate.project!(
        gi,
        (view(TAUX, :, :, t) + view(TAUX, :, :, t+1))/2.0,
        (view(TAUY, :, :, t) + view(TAUY, :, :, t+1))/2.0,
        τx, τy,
        direction=:Forward,
    )
    =#
    τx = (view(TAUX, :, :, t) + view(TAUX, :, :, t+1))/2.0
    τy = (view(TAUY, :, :, t) + view(TAUY, :, :, t+1))/2.0


    # Calculate Ekman flow
    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue

        h = (HMXL[i, j, t] + HMXL[i, j, t+1] ) / 2.0
        s̃ = ϵs[i, j] + fs[i, j] * im
        H̃ = √(K_v / s̃)
        H = abs(H̃)
        H_ek = max(h, 2*H)
        
        M̃ = (τx[i, j] + τy[i, j] * im) / s̃
        ṽ_ek =   M̃ / (ρ * H_ek)
        u_ek, v_ek = real(ṽ_ek), imag(ṽ_ek)

        u[i, j] = u_ek
        v[i, j] = v_ek
        T[i, j] = (SST[i, j, t] + SST[i, j, t+1])/2.0

    end
   
    DisplacedPoleCoordinate.hadv_upwind!(
        gi,
        tmp_T_hadvs,
        u,
        v,
        T,
        mi.mask,
    )

    DisplacedPoleCoordinate.DIV!(gi, u,  v, tmp_DIV, mi.mask)

    # Calculate each term
    #@distributed for idx in CartesianIndices((1:Nx, 1:Ny))
        i = idx[1]
        j = idx[2]


    for i=1:Nx, j=1:Ny

        isnan(SST[i, j, 1]) && continue
        
        h_mean = (HMXL[i, j, t] + HMXL[i, j, t+1]) / 2.0

        # Temperature change
        tmp_dTdts = (SST[i, j, t+1] - SST[i, j, t]) / Δts[m]
        var_dTdts[i, j, m] += tmp_dTdts^2.0
        dTdts[i, j, m] += tmp_dTdts

        # Sea-ice melting/formation
        tmp_sea_ices = ( SEAICE_HFLX[i, j, t] + SEAICE_HFLX[i, j, t+1] ) / 2.0 / h_mean / ρ / c_p
        sea_ices[i, j, m] += tmp_sea_ices
        var_sea_ices[i, j, m] += tmp_sea_ices^2.0

        # Surface downward heat flux heat flux
        tmp_Fs = (SHF[i, j, t] + SHF[i, j, t+1]) / 2.0 / h_mean / ρ / c_p
        Fs[i, j, m] += tmp_Fs
        var_Fs[i, j, m] += tmp_Fs^2.0

        # Entrainment
        Δh = HMXL[i, j, t+1] - HMXL[i, j, t]

        tmp_Ents = 0.0
        if Δh > 0.0
#=
            Ents[i, j, m] += - ρ * c_p * (
                (SST[i, j, t+1] - getTd(-HMXL[i, j, t+1]-ent_thickness, i, j, TEMP_1, BOT_TEMP_1))
              + (SST[i, j, t  ] - getTd(-HMXL[i, j, t  ]-ent_thickness, i, j, TEMP_0, BOT_TEMP_0))
            ) * Δh  / 2.0 / Δts[m] / h_mean
=#
            tmp_Ents = - (
                (SST[i, j, t] - getTd(-HMXL[i, j, t+1], i, j, TEMP_0, BOT_TEMP_0))
            ) * Δh  / Δts[m] / h_mean

            Ents[i, j, m] += tmp_Ents
            var_Ents[i, j, m] += tmp_Ents^2.0


        end

        # Ekman advection
        tmp_T_vadvs = - tmp_DIV[i, j] * getdTdz(-h_mean, i, j, TEMP_half) * 10.0

        T_hadvs[i, j, m] += tmp_T_hadvs[i ,j]
        T_vadvs[i, j, m] += tmp_T_vadvs

        var_T_hadvs[i, j, m] += tmp_T_hadvs[i, j]^2.0
        var_T_vadvs[i, j, m] += tmp_T_vadvs^2.0

        DIVs[i, j, m] += tmp_DIV[i, j]
        us[i, j, m] += u[i, j]
        vs[i, j, m] += v[i, j]

        #EKs[i, j, m] += ρ * c_p * (tmp_T_hadvs[i, j] + tmp_T_vadvs) * h_mean / h_mean


        tmp_Qs = tmp_dTdts - ( tmp_Fs + tmp_Ents + (tmp_T_hadvs[i, j] + tmp_T_vadvs))
        Qs[i, j, m]     += tmp_Qs
        var_Qs[i, j, m] += tmp_Qs^2.0

        tmp_Qs_energy = tmp_Qs * h_mean * ρ * c_p
        Qs_energy[i, j, m] += tmp_Qs_energy
        var_Qs_energy[i, j, m] += tmp_Qs_energy^2.0




        # calculate OLD qflux
        h_mean = mean(HMXL[i, j, :]) 
        
        tmp_Qs_OLD = (SST[i, j, t+1] - SST[i, j, t]) / Δts[m] - (SHF[i, j, t] + SHF[i, j, t+1]) / 2.0 / h_mean / ρ / c_p
        Qs_OLD[i, j, m]     += tmp_Qs_OLD
        var_Qs_OLD[i, j, m] += tmp_Qs_OLD^2.0


        tmp_Qs_energy_OLD = tmp_Qs_OLD * h_mean * ρ * c_p
        Qs_energy_OLD[i, j, m] += tmp_Qs_energy_OLD
        var_Qs_energy_OLD[i, j, m] += tmp_Qs_energy_OLD^2.0


    end

end

needed_vars = (
    dTdts,
    Ents,
    Fs,
    sea_ices,
    T_hadvs,
    T_vadvs,
    var_dTdts,
    var_Ents,
    var_Fs,
    var_sea_ices,
    var_T_hadvs,
    var_T_vadvs,
    DIVs,
    us,
    vs,
    Qs,
    var_Qs,
    Qs_energy,
    var_Qs_energy,
    Qs_energy_OLD,
    var_Qs_energy_OLD,
    Qs_OLD,
    var_Qs_OLD,
)
 
for var in needed_vars
    for i=1:Nx, j=1:Ny, m=1:12
        # remember we discarded the first and last year
        var[i, j, m] /= (nyears - 2.0)
    end
end


var_dTdts   -= dTdts.^2.0
var_Fs      -= Fs.^2.0
var_sea_ices-= sea_ices.^2.0
var_Ents    -= Ents.^2.0
var_T_hadvs -= T_hadvs.^2.0
var_T_vadvs -= T_vadvs.^2.0
var_Qs      -= Qs.^2.0
var_Qs_OLD  -= Qs_OLD.^2.0
var_Qs_energy -= Qs_energy.^2.0
var_Qs_energy_OLD -= Qs_energy_OLD.^2.0



for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
         
        for var in needed_vars
            var[i, j, :] .= NaN
        end
 
    else 

    #=
        for var in needed_vars
            var = view(var, i, j, :)
            var[:] = (var + circshift(var, 1) ) / 2.0
        
            if any(isnan.(var))
                println(format("Position (i, j) = ({:d}, {:d}) has NaN", i, j))
            end
        end
    =#

    end

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
            "qflux_energy_OLD", Qs_energy_OLD, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_qflux_energy_OLD", var_Qs_energy_OLD, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "qflux_energy", Qs_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_qflux_energy", var_Qs_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "qflux_OLD", Qs_OLD, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_qflux_OLD", var_Qs_OLD, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "qflux", Qs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "Fs", Fs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"SHF flux",
            "units"=>"K / s",
            )
        ], [
            "Ents", Ents, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Entrainment flux",
            "units"=>"K / s",
            )

        ], [
            "T_hadvs", T_hadvs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "T_vadvs", T_vadvs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )

        ], [
            "var_qflux", var_Qs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "var_Fs", var_Fs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"SHF flux",
            "units"=>"K / s",
            )
        ], [
            "var_Ents", var_Ents, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Entrainment flux",
            "units"=>"K / s",
            )

        ], [
            "var_T_hadvs", var_T_hadvs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "var_T_vadvs", var_T_vadvs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "sea_ices", sea_ices, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "var_sea_ices", var_sea_ices, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )

        ], [
            "DIVs", DIVs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"1 / s",
            )
        ], [
            "us", us, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"m / s",
            )
        ], [
            "vs", vs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"m / s",
            )

        ], [
            "dTdts", dTdts, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "var_dTdts", var_dTdts, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
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
