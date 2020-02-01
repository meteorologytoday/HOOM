include("../../src/share/DisplacedPoleCoordinate.jl")
include("../../src/share/MapInfo.jl")
include("../../src/share/constants.jl")

using Formatting
using NCDatasets
using Statistics
using SharedArrays

using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--years"
            help = "Processed years."
            arg_type = Int64
            default = -1
     
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))


ent_thickness = 10.0
ρ    = 1026.0  # kg / m^3
c_p  = 3996.0  # J / kg / K

dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
sum(dom) == 365 || throw(ErrorException("Sum of dom is not 365"))

Δts = (dom + circshift(dom, -1))/2.0 * 86400.0

println("Δts = ", Δts)

in_SSH  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SSH.100001-109912.nc"
in_SST  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SST.100001-109912.nc"
in_SHF  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc" 
in_TAUX = "b.e11.B1850C5CN.f09_g16.005.pop.h.TAUX.100001-109912.nc"
in_TAUY = "b.e11.B1850C5CN.f09_g16.005.pop.h.TAUY.100001-109912.nc"
in_HMXL = "b.e11.B1850C5CN.f09_g16.005.pop.h.HMXL.100001-109912.nc"
in_HBLT = "b.e11.B1850C5CN.f09_g16.005.pop.h.HBLT.100001-109912.nc"
in_TEMP = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
in_QFLUX = "b.e11.B1850C5CN.f09_g16.005.pop.h.QFLUX.100001-109912.nc"
in_MELTH_F = "b.e11.B1850C5CN.f09_g16.005.pop.h.MELTH_F.100001-109912.nc"





in_zcoord = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
gridinfo_file = "domain.ocn.gx3v7.120323.nc"
gridinfo_file = "domain.ocn.gx1v6.090206.nc"

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
    global SST = replace(ds["SST"][:, :, :], missing=>NaN)
    global Nx, Ny, Nt = size(SST)

    (Nt%12 == 0) || throw(ErrorException("Time is not multiple of 12"))

    if parsed["years"] > 0
        Nt = parsed["years"] * 12
    end

    global nyears = Int(Nt / 12)

    println("YEARS: ", nyears)
end

Dataset(in_SHF, "r") do ds
    global SHF = replace(ds["SHF"][:], missing=>NaN)
end

Dataset(in_HMXL, "r") do ds
    global HMXL = replace(ds["HMXL"][:], missing=>NaN) / 100.0
end

#Dataset(in_HBLT, "r") do ds
#    global HBLT = replace(ds["HBLT"][:], missing=>NaN) / 100.0
#end


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

Dataset(in_SSH, "r") do ds
    global SSH = replace(ds["SSH"][:], missing=>NaN) / 100.0
end


SURF = SHF + QFLUX

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
gi = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv, mi.area; angle_unit=:deg)


ϵs = mi.xc * 0.0 .+ 1e-5
fs = 2 * Ωe * sin.(mi.yc * π/180.0)

K_v = 1e-2


Q1s             = SharedArray{Float64}(Nx, Ny, 12)
var_Q1s         = SharedArray{Float64}(Nx, Ny, 12)
Q1s_energy      = SharedArray{Float64}(Nx, Ny, 12)
var_Q1s_energy  = SharedArray{Float64}(Nx, Ny, 12)

Q2s                 = SharedArray{Float64}(Nx, Ny, 12)
var_Q2s             = SharedArray{Float64}(Nx, Ny, 12)
Q2s_energy          = SharedArray{Float64}(Nx, Ny, 12)
var_Q2s_energy      = SharedArray{Float64}(Nx, Ny, 12)

Q3s                 = SharedArray{Float64}(Nx, Ny, 12)
var_Q3s             = SharedArray{Float64}(Nx, Ny, 12)
Q3s_energy          = SharedArray{Float64}(Nx, Ny, 12)
var_Q3s_energy      = SharedArray{Float64}(Nx, Ny, 12)

Q4s                 = SharedArray{Float64}(Nx, Ny, 12)
var_Q4s             = SharedArray{Float64}(Nx, Ny, 12)
Q4s_energy          = SharedArray{Float64}(Nx, Ny, 12)
var_Q4s_energy      = SharedArray{Float64}(Nx, Ny, 12)

Q5s                 = SharedArray{Float64}(Nx, Ny, 12)
var_Q5s             = SharedArray{Float64}(Nx, Ny, 12)
Q5s_energy          = SharedArray{Float64}(Nx, Ny, 12)
var_Q5s_energy      = SharedArray{Float64}(Nx, Ny, 12)






dTdts    = SharedArray{Float64}(Nx, Ny, 12)
Fs       = SharedArray{Float64}(Nx, Ny, 12)
Ents     = SharedArray{Float64}(Nx, Ny, 12)
hs       = SharedArray{Float64}(Nx, Ny, 12) 
sea_ices = SharedArray{Float64}(Nx, Ny, 12) 
diffs    = SharedArray{Float64}(Nx, Ny, 12) 

T_hadvs = SharedArray{Float64}(Nx, Ny, 12)
T_vadvs = SharedArray{Float64}(Nx, Ny, 12)
DIVs = SharedArray{Float64}(Nx, Ny, 12)
us = SharedArray{Float64}(Nx, Ny, 12)
vs = SharedArray{Float64}(Nx, Ny, 12)

ugs = SharedArray{Float64}(Nx, Ny, 12)
vgs = SharedArray{Float64}(Nx, Ny, 12)
T_gadvs = SharedArray{Float64}(Nx, Ny, 12)


var_dTdts = SharedArray{Float64}(Nx, Ny, 12)
var_Fs    = SharedArray{Float64}(Nx, Ny, 12)
var_Ents  = SharedArray{Float64}(Nx, Ny, 12)
var_hs    = SharedArray{Float64}(Nx, Ny, 12) 
var_sea_ices = SharedArray{Float64}(Nx, Ny, 12) 
var_diff    = SharedArray{Float64}(Nx, Ny, 12) 

var_T_hadvs = SharedArray{Float64}(Nx, Ny, 12)
var_T_vadvs = SharedArray{Float64}(Nx, Ny, 12)
var_T_gadvs = SharedArray{Float64}(Nx, Ny, 12)





τx = SharedArray{Float64}(Nx, Ny)
τy = SharedArray{Float64}(Nx, Ny)

u           = SharedArray{Float64}(Nx, Ny)
v           = SharedArray{Float64}(Nx, Ny)
T           = SharedArray{Float64}(Nx, Ny)
tmp_DIV     = SharedArray{Float64}(Nx, Ny)
tmp_T_hadvs = SharedArray{Float64}(Nx, Ny)

ug          = SharedArray{Float64}(Nx, Ny)
vg          = SharedArray{Float64}(Nx, Ny)
tmp_T_gadvs = SharedArray{Float64}(Nx, Ny)

wT_bot = SharedArray{Float64}(Nx, Ny)

tmp_u_dis   = SharedArray{Float64}(Nx, Ny)
tmp_v_dis   = SharedArray{Float64}(Nx, Ny)

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

    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue
        T[i, j] = (SST[i, j, t] + SST[i, j, t+1])/2.0
    end
 
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

        tmp_u_dis[i, j] = u_ek
        tmp_v_dis[i, j] = v_ek
    end
   
    DisplacedPoleCoordinate.hadv_upwind!(
        gi,
        tmp_T_hadvs,
        tmp_u_dis,
        tmp_v_dis,
        T,
        mi.mask,
    )
 
    DisplacedPoleCoordinate.DIV!(
        gi,
        tmp_u_dis,
        tmp_v_dis,
        tmp_DIV,
        mi.mask,
    )

    DisplacedPoleCoordinate.project!(
        gi,
        tmp_u_dis,
        tmp_v_dis,
        u,
        v,
        direction = :Backward,
    )

   
    # Calculate geostrophic flow
    DisplacedPoleCoordinate.GRAD!(
        gi,
        (SSH[:, :, t] + SSH[:, :, t+1]) / 2.0,
        tmp_u_dis,
        tmp_v_dis,
        mi.mask,
    )

    for i=1:Nx, j=1:Ny
        isnan(SST[i, j, 1]) && continue

        s̃ = ϵs[i, j] + fs[i, j] * im
        ṽ = - g * (tmp_u_dis[i, j] + tmp_v_dis[i, j] * im) / s̃
        tmp_u_dis[i, j], tmp_v_dis[i, j] = real(ṽ), imag(ṽ)
    end
 
    DisplacedPoleCoordinate.hadv_upwind!(
        gi,
        tmp_T_gadvs,
        tmp_u_dis,
        tmp_v_dis,
        T,
        mi.mask,
    )
 
    DisplacedPoleCoordinate.project!(
        gi,
        tmp_u_dis,
        tmp_v_dis,
        ug,
        vg,
        direction = :Backward,
    )



    # Calculate each term
    #@distributed for idx in CartesianIndices((1:Nx, 1:Ny))
#        i = idx[1]
#        j = idx[2]


    for i=1:Nx, j=1:Ny

        isnan(SST[i, j, 1]) && continue
        
        # Using the time-variate of HMXL to estimate h_mean
        #  is problematic because lacking of time resolution
        # making entrainment a source of energy to the ocean
        h_mean = (HMXL[i, j, t] + HMXL[i, j, t+1]) / 2.0

        #h_mean = mean(HMXL[i, j, :]) 

        # Temperature change
        tmp_dTdts = (SST[i, j, t+1] - SST[i, j, t]) / Δts[m]
        var_dTdts[i, j, m] += tmp_dTdts^2.0
        dTdts[i, j, m] += tmp_dTdts

        # Sea-ice melting/formation
        tmp_sea_ices = ( SEAICE_HFLX[i, j, t] + SEAICE_HFLX[i, j, t+1] ) / 2.0 / h_mean / ρ / c_p
        sea_ices[i, j, m] += tmp_sea_ices
        var_sea_ices[i, j, m] += tmp_sea_ices^2.0

        # Surface downward heat flux 
        tmp_Fs = (SURF[i, j, t] + SURF[i, j, t+1]) / 2.0 / h_mean / ρ / c_p
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
#            tmp_Ents = - (
#                (SST[i, j, t] - getTd(-HMXL[i, j, t+1], i, j, TEMP_0, BOT_TEMP_0))
#            ) * Δh  / Δts[m] / h_mean

            tmp_Ents = - (
                integrate(-HMXL[i, j, t+1], -HMXL[i, j, t], zs, SST[i, j, t] .- TEMP_0[i, j, :])
            )  / Δts[m] / h_mean


            Ents[i, j, m] += tmp_Ents
            var_Ents[i, j, m] += tmp_Ents^2.0


        end

        # Ekman advection
        tmp_T_vadvs = - tmp_DIV[i, j] * (
                (SST[i, j, t+1] - getTd(-HMXL[i, j, t+1]-ent_thickness, i, j, TEMP_1, BOT_TEMP_1))
              + (SST[i, j, t  ] - getTd(-HMXL[i, j, t  ]-ent_thickness, i, j, TEMP_0, BOT_TEMP_0))
        ) / 2.0

        T_hadvs[i, j, m] += tmp_T_hadvs[i ,j]
        T_vadvs[i, j, m] += tmp_T_vadvs

        var_T_hadvs[i, j, m] += tmp_T_hadvs[i, j]^2.0
        var_T_vadvs[i, j, m] += tmp_T_vadvs^2.0

        DIVs[i, j, m] += tmp_DIV[i, j]
        us[i, j, m] += u[i, j]
        vs[i, j, m] += v[i, j]

        # Geostrophic advection
        T_gadvs[i, j, m] += tmp_T_gadvs[i ,j]
        var_T_gadvs[i, j, m] += tmp_T_gadvs[i, j]^2.0

        ugs[i, j, m] += ug[i, j]
        vgs[i, j, m] += vg[i, j]


        #EKs[i, j, m] += ρ * c_p * (tmp_T_hadvs[i, j] + tmp_T_vadvs) * h_mean / h_mean

        tmp_Q4s = - ( tmp_dTdts - ( tmp_Fs + tmp_Ents + (tmp_T_hadvs[i, j] + tmp_T_vadvs)) )
        Q4s[i, j, m]     += tmp_Q4s
        var_Q4s[i, j, m] += tmp_Q4s^2.0

        tmp_Q4s_energy = tmp_Q4s * h_mean * ρ * c_p
        Q4s_energy[i, j, m] += tmp_Q4s_energy
        var_Q4s_energy[i, j, m] += tmp_Q4s_energy^2.0


        tmp_Q5s = - ( tmp_dTdts - ( tmp_Fs + tmp_Ents + (tmp_T_hadvs[i, j] + tmp_T_vadvs + tmp_T_gadvs[i, j])) )
        Q5s[i, j, m]     += tmp_Q5s
        var_Q5s[i, j, m] += tmp_Q5s^2.0

        tmp_Q5s_energy = tmp_Q5s * h_mean * ρ * c_p
        Q5s_energy[i, j, m] += tmp_Q5s_energy
        var_Q5s_energy[i, j, m] += tmp_Q5s_energy^2.0


        # Q2
        tmp_Q2s = - ( tmp_dTdts - ( tmp_Fs + tmp_Ents ) )
        Q2s[i, j, m]     += tmp_Q2s
        var_Q2s[i, j, m] += tmp_Q2s^2.0

        tmp_Q2s_energy = tmp_Q2s * h_mean * ρ * c_p
        Q2s_energy[i, j, m] += tmp_Q2s_energy
        var_Q2s_energy[i, j, m] += tmp_Q2s_energy^2.0




        # calculate Q1, Q3 qflux
        h_mean = mean(HMXL[i, j, :]) 
        tmp_dTdts = (SST[i, j, t+1] - SST[i, j, t]) / Δts[m]
        tmp_Fs    = (SURF[i, j, t] + SURF[i, j, t+1]) / 2.0 / h_mean / ρ / c_p 
        
        tmp_Q1s = - ( tmp_dTdts - tmp_Fs )
        Q1s[i, j, m]     += tmp_Q1s
        var_Q1s[i, j, m] += tmp_Q1s^2.0


        tmp_Q1s_energy = tmp_Q1s * h_mean * ρ * c_p
        Q1s_energy[i, j, m] += tmp_Q1s_energy
        var_Q1s_energy[i, j, m] += tmp_Q1s_energy^2.0


        tmp_Q3s = - ( tmp_dTdts - ( tmp_Fs + (tmp_T_hadvs[i, j] + tmp_T_vadvs)) )
        Q3s[i, j, m]     += tmp_Q3s
        var_Q3s[i, j, m] += tmp_Q3s^2.0

        tmp_Q3s_energy = tmp_Q3s * h_mean * ρ * c_p
        Q3s_energy[i, j, m] += tmp_Q3s_energy
        var_Q3s_energy[i, j, m] += tmp_Q3s_energy^2.0


    end

end

needed_vars = (
    dTdts,
    Ents,
    Fs,
    sea_ices,
    diffs,
    T_hadvs,
    T_vadvs,
    T_gadvs,
    var_dTdts,
    var_Ents,
    var_Fs,
    var_sea_ices,
    var_T_hadvs,
    var_T_vadvs,
    var_T_gadvs,
    DIVs,
    us,
    vs,
    Q1s,
    var_Q1s,
    Q1s_energy,
    var_Q1s_energy,
    Q2s,
    var_Q2s,
    Q2s_energy,
    var_Q2s_energy,
    Q3s,
    var_Q3s,
    Q3s_energy,
    var_Q3s_energy,
    Q4s,
    var_Q4s,
    Q4s_energy,
    var_Q4s_energy,
    Q5s,
    var_Q5s,
    Q5s_energy,
    var_Q5s_energy,

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

var_Q1s        -= Q1s.^2.0
var_Q2s        -= Q2s.^2.0
var_Q3s        -= Q3s.^2.0
var_Q4s        -= Q4s.^2.0
var_Q5s        -= Q5s.^2.0

var_Q1s_energy -= Q1s_energy.^2.0
var_Q2s_energy -= Q2s_energy.^2.0
var_Q3s_energy -= Q3s_energy.^2.0
var_Q4s_energy -= Q4s_energy.^2.0
var_Q5s_energy -= Q5s_energy.^2.0



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
            "h_ML_variate", hs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "h_ML_fixed", repeat(mean(hs, dims=(3,))[:,:,1], outer=(1,1,12)), ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "Q1s", Q1s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q1s", var_Q1s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q1s_energy", Q1s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q1s_energy", var_Q1s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q2s", Q2s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q2s", var_Q2s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q2s_energy", Q2s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q2s_energy", var_Q2s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q3s", Q3s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q3s", var_Q3s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q3s_energy", Q3s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q3s_energy", var_Q3s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q4s", Q4s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "var_Q4s", var_Q4s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "Q4s_energy", Q4s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q4s_energy", var_Q4s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )
        ], [
            "Q5s", Q5s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "var_Q5s", var_Q5s, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"K / s",
            )
        ], [
            "Q5s_energy", Q5s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
            "units"=>"W / m^2",
            )
        ], [
            "var_Q5s_energy", var_Q5s_energy, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Variance of qflux_energy",
            "units"=>"(W / m^2)^2",
            )

        ], [
            "Fs", Fs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"SURF flux",
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
            "T_gadvs", T_gadvs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"K / s",
            )
        ], [
            "var_Fs", var_Fs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"SURF flux",
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
            "var_T_gadvs", var_T_gadvs, ("Nx", "Ny", "time"), Dict(
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
            "ugs", ugs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Ekman flux",
            "units"=>"m / s",
            )
        ], [
            "vgs", vgs, ("Nx", "Ny", "time"), Dict(
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
