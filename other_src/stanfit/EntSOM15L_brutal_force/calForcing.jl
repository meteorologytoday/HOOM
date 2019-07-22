using Formatting
using NCDatasets
using Statistics
ρ    = 1026.0  # kg / m^3
c_p  = 3996.0  # J / kg / K

dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
sum(dom) == 365 || throw(ErrorException("Sum of dom is not 365"))

Δts = (dom + circshift(dom, -1))/2.0 * 86400.0

println("Δts = ", Δts)

in_SST="transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SST.100001-109912.nc"
in_SHF="transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SHF.100001-109912.nc"
in_HMXL="transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.HMXL.100001-109912.nc"
in_TEMP="transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.TEMP.100001-109912.nc"

out_file = "forcing.gx3v7.nc"

missing_value = 1e20

Dataset("transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SST.100001-109912.nc", "r") do ds
    global SST = replace(ds["SST"][:, :, 1, :], missing=>NaN)
    global Nx, Ny, Nt = size(SST)

    (Nt%12 == 0) || throw(ErrorException("Time is not multiple of 12"))

    global nyears = Int(Nt / 12)

end

Dataset("transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.SHF.100001-109912.nc", "r") do ds
    global SHF = replace(ds["SHF"][:], missing=>NaN)
end


Dataset("transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.HMXL.100001-109912.nc", "r") do ds
    global HMXL = replace(ds["HMXL"][:], missing=>NaN) / 100.0
end

Dataset("transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.TEMP.100001-109912.nc", "r") do ds
    global TEMP = replace(ds["TEMP"][:, :, :, 1], missing=>NaN)

    z_top = - replace(ds["z_w_top"][:], missing=>NaN) / 100.0
    z_bot = - replace(ds["z_w_bot"][:], missing=>NaN) / 100.0
    
    global zs = [z_top..., z_bot[end]]
    global hs = z_top - z_bot
    global Nz = length(zs) - 1


    # Construct the last value of temperature
    global BOT_TEMP = TEMP[:, :, 1]

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

function findZInd(z)
    (z == 0.0) && (return 1)
    for k=1:Nz
        (zs[k] > z >= zs[k+1]) && (return k)
    end
end

function getTd(z, i, j)
    T = TEMP[i, j, findZInd(z)]
    return ( isfinite(T) ) ? T : BOT_TEMP[i, j]
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

Qs_EntSOM15L = zeros(Float64, Nx, Ny, 12)
hs_EntSOM15L = zeros(Float64, Nx, Ny, 12) 

Qs_SOM       = zeros(Float64, Nx, Ny, 12)
hs_SOM       = zeros(Float64, Nx, Ny, 12)


Qs_EntSOM15L .= NaN
hs_EntSOM15L .= NaN
Qs_SOM       .= NaN
hs_SOM       .= NaN

# SOM
println("Doing SOM...")
for i=1:Nx, j=1:Ny

    isnan(SST[i, j, 1]) && continue

    Q = view(Qs_SOM, i, j, :)

    T      = SST[i, j, 13:Nt-12]
    T_next = SST[i, j, 14:Nt-11]

    h = mean(HMXL[i, j, 13:Nt-12])
    
    F      = SHF[i, j, 13:Nt-12]
    F_next = SHF[i, j, 14:Nt-11]

    Q .= 0.0
    for t=1:length(T)
        m = (t-1)%12+1
        
        Q[m] += - (F[t] + F_next[t]) / 2.0
        Q[m] += ρ * c_p * (T_next[t] - T[t]) * h / Δts[m]

    end

    Q[:] /= (nyears - 2)   # remember we discard the first and last year
    Q[:] = (Q + circshift(Q, 1) ) / 2.0

    if any(isnan.(Q))
        println(format("Position (i, j) = ({:d}, {:d}) has NaN Qflux!", i, j))
    end

    hs_SOM[i, j, :] .= h
end


# EntSOM15L
println("Doing EntSOM15L...")
for i=1:Nx, j=1:Ny

    isnan(SST[i, j, 1]) && continue

    Q = view(Qs_EntSOM15L, i, j, :)

    T      = SST[i, j, 13:Nt-12]
    T_next = SST[i, j, 14:Nt-11]

    #h  =  HMXL[i, j, 14:Nt-11] + HMXL[i, j, 13:Nt-12] ) / 2.0
    h       =   HMXL[i, j, 13:Nt-12]
    h_next  =   HMXL[i, j, 14:Nt-11]
    Δh      = h_next - h
    
    F      = SHF[i, j, 13:Nt-12]
    F_next = SHF[i, j, 14:Nt-11]

    Q .= 0.0
    for t=1:length(T)
        m = (t-1)%12+1
        
        Q[m] += - (F[t] + F_next[t]) / 2.0
        Q[m] += ρ * c_p * (T_next[t] - T[t]) * (h[t] + h_next[t]) / 2.0 / Δts[m]

        if Δh[t] > 0.0
            Q[m] += ρ * c_p * (
                (T_next[t] - getTd(-h_next[t], i, j))
              + (T[t]      - getTd(-h[t]     , i, j))
            ) * Δh[t]  / 2.0 / Δts[m]
        end

    end

    Q[:] /= (nyears - 2)   # remember we discard the first and last year
    Q[:] = (Q + circshift(Q, 1) ) / 2.0

    if any(isnan.(Q))
        println(format("Position (i, j) = ({:d}, {:d}) has NaN Qflux!", i, j))
    end

    hs_EntSOM15L[i, j, :] = mean(reshape(h, 12, :), dims=2)[:, 1]
end

println("Output file...")

Dataset(out_file, "c") do ds

    defDim(ds,"time", 12)
    defDim(ds,"Nx", Nx)
    defDim(ds,"Ny", Ny)

    for o in (
        [
            "h_SOM", hs_SOM, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "qflux_SOM", Qs_SOM, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            )
        ],

        [
            "h_EntSOM15L", hs_EntSOM15L, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "qflux_EntSOM15L", Qs_EntSOM15L, ("Nx", "Ny", "time"), Dict(
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
