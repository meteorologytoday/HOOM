include("../../src/share/DisplacedPoleCoordinate.jl")
include("../../src/share/MapInfo.jl")
include("../../src/share/constants.jl")

using Formatting
using NCDatasets
using Statistics

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
doy = sum(dom)
doy == 365 || throw(ErrorException("Sum of dom is not 365"))

Δts = (dom + circshift(dom, -1))/2.0 * 86400.0

println("Δts = ", Δts)

in_SST  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SST.100001-109912.nc"
in_SHF  = "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc" 
in_HMXL = "b.e11.B1850C5CN.f09_g16.005.pop.h.HMXL.100001-109912.nc"
in_QFLUX = "b.e11.B1850C5CN.f09_g16.005.pop.h.QFLUX.100001-109912.nc"
in_MELTH_F = "b.e11.B1850C5CN.f09_g16.005.pop.h.MELTH_F.100001-109912.nc"


in_zcoord = "b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc"
gridinfo_file = "domain.ocn.gx1v6.090206.nc"

out_file = "forcing.gx1v6.nc"

missing_value = 1e20

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

Dataset(in_QFLUX, "r") do ds
    global QFLUX = replace(ds["QFLUX"][:], missing=>NaN)
end

Dataset(in_MELTH_F, "r") do ds
    global MELTH_F = replace(ds["MELTH_F"][:], missing=>NaN)
end

SURF = SHF + QFLUX

mi = ModelMap.MapInfo{Float64}(gridinfo_file)
gi = DisplacedPoleCoordinate.GridInfo(Re, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv, mi.area; angle_unit=:deg)


qflux_T = zeros(Float64, Nx, Ny, 12)
hs      = zeros(Float64, Nx, Ny, 12)

println("Calculating MLD")
# MLD
for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
        hs[i, j, :] .= NaN
        continue
    end 
    hs_variate = mean(reshape(view(HMXL, i, j, :), 12, :), dims=2)[:, 1]
    hs[i, j, :] .= sum(hs_variate .* dom) / doy

end

for t = 13:Nt - 12

    print(format("\rCalculate New Qflux Progress: {:d}", t))

    # month
    m = (t-1)%12+1

    for i=1:Nx, j=1:Ny

        isnan(SST[i, j, 1]) && continue
        
        h = hs[i, j, 1]
        tmp_dTdts = (SST[i, j, t+1] - SST[i, j, t]) / Δts[m] * h * ρ * c_p
        tmp_Fs    = (SURF[i, j, t] + SURF[i, j, t+1]) / 2.0
        
        tmp_qflux_T = - ( tmp_dTdts - tmp_Fs )
        qflux_T[i, j, m]     += tmp_qflux_T

    end

end

needed_vars = (qflux_T,)

for var in needed_vars
    var ./= (nyears - 2.0)

    # What we derived is the flux between middle of months
    for i=1:Nx, j=1:Ny
        var[i, j, :] = (var[i, j, :] + circshift(var[i, j, :], 1)) / 2.0
    end
end

# mask out data
for i=1:Nx, j=1:Ny

    if isnan(SST[i, j, 1])
         
        for var in needed_vars
            var[i, j, :] .= NaN
        end
 
    end

end

println("Output file...")

Dataset(out_file, "c") do ds

    defDim(ds,"time", 12)
    defDim(ds,"Nx", Nx)
    defDim(ds,"Ny", Ny)

    for o in (
        [
            "h_ML", hs, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "qflux_T", qflux_T, ("Nx", "Ny", "time"), Dict(
            "long_name"=>"Qflux in energy flux form",
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
