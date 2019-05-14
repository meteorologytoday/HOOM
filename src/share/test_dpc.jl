include("DisplacedPoleCoordinate.jl")
include("MapInfo.jl")

using .DisplacedPoleCoordinate
using .ModelMap
using NCDatasets

R = 6371000.0

mi = ModelMap.MapInfo{Float64}("domain.ocn.gx3v7.120323.nc")
#mi = ModelMap.MapInfo{Float64}("domain.ocn.gx1v6.090206.nc")
gi = DisplacedPoleCoordinate.GridInfo(R, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv)

u_Eu = zeros(Float64, gi.Nx, gi.Ny)
v_Eu = zeros(Float64, gi.Nx, gi.Ny)

u_Dp = zeros(Float64, gi.Nx, gi.Ny)
v_Dp = zeros(Float64, gi.Nx, gi.Ny)

u2_Eu = zeros(Float64, gi.Nx, gi.Ny)
v2_Eu = zeros(Float64, gi.Nx, gi.Ny)

div = zeros(Float64, gi.Nx, gi.Ny)


for i=1:gi.Nx, j=1:gi.Ny
    lat = mi.yc[i, j]
    lon = mi.xc[i, j]
    u_Eu[i, j] = 100.0 #* cos(lat * π / 180.0)
end

DisplacedPoleCoordinate.project!(gi, u_Eu, v_Eu, u_Dp, v_Dp, direction=:Forward)
DisplacedPoleCoordinate.project!(gi, u_Dp, v_Dp, u2_Eu, v2_Eu, direction=:Backward)

DisplacedPoleCoordinate.divergence2!(gi, u_Dp, v_Dp, div)

Dataset("dpc3v7.nc", "c") do ds

    defDim(ds, "ni", mi.nx)
    defDim(ds, "nj", mi.ny)
    defDim(ds, "nv", 4)
    defDim(ds, "vec", 2)

    for (var, varname, dims, attrib) in [
        (gi.dσ, "area", ("ni", "nj"), nothing),   
        (gi.dx_w, "dx_w", ("ni", "nj"), nothing),
        (gi.dx_c, "dx_c", ("ni", "nj"), nothing),       
        (gi.dx_e, "dx_e", ("ni", "nj"), nothing),       
        (gi.dy_s, "dy_s", ("ni", "nj"), nothing),       
        (gi.dy_c, "dy_c", ("ni", "nj"), nothing),       
        (gi.dy_n, "dy_n", ("ni", "nj"), nothing),       
        (gi.α,    "alpha", ("ni", "nj"), nothing),        
        (gi.cosα,  "cos_alpha", ("ni", "nj"), nothing),
        (gi.sinα,  "sin_alpha", ("ni", "nj"), nothing),
        (u_Eu, "u_Eu", ("ni", "nj"), nothing),       
        (v_Eu, "v_Eu", ("ni", "nj"), nothing),       
        (u_Dp, "u_Dp", ("ni", "nj"), nothing),       
        (v_Dp, "v_Dp", ("ni", "nj"), nothing),       
        (u2_Eu, "u2_Eu", ("ni", "nj"), nothing),       
        (v2_Eu, "v2_Eu", ("ni", "nj"), nothing),       
        (div, "div", ("ni", "nj"), nothing),
        (mi.xc, "xc", ("ni", "nj"), Dict("long_name"=>"longitude of grid cell center", "units"=>"degrees_east")),
        (mi.yc, "yc", ("ni", "nj"), Dict("long_name"=>"latitude of grid cell center", "units"=>"degrees_north")),
    ]

        println("Output: ", varname)
        println(dims, "; ", typeof(var))
        v = defVar(ds, varname, Float64, dims)
        v.attrib["_FillValue"] = mi.missing_value

        if attrib != nothing
            for (key, val) in attrib
                v.attrib[key] = val
            end
        end

        v[:] = var
         
    end




end
