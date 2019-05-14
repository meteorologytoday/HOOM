include("DisplacedPoleCoordinate.jl")
include("MapInfo.jl")

using .DisplacedPoleCoordinate
using .ModelMap
using NCDatasets

mi = ModelMap.MapInfo{Float64}("domain.ocn.gx3v7.120323.nc")
gi = DisplacedPoleCoordinate.GridInfo(mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv)

u_Eu = zeros(Float64, gi.Nx, gi.Ny)
v_Eu = zeros(Float64, gi.Nx, gi.Ny)

u_Dp = zeros(Float64, gi.Nx, gi.Ny)
v_Dp = zeros(Float64, gi.Nx, gi.Ny)

u2_Eu = zeros(Float64, gi.Nx, gi.Ny)
v2_Eu = zeros(Float64, gi.Nx, gi.Ny)



u_Eu .= 10.0

DisplacedPoleCoordinate.project!(gi, u_Eu, v_Eu, u_Dp, v_Dp, direction=:Forward)
DisplacedPoleCoordinate.project!(gi, u_Dp, v_Dp, u2_Eu, v2_Eu, direction=:Backward)



Dataset("dpc.nc", "c") do ds

    defDim(ds, "ni", mi.nx)
    defDim(ds, "nj", mi.ny)
    defDim(ds, "nv", 4)
    defDim(ds, "vec", 2)

    for (var, varname, dims) in [
        (gi.dσ, "area", ("ni", "nj")),   
        (gi.dx_w, "dx_w", ("ni", "nj")),
        (gi.dx_c, "dx_c", ("ni", "nj")),       
        (gi.dx_e, "dx_e", ("ni", "nj")),       
        (gi.dy_s, "dy_s", ("ni", "nj")),       
        (gi.dy_c, "dy_c", ("ni", "nj")),       
        (gi.dy_n, "dy_n", ("ni", "nj")),       
        (gi.α,    "alpha", ("ni", "nj")),        
        (gi.cosα,  "cos_alpha", ("ni", "nj")),
        (gi.sinα,  "sin_alpha", ("ni", "nj")),
        (u_Eu, "u_Eu", ("ni", "nj")),       
        (v_Eu, "v_Eu", ("ni", "nj")),       
        (u_Dp, "u_Dp", ("ni", "nj")),       
        (v_Dp, "v_Dp", ("ni", "nj")),       
        (u2_Eu, "u2_Eu", ("ni", "nj")),       
        (v2_Eu, "v2_Eu", ("ni", "nj")),       
    ]

        println("Output: ", varname)
        println(dims, "; ", typeof(var))
        v = defVar(ds, varname, Float64, dims)
        v.attrib["_FillValue"] = mi.missing_value
        v[:] = var
         
    end




end
