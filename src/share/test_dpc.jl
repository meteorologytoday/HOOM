include("DisplacedPoleCoordinate.jl")
include("MapInfo.jl")

using .DisplacedPoleCoordinate
using .ModelMap
using NCDatasets

mi = ModelMap.MapInfo{Float64}("domain.ocn.gx3v7.120323.nc")
gi = DisplacedPoleCoordinate.GridInfo(mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv)

Dataset("dpc.nc", "c") do ds

    defDim(ds, "ni", mi.nx)
    defDim(ds, "nj", mi.ny)
    defDim(ds, "nv", 4)

    for (var, varname, dims) in [
        (gi.dσ, "area", ("ni", "nj")),   
        (gi.dx_w, "dx_w", ("ni", "nj")),
        (gi.dx_c, "dx_c", ("ni", "nj")),       
        (gi.dx_e, "dx_e", ("ni", "nj")),       
        (gi.dy_s, "dy_s", ("ni", "nj")),       
        (gi.dy_c, "dy_c", ("ni", "nj")),       
        (gi.dy_n, "dy_n", ("ni", "nj")),       
        (gi.α,    "alpha", ("ni", "nj")),        
    ]

        println("Output: ", varname)
        println(dims, "; ", typeof(var))
        v = defVar(ds, varname, Float64, dims)
        v.attrib["_FillValue"] = mi.missing_value
        v[:] = var
         
    end




end
