include("DisplacedPoleCoordinate.jl")
include("MapInfo.jl")

using .DisplacedPoleCoordinate
using .ModelMap
using NCDatasets

R =  1.0#6371000.0

#mi = ModelMap.MapInfo{Float64}("domain.ocn.gx3v7.120323.nc")
mi = ModelMap.MapInfo{Float64}("domain.ocn.gx1v6.090206.nc")
#mi = ModelMap.MapInfo{Float64}("domain.lnd.fv4x5_gx3v7.091218.nc")
gi = DisplacedPoleCoordinate.GridInfo(R, mi.nx, mi.ny, mi.xc, mi.yc, mi.xv, mi.yv)

u_Eu = zeros(Float64, gi.Nx, gi.Ny)
v_Eu = zeros(Float64, gi.Nx, gi.Ny)

u_Dp = zeros(Float64, gi.Nx, gi.Ny)
v_Dp = zeros(Float64, gi.Nx, gi.Ny)

u2_Eu = zeros(Float64, gi.Nx, gi.Ny)
v2_Eu = zeros(Float64, gi.Nx, gi.Ny)

div = zeros(Float64, gi.Nx, gi.Ny)

div_analytic = zeros(Float64, gi.Nx, gi.Ny)

for i=1:gi.Nx, j=1:gi.Ny

    λ = mi.xc[i, j] * π/180.0
    ϕ = mi.yc[i, j] * π/180.0

    u_Eu[i, j] = cos.(ϕ) .* sin.(λ) #sin.(    λ) .* sin.(2.0*ϕ)
    v_Eu[i, j] = cos.(ϕ)            #cos.(2.0*λ) .* sin.(    ϕ)

    div_analytic[i, j] = (cos.(λ) - 2.0 * sin.(ϕ)) / R

end

DisplacedPoleCoordinate.project!(gi, u_Eu, v_Eu, u_Dp, v_Dp, direction=:Forward)
DisplacedPoleCoordinate.project!(gi, u_Dp, v_Dp, u2_Eu, v2_Eu, direction=:Backward)


#mi.mask[ mi.mask .== 0.0 ] .= 1.0
#mi.mask[ mi.mask .== 0.0 ] .= 1.0

DisplacedPoleCoordinate.DIV!(gi, u_Dp, v_Dp, div, mi.mask)


div[mi.mask.==0] .= NaN
div_analytic[mi.mask.==0] .= NaN


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
        (div_analytic, "div_analytic", ("ni", "nj"), nothing),
        (div - div_analytic, "div_diff", ("ni", "nj"), nothing),
        (u_Eu.^2  +  v_Eu.^2, "u_Eu_square", ("ni", "nj"), nothing),       
        (u2_Eu.^2 + v2_Eu.^2, "u2_Eu_square", ("ni", "nj"), nothing),       
 
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
