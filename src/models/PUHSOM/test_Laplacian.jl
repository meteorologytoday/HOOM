include("../../share/constants.jl")
include("../../share/PolelikeCoordinate.jl")
include("../../share/MapInfo.jl")
include("MatrixOperators.jl")
include("Dyn/AdvectionSpeedUpMatrix_dyn.jl")
include("Dyn/DiffusionSolver.jl")


H  = 1000.0;
Δt = 86400.0;

n = 6;

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv0.9x1.25.nc"
hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx1v6.nc"



mi = ModelMap.MapInfo{Float64}(hrgrid_file)
#mi.mask = 1.0 .- mi.mask


gi = PolelikeCoordinate.CurvilinearSphericalGridInfo(;
    R=Re,
    Ω=Ωe,
    Nx=mi.nx,
    Ny=mi.ny,
    c_lon=mi.xc,
    c_lat=mi.yc,
    vs_lon=mi.xv,
    vs_lat=mi.yv,
    area=mi.area,
    angle_unit=:deg,
)

#=
println("Create Gridinfo");
gi = PolelikeCoordinate.RegularCylindricalGridInfo(;
    R = 6371e3,
    Ω = 2π/86400.0,
    Nx = Int(360 / n),
    Ny = Int(180 / n),
    Ly = 100e3 * 100,
    lat0 = 0.0,
    β    = 100/(100e3*100),
);
=#



println("Making DynamicAdvSpeedUpMatrix")
M = DynamicAdvSpeedUpMatrix(;
    gi    = gi,
    Nz    = 1,
    mask2 = mi.mask,
);

println("Making DiffusionSolver") 
DS = DiffusionSolver(;
    gi = gi,
    M  = M,
    D  = 30000.0,
    Δt = 86400.0,
    mask2 = mi.mask,
)



f          =  sin.(2*gi.c_lon) .* cos.(gi.c_lat);
Lapf_true  =  - sin.(2 * gi.c_lon) .* ( cos.(2 * gi.c_lat) .+ 4.0 ) ./ ( gi.R^2 .* cos.(gi.c_lat));

f         = reshape(M.filter_T * f[:]        , size(f)...)
Lapf_true = reshape(M.filter_T * Lapf_true[:], size(f)...)

Lapf  = reshape(M.T_Lap_T * f[:], size(f)...)


f_U = copy(f)
f_V = hcat( copy(f), f[:, end] )

Lapf_U  = reshape(M.U_Lap_U * f_U[:], size(f_U)...)
Lapf_V  = reshape(M.V_Lap_V * f_V[:], size(f_V)...)

steps=300

rec_f_U = zeros(Float64, size(f_U)..., steps)
rec_f_U[:, :, 1] = f_U
for step=2:steps
    solveDiffusion!(
        DS, :U,
        view(rec_f_U, :, :, step-1),
        view(rec_f_U, :, :, step),
    )
end





using NCDatasets
Dataset("output_test_Laplacian.nc", "c") do ds

    defDim(ds, "Nx", gi.Nx)
    defDim(ds, "Ny", gi.Ny)
    defDim(ds, "Nyp1", gi.Ny+1)
    defDim(ds, "steps", steps)

    for (varname, vardata, vardim, attrib) in [
        ("f",          f,         ("Nx", "Ny"), Dict()),
        ("Lapf_true",  Lapf_true, ("Nx", "Ny"), Dict()),
        ("Lapf",       Lapf,      ("Nx", "Ny"), Dict()),
        ("Lapf_U",     Lapf_U,    ("Nx", "Ny"), Dict()),
        ("Lapf_V",     Lapf_V,    ("Nx", "Nyp1"), Dict()),
        ("f_U_diffuse", rec_f_U,  ("Nx", "Ny", "steps"), Dict()),
    ]

        if ! haskey(ds, varname)
            var = defVar(ds, varname, Float64, vardim)
            var.attrib["_FillValue"] = 1e20
        end

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end

#=



using LinearAlgebra
@time F = lu(cM.cT_Lap_cT)

f         =   cos.(π/gi.Ly * gi.c_y) .* sin.(2*gi.c_lon);
dfdx_true =   cos.(π/gi.Ly * gi.c_y) .* cos.(2*gi.c_lon) * 2 / gi.R;
dfdy_true = - sin.(π/gi.Ly * gi.c_y) .* sin.(2*gi.c_lon) * π / gi.Ly;
Lapf_true =   - f .* ( (2/ gi.R)^2 + (π/gi.Ly)^2 );

f = f[:]

dfdx = cM.M.U_∂x_T * f
dfdy = cM.M.V_∂y_T * f

Lapf = cM.T_Lap_T * f

solve_Lapf = cM.T_send_cT * (F \ (cM.cT_send_T * Lapf))

reshape2 = (m,) -> reshape(m, gi.Nx, :)

f    = reshape2(f)
Lapf = reshape2(Lapf)
solve_Lapf = reshape2(solve_Lapf)
dfdx = reshape2(dfdx)
dfdy = reshape2(dfdy)

using NCDatasets
Dataset("output.nc", "c") do ds

    defDim(ds, "Nx", gi.Nx)
    defDim(ds, "Ny", gi.Ny)
    defDim(ds, "Nyp1", gi.Ny+1)

    for (varname, vardata, vardim, attrib) in [
        ("f",  f, ("Nx", "Ny"), Dict()),
        ("dfdx",  dfdx, ("Nx", "Ny"), Dict()),
        ("dfdy",  dfdy, ("Nx", "Nyp1"), Dict()),
        ("Lapf",  Lapf, ("Nx", "Ny"), Dict()),
        ("solve_Lapf",  solve_Lapf, ("Nx", "Ny"), Dict()),
        ("dfdx_true",  dfdx_true, ("Nx", "Ny"), Dict()),
        ("dfdy_true",  dfdy_true, ("Nx", "Ny"), Dict()),
        ("Lapf_true",  Lapf_true, ("Nx", "Ny"), Dict()),

    ]

        if ! haskey(ds, varname)
            var = defVar(ds, varname, Float64, vardim)
            var.attrib["_FillValue"] = 1e20
        end

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end
=#
