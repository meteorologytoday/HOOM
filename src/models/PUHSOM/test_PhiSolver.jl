include("../../../share/PolelikeCoordinate.jl")
include("AdvectionSpeedUpMatrix_dyn.jl")
include("PhiSolver.jl")

H  = 1000.0;
Δt = 86400.0;

n = 6;

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

mask2 = ones(Float64, gi.Nx, gi.Ny);

#mask2[[20, 50], :] .= 0

#=
println("Making DynamicAdvSpeedUpMatrix")
M = DynamicAdvSpeedUpMatrix(;
    gi    = gi,
    Nz    = 3,
    mask2 = mask2,
);
=#

println("Making ΦSolver")
@time cM = ΦSolver(;
    gi    = gi,
    mask2 = mask2,
    α     = 1 / (Δt^2 * H)
);




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

