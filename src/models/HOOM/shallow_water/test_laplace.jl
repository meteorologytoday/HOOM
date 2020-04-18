include("../../../share/PolelikeCoordinate.jl")
include("AdvectionSpeedUpMatrix_dyn.jl")

gi = PolelikeCoordinate.RegularCylindricalGridInfo(;
    R = 6371e3,
    Ω = 2π/86400.0,
    Nx = 12,
    Ny = 10,
    Ly = 100e3 * 100,
    lat0 = 0.0,
    β    = 100/(100e3*100),
)

M = DynamicAdvSpeedUpMatrix(;
    gi    = gi,
    Nz    = 1,
    mask2 = ones(Float64, gi.Nx, gi.Ny),
)




