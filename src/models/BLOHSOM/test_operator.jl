include("../../share/constants.jl")
include("../../share/GridFiles.jl")
include("../../share/PolelikeCoordinate.jl")

include("MatrixOperators.jl")
include("Dyn/AdvectionSpeedUpMatrix_dyn.jl")



gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = 20,
        Ny   = 10,
        Ly   = 100e3 * 60,
        lat0 = 0.0 |> deg2rad,
        β    = Ωe / Re,
)

gi = PolelikeCoordinate.genGridInfo(gf);

xcutoff = 2
ycutoff = 3

gf.mask                       .= 1
gf.mask[:, 1:ycutoff]         .= 0 
gf.mask[:, end-ycutoff+1:end] .= 0 

gf.mask[1:xcutoff, :]         .= 0 
gf.mask[end-xcutoff+1:end, :] .= 0 

println("Making DynamicAdvSpeedUpMatrix")
M = DynamicAdvSpeedUpMatrix(;
    gi    = gi,
    Nz    = 1,
    mask2 = gf.mask,
);

f_U = rand(gi.Nx, gi.Ny)
f_V = rand(gi.Nx, gi.Ny+1)

half = floor(Int64, gi.Ny/2)
f_V[:, 1:half] = f_V[:, end:-1:end-half+1]


f_T = cos.(gi.c_y / gi.Ly * 2*π);

#f_U = reshape( M.filter_U * f_U[:], gi.Nx, gi.Ny  )
#f_V = reshape( M.V_∂filter_V * f_V[:], gi.Nx, gi.Ny+1)
#f_T = reshape( M.T_DIVy_V * f_V[:], gi.Nx, gi.Ny  )
f_V = reshape( M.V_∂y_T * f_T[:], gi.Nx, gi.Ny+1 )


