include("../CoordTrans/CoordTrans.jl")

using .CoordTrans

if length(ARGS) >= 3
    in_filename = ARGS[1]
    out_filename = ARGS[2]
    wgt_filename = ARGS[3]
end



CoordTrans.convertFile(
    in_filename,
    out_filename,
    wgt_filename,
    varnames=("T_ML",);
    xdim = "Nx",
    ydim = "Ny",
    tdim = "time",
)
