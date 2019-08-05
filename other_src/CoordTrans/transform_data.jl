include(joinpath(@__DIR__, "..", "CoordTrans", "CoordTrans.jl"))

using .CoordTrans
using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--w-file"
            help = "Weighting file."
            arg_type = String
            required = true

        "--s-file"
            help = "Input source file."
            arg_type = String
            required = true
       
        "--d-file"
            help = "Output destination file."
            arg_type = String
            required = true

        "--vars"
            help = "Variable names list. They should by dimension 2 (x, y) or 3 (x, y, z) with or without record (time) dimension. Ex: --vars=Ts,Ss,MLD"
            arg_type = String
            required = true

        "--copy-vars"
            help = "Variable names list."
            arg_type = String


        "--x-dim"
            help = "Variable name of x-dimension."
            arg_type = String
            required = true

        "--y-dim"
            help = "Variable name of y-dimension."
            arg_type = String
            required = true

        "--z-dim"
            help = "Variable name of z-dimension."
            arg_type = String

        "--t-dim"
            help = "Variable name of time-dimension."
            arg_type = String
            default = "time"

    end

    return parse_args(ARGS, s)
end

println("Running ", @__FILE__)

parsed = parse_commandline()
print(json(parsed, 4))
    
varnames=collect(split(parsed["vars"], ","; keepempty=false))

CoordTrans.convertFile(
    parsed["s-file"],
    parsed["d-file"],
    parsed["w-file"],
    varnames=varnames;
    xdim = parsed["x-dim"],
    ydim = parsed["y-dim"],
    zdim = parsed["z-dim"],
    tdim = parsed["t-dim"],
    copy_varnames = (parsed["copy-vars"] == nothing) ? (:,) : split(parsed["copy-vars"], ",") ,
)
