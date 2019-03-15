

#=

# Example:
caseroot    = "/home/tienyiah/projects/cesm1_test/CICE_f45"
wdir        = "/home/tienyiah/cesm1/scratch/CICE_f45/run"
domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5
max_try = 60
output_record_length = 365

=#


# ===== Configuration BEGIN =====
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "Configuration file."
            arg_type = String
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println("##### SSM parsed args beg #####")
for (arg,val) in parsed_args
    println("$arg  =>  ", repr(val))
end
println("##### SSM parsed args end #####")


config_file = parsed_args["config"]
if config_file != nothing
    config_file = normpath( (isabspath(config_file)) ? config_file : joinpath(pwd(), config_file) )
    println("Load config file: ", config_file)
    include(config_file)
end
# ===== Configuration END =====

include("../MLMML/SSM.jl")
OMMODULE = Main.SSM
OMDATA   = OMMODULE.MLMML_DATA()

include("driver.jl")
