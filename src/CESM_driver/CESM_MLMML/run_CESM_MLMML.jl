using ArgParse
using Printf

# ===== READ ALL CONFIG BEGIN =====

include("../default_config.jl")
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "Configuration file."
            arg_type = String

        "--init-file"
            help = "Initial profile of ocean. If not provided then a default profile will be used."
            arg_type = String

    end

    return parse_args(s)
end

parsed_args = parse_commandline()
#for (arg,val) in parsed_args
#    println("$arg  =>  ", repr(val))
#end

config_file = parsed_args["config"]
if config_file != nothing
    config_file = normpath( (isabspath(config_file)) ? config_file : joinpath(pwd(), config_file) )
    println("Load config file: ", config_file)
    include(config_file)
end

init_file = parsed_args["init-file"]
if init_file != nothing
    init_file = normpath( (isabspath(init_file)) ? init_file : joinpath(pwd(), init_file) )
    println("Ocean init file : ", init_file)

    if !isfile(init_file)
        throw(ErrorException("File missing: ", init_file))
    end
end

# ===== READ ALL CONFIG END =====

# ===== Load Core Module =====
include("CESM_CORE_MLMML.jl")
OMMODULE = Main.CESM_CORE_MLMML

# ===== Load General Driver =====
include("../driver.jl")
