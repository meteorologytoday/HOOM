include("../MLMML/SSM.jl")
include("julia_lib/Mailbox.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")


using ArgParse
using Formatting
using Printf
using JSON
using .Mailbox
using .BinaryIO
using .SSM
using .NetCDFIO
using Statistics: std, mean

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

include("SSM_config.jl")
include("SSM_driver.jl")
