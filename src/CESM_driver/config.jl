configs = Dict(
    "wdir"        => pwd(),
    "caseroot"    => pwd(),
    "domain_file" => "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc",
    "short_term_archive_dir" => pwd(),
    "long_term_archive_dir"  => pwd(),
    "enable_short_term_archive" => false,
    "enable_long_term_archive"  => false,
    "monthly_record"  => true,
    "yearly_snapshot"  => true,
    "short_term_archive_list" => "SMARTSLAB_short_term_archive_list.txt",
)

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "Configuration file."
            arg_type = String
       
        "--core"
            help = "Core of the model."
            arg_type = String
        #=
        "--init-file"
            help = "Initial profile of ocean. If not provided then a default profile will be used."
            arg_type = String
        =#
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
    println("===== Load config file: ", config_file, " =====")
    include(config_file)
end

if !(@isdefined OMMODULE)
    core = parsed_args["core"]
    if core == nothing
        throw(ErrorException("Core ocean module is not provided. Please set --core option, or define OMMODULE in configuration file."))
    else
        modulename = "CESMCORE_" * core
        core_file = joinpath(dirname(@__FILE__), "cores", core, modulename * ".jl")
        println("Selected core: ", core, " => ", core_file )
        include(core_file)
        OMMODULE = getfield(Main, Symbol(modulename))
    end 
end

println("===== Defining variables BEGIN =====")
print(json(configs, 4))
println("===== Defining variables END =====")



#=
init_file = parsed_args["init-file"]
if init_file != nothing
    init_file = normpath( (isabspath(init_file)) ? init_file : joinpath(pwd(), init_file) )
    println("Ocean init file : ", init_file)

    if !isfile(init_file)
        throw(ErrorException("File missing: ", init_file))
    end
end
=#


function appendLine(filename, content)
    open(filename, "a") do io
        write(io, content)
        write(io, "\n")
    end
end

