overwrite_configs = Dict()
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
    println("===== Load config file: ", config_file, " =====")
    include(config_file)
end

if !(@isdefined OMMODULE)
    core_name = parsed_args["core"]
    if core_name == nothing
        throw(ErrorException("Core ocean module is not provided. Please set --core option, or define OMMODULE in configuration file."))
    else
        module_name = "CESMCORE_" * core_name
        #module_symb = Symbol(module_name)
        core_file = joinpath(dirname(@__FILE__), "cores", core_name, module_name * ".jl")
        println("Selected core: ", core_name, " => ", core_file )
        #=
        for pid in procs()
            remotecall_fetch(include, pid, core_file)
            remotecall_fetch(eval, pid, :(using .$(module_symb))) 
        end
        =#
        include(core_file)
        OMMODULE = getfield(Main, Symbol(module_name))
    end 
end

println("===== Defining variables BEGIN =====")
for (k, v) in overwrite_configs
    if k in keys(configs)
        println("Overwrite config ", k, "...")
    else
        println("Add config ", k, "...")
    end

    configs[k] = v
end
print(json(configs, 4))
println("===== Defining variables END =====")



init_file = parsed_args["init-file"]
if init_file != nothing
    init_file = normpath( (isabspath(init_file)) ? init_file : joinpath(pwd(), init_file) )
    println("Ocean init file : ", init_file)

    if !isfile(init_file)
        throw(ErrorException("File missing: ", init_file))
    end
end



