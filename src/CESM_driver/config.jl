overwrite_configs = Dict()
configs = Dict(
    :casename    => "casename",
    :substeps    => 1,                 # This controls how many steps will occur for each CESM coupling. Example: ocean couple to atmosphere every 24 hours but itself steps every 3 hours. This means we would expect `Δt` = 86400, and we set `substeps` = 8.
    :caseroot                  => pwd(),
    :domain_file               => "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc",
    :short_term_archive_dir    => pwd(),
    :long_term_archive_dir     => pwd(),
    :enable_short_term_archive => false,
    :enable_long_term_archive  => false,
    :daily_record              => false,
    :monthly_record            => true,
    :yearly_snapshot           => true,
    :short_term_archive_list   => "SMARTSLAB_short_term_archive_list.txt",
    :rpointer_file             => "rpointer.ssm_ocn",
    :wdir                      => pwd(),
    :timeout                   => 60.0 * 20, 
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

if ! ( "tmp_folder" in keys(configs) )
    
    configs[:tmp_folder] = joinpath(configs[:caserun], "x_tmp")
    
end


print(json(configs, 4))
println("===== Defining variables END =====")
