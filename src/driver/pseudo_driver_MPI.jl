using CFTime
using Dates
using Formatting
using JSON
using Distributed
using SharedArrays
using MPI

if !(:ModelClockSystem in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "..", "share", "ModelClockSystem.jl")))
end
using .ModelClockSystem

if !(:ConfigCheck in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "..", "share", "ConfigCheck.jl")))
end
using .ConfigCheck



function runModel(
    OMMODULE      :: Any,
    coupler_funcs :: Any, 
)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    writeLog("===== [ Master Created ] =====")

    MPI.Barrier(comm)

    is_master = rank == 0

    t_start, read_restart, config = coupler_funcs.before_model_init!()

    if is_master
        writeLog("Validate driver config.")
        config[:DRIVER] = validateConfigEntries(config[:DRIVER], getDriverConfigDescriptor())
    end

    writeLog("Setting working directory: {:s}", config[:DRIVER][:caserun])
    cd(config[:DRIVER][:caserun])

    # copy the start time
    beg_datetime = t_start + Dates.Second(0)

    if is_master
        println(format("Begin datetime: {:s}", dt2str(beg_datetime) ))
        println(format("Read restart  : {}", read_restart))
    end

    # Construct model clock
    clock = ModelClock("Model", beg_datetime)

    # initializing
    is_master && println("===== INITIALIZING MODEL: ", OMMODULE.name , " =====")
    
    OMDATA = OMMODULE.init(
        casename     = config[:DRIVER][:casename],
        clock        = clock,
        config      = config,
        read_restart = read_restart,
    )

    coupler_funcs.after_model_init!(OMMODULE, OMDATA)
    
    is_master && println("Ready to run the model.")
    step = 0
    while true

        step += 1
        
        is_master && println(format("Current time: {:s}", clock2str(clock)))

        stage, Δt, write_restart = coupler_funcs.before_model_run!(OMMODULE, OMDATA)

        if stage == :RUN 

            cost = @elapsed let

                OMMODULE.run!(
                    OMDATA;
                    Δt = Δt,
                    write_restart = write_restart,
                )

                MPI.Barrier(comm)

            end
            is_master && println(format("Computation cost: {:f} secs.", cost))
            coupler_funcs.after_model_run!(OMMODULE, OMDATA)
            
            advanceClock!(clock, Δt)

            

        elseif stage == :END
            is_master && println("stage == :END. Break loop now.")
            break
        end
    end
    
    OMMODULE.final(OMDATA) 
    coupler_funcs.final(OMMODULE, OMDATA)
  
    if is_master
        archive(joinpath(
            config[:DRIVER][:caserun],
            config[:DRIVER][:archive_list],
        ))
    end
 
    is_master && println("Program Ends.")

end

function getDriverConfigDescriptor()

    return [
            ConfigEntry(
                :casename,
                :required,
                [String,],
            ),

            ConfigEntry(
                :caseroot,
                :required,
                [String,],
            ),

            ConfigEntry(
                :caserun,
                :required,
                [String,],
            ),

            ConfigEntry(
                :archive_root,
                :required,
                [String,],
            ),

            ConfigEntry(
                :archive_list,
                :optional,
                [String,],
                "archive_list.txt",
            ),
   ]
end

function archive(
    archive_list_file :: String,
)

    println("===== Archiving files BEGIN =====")
    
    for line in eachline(archive_list_file)

        args = split(line, ",")

        if length(args) == 0
            continue
        end
      
        action = args[1]
        args = args[2:end]

        if action in ["mv", "cp"]

            fname, src_dir, dst_dir = args

            if ! isdir(dst_dir)
                mkpath(dst_dir)
            end
 
            src_file = joinpath(src_dir, fname)
            dst_file = joinpath(dst_dir, fname)

            if isfile(src_file)

                if action == "mv"
                    mv(src_file, dst_file, force=true)
                    println(format("Moving file: {:s} ( {:s} => {:s} )", fname, src_dir, dst_dir))
                elseif action == "cp"
                    cp(src_file, dst_file, force=true)
                    println(format("Copying file: {:s} ( {:s} => {:s} )", fname, src_dir, dst_dir))
                end

            else
                println("File does not exist: ", src_file)
            end

        elseif action == "rm"
            fname, fdir = args
            rm(joinpath(fdir, fname), force=true)
            println(format("Removing file: {:s} in {:s}", fname, fdir))
        else
            throw(ErrorException(format("Unknown action in archive list: {:s}", action)))
        end

    end

    println("===== Archiving files END =====")

end
