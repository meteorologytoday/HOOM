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


function runModel(
    OMMODULE      :: Any,
    coupler_funcs :: Any, 
)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println("Hello world, I am $(rank) of $(MPI.Comm_size(comm))")
    MPI.Barrier(comm)

    is_master = rank == 0

    t_start, read_restart, configs = coupler_funcs.before_model_init!()

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
        casename     = configs[:casename],
        clock        = clock,
        configs      = configs,
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

            end
            coupler_funcs.after_model_run!(OMMODULE, OMDATA)
            
            advanceClock!(clock, Δt)

        elseif stage == :END
            is_master && println("stage == :END. Break loop now.")
            break
        end
    end
        
    coupler_funcs.finalize!(OMMODULE, OMDATA)
    
    is_master && println("Program Ends.")

end
  
