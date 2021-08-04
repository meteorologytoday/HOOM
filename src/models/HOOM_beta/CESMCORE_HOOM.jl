
if ! ( :HOOM in names(Main) )
    include(joinpath(@__DIR__, "HOOM.jl"))
end

module CESMCORE_HOOM

    include(joinpath(@__DIR__, "..", "..", "share", "CheckDict.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "AppendLine.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "Log.jl"))
    include(joinpath(@__DIR__, "DataManager.jl"))
    include(joinpath(@__DIR__, "Parallelization.jl"))

    using MPI
    using Dates
    using Formatting

    using ..HOOM
    using .DataManager
    using .Parallization

    using ..ModelClockSystem
    using NCDatasets


    name = "HOOM_beta"

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct HOOM_DATA
        casename    :: AbstractString
        mb          :: HOOM.ModelBlock
        clock       :: ModelClock

        x2o         :: Dict
        o2x         :: Dict

        config     :: Dict

        data_table  :: DataTable
        recorders   :: Union{Dict, Nothing}

        jdi        :: JobDistributionInfo

        sync_data   :: Dict
    end


    function init(;
        casename     :: AbstractString,
        clock        :: ModelClock,
        config      :: Dict,
        read_restart :: Bool,
    )

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        comm_size = MPI.Comm_size(comm)
        
        if comm_size == 1
            throw(ErrorException("Need at least 2 processes to work"))
        end

        is_master = ( rank == 0 )

        local master_mb = nothing
        local master_ev = nothing

        if is_master

            cfg_desc = HOOM.getConfigDescriptor()

            misc_config = HOOM.validateConfigEntries(config[:MODEL_MISC], cfg_desc[:MODEL_MISC])
            init_file = misc_config[:init_file]

            # If `read_restart` is true then read restart file: config[:rpointer_file]
            # If not then initialize ocean with default profile if `initial_file`
            # is "nothing", with `init_file` if it is nonempty.

            if read_restart

                println("`read_restart` is on. Look for rpointer file...")

                if !isfile(init_file_config[:rpointer_file])
                    throw(ErrorException(init_file_config[:rpointer_file] * " does not exist!"))
                end
                
                println("Going to read restart pointer file: ", init_file_config[:rpointer_file])

                open(init_file_config[:rpointer_file], "r") do file
                    init_file = readline(file)
                end

                if !isfile(init_file)
                    throw(ErrorException(format("Initial file \"{:s}\" does not exist!", init_file)))
                end
     
            end

            core_config = HOOM.validateConfigEntries(config[:MODEL_CORE], cfg_desc[:MODEL_CORE])

            if init_file != ""

                println("Initial ocean with profile: ", init_file)
                println("Initial ocean with domain file: ", core_config[:domain_file])
                master_mb = HOOM.loadSnapshot(init_file; gridinfo_file=core_config[:domain_file])
            
            else

                writeLog("Debugging status. Initialize an empty ocean.")

                

                master_ev = HOOM.Env(core_config)

                
                master_mb = HOOM.ModelBlock(master_ev; init_core=false)
                master_mb.fi.sv[:TEMP][:, :, :]  .= 10
                master_mb.fi.sv[:TEMP][1, 20, :] .= 30

                master_mb.fi.sv[:TEMP][1, :, 21] .= 21
                master_mb.fi.sv[:TEMP][1, :, 22] .= 22
                master_mb.fi.sv[:TEMP][1, :, 23] .= 23
                master_mb.fi.sv[:TEMP][1, :, 24] .= 24
                master_mb.fi.sv[:TEMP][1, :, 25] .= 25
                master_mb.fi.sv[:TEMP][1, :, 26] .= 26
                master_mb.fi.sv[:TEMP][1, :, 27] .= 27
                master_mb.fi.sv[:TEMP][1, :, 28] .= 28

                master_mb.fi.HMXL .= 50.0
                master_mb.fi.HMXL[1, 20:25, 21:24] .= 75.0

                #throw(ErrorException("Variable `init_file` is absent in `config`."))

            end

            #in_flds = ocn.in_flds
            #
            # If it is "datastream", entrainment speed w_e would be 
            # calculated from h given. In fact there is no need
            # to calculate w_e.
            #
            # If it is "prognostic", entrainment speed w_e would be
            # calculated accroding to Niiler-Kraus dynamics.
            #

        end

        # Create plans
        if is_master
            jdi = JobDistributionInfo(nworkers = comm_size - 1, Ny = master_ev.Ny; overlap=3)
            printJobDistributionInfo(jdi)
        end
     
        # First, broadcast ev and create plan 
        master_ev = MPI.bcast(master_ev, 0, comm)
        jdi = JobDistributionInfo(nworkers = comm_size - 1, Ny = master_ev.Ny; overlap=3)



        # Second, create ModelBlocks based on ysplit_info
        if is_master
            my_ev = master_ev
            my_mb = master_mb
        else
            my_ev          = deepcopy(master_ev) 
            my_ev.sub_yrng = getYsplitInfoByRank(jdi, rank).pull_fr_rng
            my_ev.Ny       = length(my_ev.sub_yrng)
            my_mb          = HOOM.ModelBlock(my_ev; init_core = true)
        end

        MPI.Barrier(comm)
        # Third, register all the variables.
        # Fourth, weaving MPI sending relation 

        empty_arr_sT = zeros(Float64, 1, my_ev.Nx, my_ev.Ny)
        x2o = Dict(
            "SWFLX"  => my_mb.fi.SWFLX,
            "NSWFLX" => my_mb.fi.NSWFLX,
            "TAUX_east"    => my_mb.fi.TAUX_east,
            "TAUY_north"   => my_mb.fi.TAUY_north,
            "IFRAC"  => copy(empty_arr_sT),
            "FRWFLX" => copy(empty_arr_sT),
            "VSFLX"  => copy(empty_arr_sT),
            "QFLX_T" => copy(empty_arr_sT),
            "QFLX_S" => copy(empty_arr_sT),
            "MLD"    => copy(empty_arr_sT),
        )


        o2x = Dict(
            "SST"      => view(my_mb.fi.sv[:TEMP], 1, :, :),
            "QFLX2ATM" => my_mb.fi.QFLX2ATM,
        )

        # Create DataTable
        data_table = DataTable(Nz = my_ev.Nz, Nx = my_ev.Nx, Ny = my_ev.Ny)

        
        for (k, (varref, grid_type)) in HOOM.getVariableList(my_mb, :ALL)
            regVariable!(data_table, k, grid_type, varref) 
        end

        # Synchronizing Data
        sync_data_list = Dict(
            :forcing => (
                "SWFLX",
                "NSWFLX",
                "TAUX_east",
                "TAUY_north",
            ),

            :state   => (
                "TEMP",
                "SALT",
                "UVEL",
                "VVEL",
                "WVEL",
                "CHKTEMP",
                "CHKSALT",
                "TAUX",
                "TAUY",
                "ADVT",
                "HMXL",
            ),

            :bnd_state   => (
                "TEMP",
                "SALT",
            )

        )
        
        sync_data = Dict()
        for (k, l) in sync_data_list
            sync_data[k] = Array{DataUnit, 1}(undef, length(l))
            for (n, varname) in enumerate(l)
                sync_data[k][n] = data_table.data_units[varname]
            end
        end


        MD = HOOM_DATA(
            casename,
            my_mb,
            clock,
            x2o,
            o2x,
            config,
            data_table,
            nothing,
            jdi,
            sync_data,
        )


        if is_master

            activated_record = []
            MD.recorders = Dict()
            complete_variable_list = HOOM.getVariableList(my_mb, :ALL)

#            additional_variable_list = HOOM.getVariableList(ocn, :COORDINATE)

            for rec_key in [:daily_record, :monthly_record]
       
                activated = false
 
                println("# For record key: " * string(rec_key))

                var_list = []
                
                if typeof(misc_config[rec_key]) <: Symbol 
                    misc_config[rec_key] = HOOM.getVariableList(my_mb, misc_config[rec_key]) |> keys |> collect
                end

                # Qflux_finding mode requires certain output
                if my_ev.config[:Qflx_finding] == :on
                    append!(misc_config[rec_key], HOOM.getVariableList(my_mb, :QFLX_FINDING) |> keys )
                end
     
                # Load variables information as a list
                for varname in misc_config[rec_key]

                    println(format("Request output variable: {:s}", varname))
                    if haskey(complete_variable_list, varname)
                        println(format("Using varaible: {:s}", varname))
                        push!(var_list, varname)#, complete_variable_list[varname]... ) )
                    else
                        throw(ErrorException("Unknown varname in " * string(rec_key) * ": " * varname))
                    end
                end

                # additional variables
                add_var_list = []
                #for (k, v) in additional_variable_list
                #   push!(add_var_list, ( k, v... ) )
                #end

                     
                MD.recorders[rec_key] = Recorder(
                    data_table,
                    var_list,
                    HOOM.var_desc;
                    other_varnames=add_var_list,
                )
                
                if length(var_list) != 0
                    push!(activated_record, rec_key) 
                end
            end
 
            # Record (not output) happens AFTER the simulation.
            # Output of the current simulation happens at the
            # BEGINNING of the next simulation.
            #
            # Reason 1:
            # CESM does not simulate the first day of a `continue` run.
            # The first day has been simulated which is the last day of
            # the last run which is exactly the restart file. This is 
            # also why we have to call archive_record! function in the 
            # end of initialization.
            #
            # Reason 2:
            # Output happens at the beginning the next simulation. By
            # doing this we can get rid of the problem of deciding which
            # day is the end of month.
            #
            # This is also the way CAM chooses to do detect the end of
            # current month. 
            # See: http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesmBbrowser/html_code/cam/time_manager.F90.html
            #      is_end_curr_month
            #

            # Must create the record file first because the
            # run of the first day is not called in CESM
            if misc_config[:enable_archive]

                if :daily_record in activated_record
                    recorder_day = MD.recorders[:daily_record]
                    addAlarm!(
                        clock,
                        "[Daily] Create daily output file.",
                        clock.time,
                        2;
                        callback = function (clk, alm)
                            createRecordFile!(MD, "h1.day", recorder_day)
                        end,
                        recurring = Month(1),
                    )


                    addAlarm!(
                        clock,
                        "[Daily] Daily output",
                        clock.time,
                        1;
                        callback = function (clk, alm)
                            record!(recorder_day)
                            avgAndOutput!(recorder_day) # This is important
                        end,
                        recurring = Day(1),
                    )


                end

                if :monthly_record in activated_record

                    # Design alarm such that
                    # (1) Create output file first
                    # (2) Record the initial condition
                    # (3) Record simulation after one day is stepped
                    # (4) If it is the first day of next month
                    #     (i)   avg and output data
                    #     (ii)  create next monthly file
                    #     (iii) record this step in the new file in (ii)
                    #     

                    recorder_mon = MD.recorders[:monthly_record]

                    addAlarm!(
                        clock,
                        "[Monthly] Create monthly output file.",
                        clock.time, # Rings immediately
                        2;
                        callback = function (clk, alm)
                            createRecordFile!(MD, "h0.mon", recorder_mon)
                        end,
                        recurring = Month(1),
                    )
                    
                    addAlarm!(
                        clock,
                        "[Monthly] Daily accumulation using record!",
                        clock.time, # Remember we need to record the first one. So alarm rings immediately.
                        1;
                        callback = function (clk, alm)
                            record!(recorder_mon)
                        end,
                        recurring = Day(1),
                    )

                    addAlarm!(
                        clock,
                        "[Monthly] Average and output monthly data.",
                        clock.time + Month(1), # Start from next month
                        3;  # Higher priority so it outputs data before creating next new monthly file
                        callback = function (clk, alm)
                            avgAndOutput!(recorder_mon)
                        end,
                        recurring = Month(1),
                    )
     
                end

            end
            
        end

        MPI.Barrier(comm)

        syncField!(
            MD.sync_data[:state],
            MD.jdi,
            :M2S,
            :BLOCK,
        ) 
        
        return MD

    end

    function run!(
        MD            :: HOOM_DATA;
        Δt            :: Second,
        write_restart :: Bool,
    )

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        is_master = rank == 0

        syncField!(
            MD.sync_data[:forcing],
            MD.jdi,
            :M2S,
            :BLOCK,
        )

        syncField!(
            MD.sync_data[:bnd_state],
            MD.jdi,
            :M2S,
            :BND,
        ) 

        Δt_float = Float64(Δt.value)
        
        if ! is_master
            #MD.mb.fi.sv[:UVEL] .= (2π / 86400 / 20) * MD.mb.co.gd.R * cos.(MD.mb.co.gd.ϕ_T)
            HOOM.checkBudget!(MD.mb, Δt_float; stage=:BEFORE_STEPPING)
            HOOM.setupForcing!(MD.mb)
        end
        
        substeps = 1
        for substep = 1:substeps
        
            if ! is_master
                HOOM.stepAdvection!(MD.mb, Δt_float/substeps)
                HOOM.checkBudget!(MD.mb, Δt_float; stage=:SUBSTEPS, substeps=substeps)
            end

            syncField!(
                MD.sync_data[:bnd_state],
                MD.jdi,
                :S2M,
                :BND,
            )

            syncField!(
                MD.sync_data[:bnd_state],
                MD.jdi,
                :M2S,
                :BND,
            )
        end

        if ! is_master
            HOOM.stepColumn!(MD.mb, Δt_float)
            HOOM.checkBudget!(MD.mb, Δt_float; stage=:AFTER_STEPPING)
        end
        
        #if write_restart 
        #    writeRestart(MD)
        #end
        
        syncField!(
            MD.sync_data[:state],
            MD.jdi,
            :S2M,
            :BLOCK,
        ) 
    end

    function final(MD::HOOM_DATA)
        
    end


    function createRecordFile!(
        MD     :: HOOM_DATA, 
        group  :: String,
        recorder :: Recorder,
    )

        t = dt2tuple(MD.clock.time)

        filename = format("{}.HOOM.{}.{:04d}-{:02d}.nc", MD.casename, group, t[1], t[2])

        setNewNCFile!(
            recorder,
            joinpath(MD.config[:DRIVER][:caserun], filename)
        )
            
        appendLine(MD.config[:DRIVER][:archive_list], 
            format("mv,{:s},{:s},{:s}",
                filename,
                MD.config[:DRIVER][:caserun],
                joinpath(MD.config[:DRIVER][:archive_root], "ocn", "hist"),
            )
        )

    end


    function output!(
        
        MD               :: HOOM_DATA;
        force_output     :: Bool = false,
    )

        # Daily record block
        if (length(MD.config[:daily_record]) != 0) && t_flags[:new_day]

            #println("##### Avg and Output DAILY!")
            avgAndOutput!(MD.recorders[:daily_record])
        
        end
 
 
        if (length(MD.config[:monthly_record]) != 0) && t_flags[:new_month]
            
            #println("##### Avg and Output MONTHLY!")
            avgAndOutput!(MD.recorders[:monthly_record])

        end
             
    end

    function writeRestart(
        MD :: HOOM_DATA,
    )

        t = MD.timeinfo.t

        restart_file = format("{}.ocn.r.{:04d}{:02d}{:02d}_{:05d}.nc", MD.config[:casename], t[1], t[2], t[3], t[4])
        HOOM.takeSnapshot(MD.ocn, restart_file)
         
        open(MD.config[:rpointer_file], "w") do file
            write(file, restart_file)
        end

        println("(Over)write restart pointer file: ", MD.config[:rpointer_file])
        println("Output restart file: ", restart_file)

        src_dir = MD.config[:caserun]
        dst_dir = joinpath(MD.config[:archive_root], "rest", format("{:04d}-{:02d}-{:02d}-{:05d}", t[1], t[2], t[3], t[4]))

        appendLine(MD.config[:archive_list], 
            format("cp,{:s},{:s},{:s}",
                restart_file,
                src_dir,
                dst_dir,
            )
        )

        appendLine(MD.config[:archive_list], 
            format("cp,{:s},{:s},{:s}",
                MD.config[:rpointer_file],
                src_dir,
                dst_dir,
            )
        )

    end
end
