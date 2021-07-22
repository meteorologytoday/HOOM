
if ! ( :HOOM in names(Main) )
    include(joinpath(@__DIR__, "HOOM.jl"))
end

module CESMCORE_HOOM

    include(joinpath(@__DIR__, "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "CheckDict.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "AppendLine.jl"))

    using Dates
    using Formatting
    using ..HOOM
    using ..ModelClockSystem
    using NCDatasets
    using .RecordTool


    name = "HOOM"

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct HOOM_DATA
        casename    :: AbstractString
        ocn         :: HOOM.Ocean
        clock       :: ModelClock

        x2o         :: Dict
        o2x         :: Dict

        configs       :: Dict

        recorders   :: Dict
    end


    function init(;
        casename     :: AbstractString,
        clock        :: ModelClock,
        configs      :: Dict,
        read_restart :: Bool,
    )

        checkDict!(configs, [
            (:init_file,                    false, (nothing, String,),          nothing),
            (:advection_scheme,              true, (:static, :ekman_all_in_ML, :ekman_simple_partition, :ekman_codron2012_partition,),  nothing),
            (:MLD_scheme,                    true, (:prognostic, :datastream,), nothing),
            (:Qflux_scheme,                  true, (:on, :off,),                nothing),
            (:Qflux_finding,                 true, (:on, :off,),                nothing),
            (:seaice_nudging,                 true, (:on, :off,),                nothing),
            (:vertical_diffusion_scheme,     true, (:on, :off,),                nothing),
            (:horizontal_diffusion_scheme,   true, (:on, :off,),                nothing),
            (:relaxation_scheme,             true, (:on, :off,),                nothing),
            (:convective_adjustment_scheme,  true, (:on, :off,),                nothing),
            (:radiation_scheme,              true, (:exponential, :step,),      nothing),
            (:daily_record,                  true, (AbstractArray, Symbol),                 []),
            (:monthly_record,                true, (AbstractArray, Symbol),                 []),
        ])


        if configs[:MLD_scheme] == :prognostic && configs[:radiation_scheme] == :step
            throw(ErrorException(":prognostic MLD_scheme cannot use :step in radiation_scheme."))
        end


        init_file = configs[:init_file]

        # If `read_restart` is true then read restart file: configs[:rpointer_file]
        # If not then initialize ocean with default profile if `initial_file`
        # is "nothing", with `init_file` if it is nonempty.

        if read_restart

            println("`read_restart` is on")
            checkDict!( configs, [
                (:rpointer_file, true, (String,), nothing),
            ])

            if !isfile(configs[:rpointer_file])
                throw(ErrorException(configs[:rpointer_file] * " does not exist!"))
            end
            
            println("Going to read restart pointer file: ", configs[:rpointer_file])

            open(configs[:rpointer_file], "r") do file
                init_file = readline(file)
            end

            if !isfile(init_file)
                throw(ErrorException(format("Initial file \"{:s}\" does not exist!", init_file)))
            end
 
        end


        if typeof(init_file) <: AbstractString

            println("Initial ocean with profile: ", init_file)
            println("Initial ocean with domain file: ", configs[:domain_file])
            ocn = HOOM.loadSnapshot(init_file; gridinfo_file=configs[:domain_file])
        
        else

            throw(ErrorException("Variable `init_file` is absent in `configs`."))

        end

        println("Initializing parallization...")
        HOOM.init(ocn)

        in_flds = ocn.in_flds
        #
        # If it is "datastream", entrainment speed w_e would be 
        # calculated from h given. In fact there is no need
        # to calculate w_e.
        #
        # If it is "prognostic", entrainment speed w_e would be
        # calculated accroding to Niiler-Kraus dynamics.
        #

        x2o = Dict(
            "SWFLX"  => in_flds.swflx,
            "NSWFLX" => in_flds.nswflx,
            "TAUX"   => in_flds.taux,
            "TAUY"   => in_flds.tauy,
            "IFRAC"  => in_flds.ifrac,
            "FRWFLX" => in_flds.frwflx,
            "VSFLX"  => in_flds.vsflx,
            "QFLX_T" => in_flds.qflx_T,
            "QFLX_S" => in_flds.qflx_S,
            "T_CLIM"      => in_flds.Tclim,
            "S_CLIM"      => in_flds.Sclim,
            "IFRAC_CLIM"  => in_flds.IFRACclim,
            "MLD"    => in_flds.h_ML,
        )


        o2x = Dict(
            "SST"      => ocn.T_ML,
            "QFLX2ATM" => ocn.qflx2atm,
        )
        
        recorders = Dict()
        complete_variable_list   = HOOM.getVariableList(ocn, :ALL)
        additional_variable_list = HOOM.getVariableList(ocn, :COORDINATE)

        for rec_key in [:daily_record, :monthly_record]
    
            println("# For record key: " * string(rec_key))

            var_list = []
            
            if typeof(configs[rec_key]) <: Symbol
                configs[rec_key] = HOOM.getVariableList(ocn, configs[rec_key]) |> keys |> collect
            else
                println("Using customized output variables.")
            end

            # Qflux_finding mode requires certain output
            #if (length(configs[rec_key]) != 0 ) && (configs[:Qflux_finding] == :on)
            if configs[:Qflux_finding] == :on
                append!(configs[rec_key], HOOM.getVariableList(ocn, :QFLX_FINDING) |> keys )
            end
 
            # Load variables information as a list
            for varname in configs[rec_key]

                println(format("Request output variable: {:s}", varname))
                if haskey(complete_variable_list, varname)
                    println(format("Using varaible: {:s}", varname))
                    push!(var_list, ( varname, complete_variable_list[varname]... ) )
                else
                    throw(ErrorException("Unknown varname in " * string(rec_key) * ": " * varname))
                end
            end

            # additional variables
            add_var_list = []
            for (k, v) in additional_variable_list
               push!(add_var_list, ( k, v... ) )
            end

 
            recorders[rec_key] = RecordTool.Recorder(
                Dict(
                    "Nx" => ocn.Nx,
                    "Ny" => ocn.Ny,
                    "Nz_bone" => ocn.Nz_bone,
                    "NP_zs_bone" => length(ocn.zs_bone),
                ),
                var_list,
                HOOM.var_desc;
                other_vars=add_var_list
            )
               
        end

        MD = HOOM_DATA(
            casename,
            ocn,
            clock,
            x2o,
            o2x,
            configs,
            recorders,
        )


        # Must create the record file first because the
        # run of the first day is not called in CESM
        if MD.configs[:enable_archive]

            if length(MD.configs[:daily_record]) != 0
                recorder_day = recorders[:daily_record]
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

            if length(MD.configs[:monthly_record]) != 0

                # Design alarm such that
                # (1) Create output file first
                # (2) Record the initial condition
                # (3) Record simulation after one day is stepped
                # (4) If it is the first day of next month
                #     (i)   avg and output data
                #     (ii)  create next monthly file
                #     (iii) record this step in the new file in (ii)
                #     

                recorder_mon = recorders[:monthly_record]

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





#=
        run!(
            MD;
            Δt            = 0.0,
            write_restart = false,
            first_run     = true,
        )
=#

        return MD

    end

    function run!(
        MD            :: HOOM_DATA;
        Δt            :: Float64,
        write_restart :: Bool,
        first_run     :: Bool = false,
    )

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
#            archive_outputIfNeeded!(MD)

            # File must be created AFTER it is output.
#            archive_createFileIfNeeded!(MD)

        HOOM.run!(
            MD.ocn;
            substeps         = MD.configs[:substeps],
            use_h_ML         = MD.configs[:MLD_scheme] == :datastream,
            Δt               = Δt,
            do_vert_diff     = MD.configs[:vertical_diffusion_scheme] == :on,
            do_horz_diff     = MD.configs[:horizontal_diffusion_scheme] == :on,
            do_relaxation    = MD.configs[:relaxation_scheme] == :on,
            do_convadjust    = MD.configs[:convective_adjustment_scheme] == :on,
            rad_scheme       = MD.configs[:radiation_scheme],
            adv_scheme       = MD.configs[:advection_scheme],
            do_qflx          = MD.configs[:Qflux_scheme] == :on,
            do_qflx_finding  = MD.configs[:Qflux_finding] == :on,
            do_seaice_nudging = MD.configs[:seaice_nudging] == :on,
        )

        archive_record!(MD)
        
        if write_restart #|| MD.timeinfo.t_flags[:new_month]
            writeRestart(MD)
        end

    end

    function final(MD::HOOM_DATA)
        
    end


    function createRecordFile!(
        MD     :: HOOM_DATA, 
        group  :: String,
        recorder :: RecordTool.Recorder,
    )

        t = dt2tuple(MD.clock.time)

        filename = format("{}.HOOM.{}.{:04d}-{:02d}.nc", MD.casename, group, t[1], t[2])

        setNewNCFile!(
            recorder,
            joinpath(MD.configs[:caserun], filename)
        )
            
        appendLine(MD.configs[:archive_list], 
            format("mv,{:s},{:s},{:s}",
                filename,
                MD.configs[:caserun],
                joinpath(MD.configs[:archive_root], "ocn", "hist"),
            )
        )

    end

#=
    function archive_createFileIfNeeded!(
        MD :: HOOM_DATA;
    )

        if ! MD.configs[:enable_archive]
            return
        end

        t_flags = MD.timeinfo.t_flags
        t = MD.timeinfo.t
 
        if (length(MD.configs[:daily_record]) != 0) && t_flags[:new_month]
            filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
            
            #println("CREATE NEW DAILY FILE:", filename)
            RecordTool.setNewNCFile!(
                MD.recorders[:daily_record],
                joinpath(MD.configs[:caserun], filename)
            )
            
            appendLine(MD.configs[:archive_list], 
                format("mv,{:s},{:s},{:s}",
                    filename,
                    MD.configs[:caserun],
                    joinpath(MD.configs[:archive_root], "ocn", "hist"),
                )
            )


        end
 
        # Monthly record block
        if (length(MD.configs[:monthly_record]) != 0) && t_flags[:new_month]

            filename = format("{}.ocn.h.monthly.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
            RecordTool.setNewNCFile!(
                MD.recorders[:monthly_record],
                joinpath(MD.configs[:caserun], filename)
            )

            appendLine(MD.configs[:archive_list], 
                format("mv,{:s},{:s},{:s}",
                    filename,
                    MD.configs[:caserun],
                    joinpath(MD.configs[:archive_root], "ocn", "hist"),
                )
            )

        end

           
    end
=#


    function output!(
        
        MD               :: HOOM_DATA;
        force_output     :: Bool = false,
    )

        # Daily record block
        if (length(MD.configs[:daily_record]) != 0) && t_flags[:new_day]

            #println("##### Avg and Output DAILY!")
            RecordTool.avgAndOutput!(MD.recorders[:daily_record])
        
        end
 
 
        if (length(MD.configs[:monthly_record]) != 0) && t_flags[:new_month]
            
            #println("##### Avg and Output MONTHLY!")
            RecordTool.avgAndOutput!(MD.recorders[:monthly_record])

        end
             
    end

    function writeRestart(
        MD :: HOOM_DATA,
    )

        t = MD.timeinfo.t

        restart_file = format("{}.ocn.r.{:04d}{:02d}{:02d}_{:05d}.nc", MD.configs[:casename], t[1], t[2], t[3], t[4])
        HOOM.takeSnapshot(MD.ocn, restart_file)
         
        open(MD.configs[:rpointer_file], "w") do file
            write(file, restart_file)
        end

        println("(Over)write restart pointer file: ", MD.configs[:rpointer_file])
        println("Output restart file: ", restart_file)

        src_dir = MD.configs[:caserun]
        dst_dir = joinpath(MD.configs[:archive_root], "rest", format("{:04d}-{:02d}-{:02d}-{:05d}", t[1], t[2], t[3], t[4]))

        appendLine(MD.configs[:archive_list], 
            format("cp,{:s},{:s},{:s}",
                restart_file,
                src_dir,
                dst_dir,
            )
        )

        appendLine(MD.configs[:archive_list], 
            format("cp,{:s},{:s},{:s}",
                MD.configs[:rpointer_file],
                src_dir,
                dst_dir,
            )
        )

    end
end
