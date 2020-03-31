include(joinpath(@__DIR__, "HOOM.jl"))
module CESMCORE_HOOM

    include(joinpath(@__DIR__, "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "CheckDict.jl"))
    include(joinpath(@__DIR__, "..", "..", "share", "AppendLine.jl"))

    using Formatting
    using ..NetCDFIO
    using ..HOOM
    using NCDatasets
    using .RecordTool

    name = "HOOM"

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct HOOM_DATA
        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        ocn         :: HOOM.Ocean

        x2o         :: Dict
        o2x         :: Dict

        configs       :: Dict

        recorders   :: Dict

    end


    function init(;
        casename     :: AbstractString,
        map          :: NetCDFIO.MapInfo,
        t            :: AbstractArray{Integer},
        configs      :: Dict,
        read_restart :: Bool,
    )

        checkDict!(configs, [
            (:init_file,                    false, (nothing, String,),          nothing),
            (:advection_scheme,              true, (:static, :ekman_all_in_ML, :ekman_simple_partition, :ekman_codron2012_partition,),  nothing),
            (:MLD_scheme,                    true, (:prognostic, :datastream,), nothing),
            (:Qflux_scheme,                  true, (:on, :off,),                nothing),
            (:Qflux_finding,                 true, (:on, :off,),                nothing),
            (:vertical_diffusion_scheme,     true, (:on, :off,),                nothing),
            (:horizontal_diffusion_scheme,   true, (:on, :off,),                nothing),
            (:relaxation_scheme,             true, (:on, :off,),                nothing),
            (:convective_adjustment_scheme,  true, (:on, :off,),                nothing),
            (:radiation_scheme,              true, (:exponential, :step,),      nothing),
            (:daily_record,                  true, (AbstractArray, Symbol),                 []),
            (:monthly_record,                true, (AbstractArray, Symbol),                 []),
            (:turn_off_frwflx,              false, (Bool,),                       false),
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
            ocn = HOOM.loadSnapshot(init_file)
        
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
            "T_CLIM"  => in_flds.Tclim,
            "S_CLIM"  => in_flds.Sclim,
            "MLD"    => in_flds.h_ML,
        )


        o2x = Dict(
            "SST"      => ocn.T_ML,
            "QFLX2ATM" => ocn.qflx2atm,
        )
        
        recorders = Dict()
        complete_variable_list   = HOOM.getCompleteVariableList(ocn)
        additional_variable_list = HOOM.getAdditionalVariableList(ocn)

        for rec_key in [:daily_record, :monthly_record]
    
            println("# For record key: " * string(rec_key))

            var_list = []
           
            if configs[rec_key] == :ALL
                configs[rec_key] = keys(complete_variable_list)
            elseif configs[rec_key] == :ESSENTIAL
                configs[rec_key] = keys(HOOM.getSubVariableList(ocn, :ESSENTIAL))
            else
                println("Using customized output variables.")
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
                    "zs_bone" => length(ocn.zs_bone),
                ),
                var_list,
                HOOM.var_desc;
                other_vars=add_var_list
            )
               
        end

        MD = HOOM_DATA(
            casename,
            map,
            ocn,
            x2o,
            o2x,
            configs,
            recorders,
        )

        # Must record the first day
        record!(
            MD;
            t                = t,
            called_from_init = true,
        ) 


        return MD

    end

    function run(
        MD            :: HOOM_DATA;
        t             :: AbstractArray{Integer},
        t_cnt         :: Integer,
        t_flags       :: Dict,
        Δt            :: Float64,
        write_restart :: Bool,
    )

        # process input fields before record
        in_flds = MD.ocn.in_flds

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
        )


        record!(
            MD;
            t       = t,
            t_cnt   = t_cnt,
            t_flags = t_flags,
        ) 

        
        if write_restart

            restart_file = format("restart.ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
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

    function final(MD::HOOM_DATA)
        
    end




    # Record happens at the end of day.
    # This must be that way because CESM does not
    # simulate the first day of a `continue` run.
    # The first day has been simulated which is the
    # last day of the last run. So that is exactly
    # the restart file.
    function record!(
        MD               :: HOOM_DATA;
        t                :: AbstractArray{Integer},
        t_cnt            :: Union{Integer, Nothing} = nothing,
        t_flags          :: Union{Dict, Nothing}    = nothing,
        called_from_init :: Bool                    = false,
    )
        if MD.configs[:enable_archive]

            substeps = MD.configs[:substeps]


            # Daily record block
            if length(MD.configs[:daily_record]) != 0

                if called_from_init || t_flags[:new_month]

                        filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
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

                 RecordTool.record!(
                    MD.recorders[:daily_record];
                    avg_and_output =  called_from_init || ( t_flags[:new_day] && t_cnt == substeps )
                 )
            
            end
            
            # Monthly record block
            if length(MD.configs[:monthly_record]) != 0


                if called_from_init || t_flags[:new_month]

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

                RecordTool.record!(
                    MD.recorders[:monthly_record];
                    avg_and_output = called_from_init || ( t_flags[:new_month] && t_cnt == substeps)
                )

            end

        end

    end


end
