include(joinpath(@__DIR__, "..", "..", "..", "NKOM", "NKOM.jl"))
module CESMCORE_NKOM

    include(joinpath(@__DIR__, "..", "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "CheckDict.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "AppendLine.jl"))

    using Formatting
    using ..NetCDFIO
    using ..NKOM
    using NCDatasets
    using .RecordTool

    name = "NKOM"

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct NKOM_DATA
        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        occ         :: NKOM.OceanColumnCollection

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
            (:MLD_scheme,                    true, (:prognostic, :datastream,), nothing),
            (:Qflux_scheme,                  true, (:on, :off,),                nothing),
            (:diffusion_scheme,              true, (:on, :off,),                nothing),
            (:relaxation_scheme,             true, (:on, :off,),                nothing),
            (:convective_adjustment_scheme,  true, (:on, :off,),                nothing),
            (:radiation_scheme,              true, (:exponential, :step,),      nothing),
            (:daily_record,                  true, (Bool,),                     nothing),
            (:monthly_record,                true, (Bool,),                     nothing),
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
            occ = NKOM.loadSnapshot(init_file)
        
        else

            throw(ErrorException("Variable `init_file` is absent in `configs`."))

        end

        println("Initializing parallization...")
        NKOM.init(occ)

        in_flds = occ.in_flds

        #
        # If it is "datastream", entrainment speed w_e would be 
        # calculated from h given. In fact there is no need
        # to calculate w_e.
        #
        # If it is "prognostic", entrainment speed w_e would be
        # calculated accroding to Niiler and Kraus dynamics.
        #

        if configs[:MLD_scheme] == :datastream

            x2o = Dict(
                "SWFLX"  => in_flds.swflx,
                "NSWFLX" => in_flds.nswflx,
                "FRWFLX" => in_flds.frwflx,
                "TFDIV"  => in_flds.qflx,
                "MLD"    => in_flds.h_ML,
            )

        elseif configs[:MLD_scheme] == :prognostic

            x2o = Dict(
                "SWFLX"  => in_flds.swflx,
                "NSWFLX" => in_flds.nswflx,
                "TAUX"  => in_flds.taux,
                "TAUY"  => in_flds.tauy,
                "IFRAC" => in_flds.ifrac,
                "FRWFLX" => in_flds.frwflx,
                "TFDIV"  => in_flds.qflx,
            )

        end

        o2x = Dict(
            "SST"      => occ.T_ML,
            "QFLX2ATM" => occ.qflx2atm,
        )
        
        recorders = Dict()

        for rec_key in [:daily_record, :monthly_record]
            if configs[rec_key]
                 recorder = RecordTool.Recorder(
                    Dict(
                        "Nx" => occ.Nx,
                        "Ny" => occ.Ny,
                        "Nz_bone" => occ.Nz_bone,
                    ), [
                        ("T",     NKOM.toXYZ(occ.Ts, :zxy), ("Nx", "Ny", "Nz_bone")),
                        ("S",     NKOM.toXYZ(occ.Ss, :zxy), ("Nx", "Ny", "Nz_bone")),
                        ("T_ML",  occ.T_ML, ("Nx", "Ny",)),
                        ("S_ML",  occ.S_ML, ("Nx", "Ny",)),
                        ("h_ML",  occ.h_ML, ("Nx", "Ny")),
                        ("weighted_fric_u",  occ.in_flds.weighted_fric_u, ("Nx", "Ny")),
                        ("fric_u",  occ.in_flds.fric_u, ("Nx", "Ny")),
                    ],
                )

                recorders[rec_key] = recorder
               
            end
        end

        return NKOM_DATA(
            casename,
            map,
            occ,
            x2o,
            o2x,
            configs,
            recorders,
        )

    end

    function run(
        MD            :: NKOM_DATA;
        t             :: AbstractArray{Integer},
        t_cnt         :: Integer,
        t_flags       :: Dict,
        Δt            :: Float64,
        substep       :: Integer,
        write_restart :: Bool,
    )

        if MD.configs[:enable_short_term_archive] && substep == 1

            if MD.configs[:daily_record]
 
                RecordTool.record!(
                    MD.recorders[:daily_record];
                    avg_and_output = ( t_flags[:new_day] && t_cnt != 1)
                )
               
                if t_flags[:new_month]
                        filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
                        RecordTool.setNewNCFile!(
                            MD.recorders[:daily_record],
                            joinpath(MD.configs[:short_term_archive_dir], filename)
                        )
                        appendLine(MD.configs[:short_term_archive_list], filename)
                end


            end
            
            if MD.configs[:monthly_record]

                RecordTool.record!(
                    MD.recorders[:monthly_record];
                    avg_and_output = ( t_flags[:new_month] && t_cnt != 1)
                )

                if t_flags[:new_year]
                        filename = format("{}.ocn.h.monthly.{:04d}.nc", MD.casename, t[1])
                        RecordTool.setNewNCFile!(
                            MD.recorders[:monthly_record],
                            joinpath(MD.configs[:short_term_archive_dir], filename)
                        )
                        appendLine(MD.configs[:short_term_archive_list], filename)
                end


            end

        end

        in_flds = MD.occ.in_flds

        # Only process incoming data for the first time!! 
        if substep == 1

            in_flds.nswflx .*= -1.0
            in_flds.swflx  .*= -1.0
            #in_flds.nswflx  .*= 0.0
            #in_flds.swflx  .=-1000.0
            #if MD.configs[:MLD_scheme] == :prognostic
            #in_flds.fric_u .= sqrt.(sqrt.((in_flds.taux).^2.0 .+ (in_flds.tauy).^2.0) / NKOM.ρ)
            #in_flds.weighted_fric_u .*= (1.0 .- in_flds.ifrac)
            #end

        end
       
        NKOM.run!(
            MD.occ;
            use_qflx      = MD.configs[:Qflux_scheme] == :on,
            use_h_ML      = MD.configs[:MLD_scheme] == :datastream,
            Δt            = Δt,
            diffusion_Δt  = Δt * MD.configs[:substeps],
            relaxation_Δt = Δt * MD.configs[:substeps],
            do_diffusion  = (MD.configs[:diffusion_scheme] == :on && substep == MD.configs[:substeps]),
            do_relaxation = (MD.configs[:relaxation_scheme] == :on && substep == MD.configs[:substeps]),
            do_convadjust = MD.configs[:convective_adjustment_scheme] == :on,
            rad_scheme    = MD.configs[:radiation_scheme],
            copy_in_flds  = substep == 1,
        )

        # Force workers to update master's profile
        if substep == MD.configs[:substeps]
            NKOM.sync!(MD.occ)
        end
        
        if write_restart
            restart_file = format("restart.ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
            NKOM.takeSnapshot(MD.occ, restart_file)
             
            open(MD.configs[:rpointer_file], "w") do file
                write(file, restart_file)
            end

            println("(Over)write restart pointer file: ", MD.configs[:rpointer_file])
            println("Output restart file: ", restart_file)
        end

    end

    function final(MD::NKOM_DATA)
        
    end

end
