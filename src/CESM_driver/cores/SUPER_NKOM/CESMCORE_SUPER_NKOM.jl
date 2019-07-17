include(joinpath(@__DIR__, "..", "..", "..", "SUPER_NKOM", "SUPER_NKOM.jl"))
module CESMCORE_SUPER_NKOM

    include("Workspace_SUPER_NKOM.jl")
    include(joinpath(@__DIR__, "..", "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "CheckDict.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "AppendLine.jl"))

    using Formatting
    using ..NetCDFIO
    using ..SUPER_NKOM
    using NCDatasets
    using .RecordTool

    name = "SUPER_NKOM"

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct SUPER_NKOM_DATA
        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        occ         :: SUPER_NKOM.OceanColumnCollection

        x2o         :: Dict
        o2x         :: Dict

        configs       :: Dict

        wksp        :: Workspace

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
            ("init_file",                    false, (nothing, String,),             nothing),
            ("MLD_scheme",                    true, ("prognostic", "datastream",),  nothing),
            ("Qflux_scheme",                  true, ("on", "off",),                 nothing),
            ("diffusion_scheme",              true, ("on", "off",),                 nothing),
            ("relaxation_scheme",             true, ("on", "off",),                 nothing),
            ("convective_adjustment_scheme",  true, ("on", "off",),                 nothing),
            ("daily_record",                  true, (Bool,),                        nothing),
            ("monthly_record",                true, (Bool,),                        nothing),
        ])

        init_file = configs["init_file"]

        # If `read_restart` is true then read restart file: configs["rpointer_file"]
        # If not then initialize ocean with default profile if `initial_file`
        # is "nothing", with `init_file` if it is nonempty.

        if read_restart

            println("`read_restart` is on")
            checkDict!( configs, [
                ("rpointer_file", true, (String,), nothing),
            ])

            if !isfile(configs["rpointer_file"])
                throw(ErrorException(configs["rpointer_file"] * " does not exist!"))
            end
            
            println("Going to read restart pointer file: ", configs["rpointer_file"])

            open(configs["rpointer_file"], "r") do file
                init_file = readline(file)
            end

            if !isfile(init_file)
                throw(ErrorException(format("Initial file \"{:s}\" does not exist!", init_file)))
            end
 
        end


        if typeof(init_file) <: AbstractString

            println("Initial ocean with profile: ", init_file)
            occ = SUPER_NKOM.loadSnapshot(init_file)
        
        else

            throw(ErrorException("Variable `init_file` is absent in `configs`."))

        end

        wksp = Workspace(occ.Nx, occ.Ny, occ.Nz_bone)

        #
        # If it is "datastream", entrainment speed w_e would be 
        # calculated from h given. In fact there is no need
        # to calculate w_e.
        #
        # If it is "prognostic", entrainment speed w_e would be
        # calculated accroding to Niiler and Kraus dynamics.
        #

        if configs["MLD_scheme"] == "datastream"

            wksp.h_ML = nothing
            x2o = Dict(
                "SWFLX"  => wksp.swflx,
                "NSWFLX" => wksp.nswflx,
                "FRWFLX" => wksp.frwflx,
                "TFDIV"  => wksp.qflx,
                "MLD"    => wksp.h_ML,
            )

        elseif configs["MLD_scheme"] == "prognostic"

            x2o = Dict(
                "SWFLX"  => wksp.swflx,
                "NSWFLX" => wksp.nswflx,
                "TAUX"  => wksp.taux,
                "TAUY"  => wksp.tauy,
                "IFRAC" => wksp.ifrac,
                "FRWFLX" => wksp.frwflx,
                "TFDIV"  => wksp.qflux,
            )

        end


        if configs["Qflux_scheme"] == "off"
            wksp.qflx = nothing
        else
            ErrorException("Unknown Qflux_scheme: ", configs["Qflux_scheme"]) |> throw
        end

        o2x = Dict(
            "SST"      => occ.T_ML,
            "QFLX2ATM" => occ.qflx2atm,
        )
        
        recorders = Dict()

        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]
                 recorder = RecordTool.Recorder(
                    Dict(
                        "Nx" => occ.Nx,
                        "Ny" => occ.Ny,
                        "Nz_bone" => occ.Nz_bone,
                    ), [
                        ("T",     occ.Ts, ("Nx", "Ny", "Nz_bone")),
                        ("S",     occ.Ss, ("Nx", "Ny", "Nz_bone")),
                        ("MLD",   occ.h_ML, ("Nx", "Ny")),
                    ],
                )

                recorders[rec_key] = recorder
               
            end
        end

        return SUPER_NKOM_DATA(
            casename,
            map,
            occ,
            x2o,
            o2x,
            configs,
            wksp,
            recorders,
        )

    end

    function run(
        MD            :: SUPER_NKOM_DATA;
        t             :: AbstractArray{Integer},
        t_cnt         :: Integer,
        t_flags       :: Dict,
        Δt            :: Float64,
        substep       :: Integer,
        write_restart :: Bool,
    )

        if MD.configs["enable_short_term_archive"]

            if MD.configs["daily_record"]
 
                RecordTool.record!(
                    MD.recorders["daily_record"];
                    avg_and_output = ( t_flags["new_day"] && substep==1 && t_cnt != 1)
                )
               
                if t_flags["new_month"] && substep == 1
                        filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
                        RecordTool.setNewNCFile!(
                            MD.recorders["daily_record"],
                            joinpath(MD.configs["short_term_archive_dir"], filename)
                        )
                        appendLine(MD.configs["short_term_archive_list"], filename)
                end


            end
            
            if MD.configs["monthly_record"]

                RecordTool.record!(
                    MD.recorders["monthly_record"];
                    avg_and_output = ( t_flags["new_month"] && substep==1 && t_cnt != 1)
                )

                if t_flags["new_year"] && substep == 1
                        filename = format("{}.ocn.h.monthly.{:04d}.nc", MD.casename, t[1])
                        RecordTool.setNewNCFile!(
                            MD.recorders["monthly_record"],
                            joinpath(MD.configs["short_term_archive_dir"], filename)
                        )
                        appendLine(MD.configs["short_term_archive_list"], filename)
                end


            end

        end

        wksp = MD.wksp

        # Only process incoming data for the first time!! 
        if substep == 1

            wksp.nswflx .*= -1.0
            wksp.swflx  .*= -1.0

            if configs["MLD_scheme"] == "prognostic"
                wksp.fric_u .= sqrt.(sqrt.((wksp.taux).^2.0 .+ (wksp.tauy).^2.0) / SUPER_NKOM.ρ)
                wksp.weighted_fric_u .*= (1.0 .- wksp.ifrac)
            end

        end



        SUPER_NKOM.stepOceanColumnCollection!(
            MD.occ;
            fric_u = wksp.weighted_fric_u,
            swflx  = wksp.swflx,
            nswflx = wksp.nswflx,
            frwflx = wksp.frwflx,
            qflx   = wksp.qflx,
            h_ML   = wksp.h_ML,
            Δt     = Δt,
            diffusion_Δt  = Δt * MD.configs["substeps"],
            relaxation_Δt = Δt * MD.configs["substeps"],
            do_diffusion  = (MD.configs["diffusion_scheme"] == "on" && substep == MD.configs["substeps"]),
            do_relaxation = (MD.configs["relaxation_scheme"] == "on" && substep == MD.configs["substeps"]),
            do_convadjust = MD.configs["convective_adjustment_scheme"] == "on",
        )

        if write_restart
            restart_file = format("restart.ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
            SUPER_NKOM.takeSnapshot(MD.occ, restart_file)
             
            open(MD.configs["rpointer_file"], "w") do file
                write(file, restart_file)
            end

            println("(Over)write restart pointer file: ", MD.configs["rpointer_file"])
            println("Output restart file: ", restart_file)
        end

    end

    function final(MD::SUPER_NKOM_DATA)
        
    end

end