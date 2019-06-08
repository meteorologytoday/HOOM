include(joinpath(@__DIR__, "..", "..", "..", "SOM", "SOM.jl"))
module CESMCORE_SOM

    include("Workspace_SOM.jl")
    include(joinpath(@__DIR__, "..", "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "AppendLine.jl"))

    using Formatting
    using ..NetCDFIO
    using ..SOM
    using NCDatasets
    using .RecordTool

    name = "SOM"
    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct SOM_DATA
        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        occ         :: SOM.OceanColumnCollection

        x2o         :: Dict
        o2x         :: Dict

        configs     :: Dict

        output_vars :: Dict
        wksp        :: Workspace

        recorders   :: Dict
    end


    function init(;
        casename     :: AbstractString,
        map          :: NetCDFIO.MapInfo,
        init_file    :: Union{Nothing, AbstractString}=nothing,
        t            :: AbstractArray{Integer},
        configs      :: Dict,
        read_restart :: Bool,
    )

        # If `read_restart` is true then read restart file: configs["rpointer_file"]
        # If not then initialize ocean with default profile if `initial_file`
        # is "nothing", with `init_file` if it is nonempty.

        if read_restart

            println("`read_restart` is on")

            if ! ("rpointer_file" in keys(configs))
                throw(ErrorException("Cannot find `rpointer_file` in configs!"))
            end

            if !isfile(configs["rpointer_file"])
                throw(ErrorException(configs["rpointer_file"] * " does not exist!"))
            end
            
            println("Going to read restart pointer file: ", configs["rpointer_file"])

            open(configs["rpointer_file"], "r") do file
                init_file = readline(file)
            end

            if !isfile(init_file)
                throw(ErrorException(init_file * " does not exist!"))
            end
 
        end


        if typeof(init_file) <: AbstractString
            println("Initial ocean with profile: ", init_file)
            occ = SOM.loadSnapshot(init_file)
        else
            println("No initial ocean profile. Using the naive one.")
            occ = SOM.makeBlankOceanColumnCollection(map.nx, map.ny; mask=map.mask)
            occ.Kh_T = 0.0
            occ.h_ML .= 30.0
            occ.T_ML .= 288.0

            snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
            println("Output snapshot: ", snapshot_file)
            SOM.takeSnapshot(occ, joinpath(configs["short_term_archive_dir"], snapshot_file))
            appendLine(configs["short_term_archive_list"], snapshot_file)
        end

        wksp = Workspace(occ.Nx, occ.Ny)

        x2o = Dict(
            "SWFLX"  => wksp.swflx,
            "NSWFLX" => wksp.nswflx,
            "TFDIV"  => wksp.tfdiv,
            "MLD"    => occ.h_ML,
        )

        o2x = Dict(
            "T_ML"      => occ.T_ML,
            "h_ML"      => occ.h_ML,
            "QFLX2ATM"  => occ.qflx2atm,
        )

        output_vars = Dict(#=
            "T_ML"      => occ.T_ML,
            "h_ML"      => occ.h_ML,
            "eflx"      => wksp.eflx,
            "tfdiv"     => wksp.tfdiv,
            =#
        )
 
        recorders = Dict()

        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]
                 recorder = RecordTool.Recorder(
                    Dict(
                        "Nx" => occ.Nx,
                        "Ny" => occ.Ny,
                    ), [
                        ("T",      occ.T_ML,   ("Nx", "Ny")),
                        ("MLD",    occ.h_ML,   ("Nx", "Ny")),
                        ("swflx",  wksp.swflx, ("Nx", "Ny")),
                        ("nswflx", wksp.swflx, ("Nx", "Ny")),
                        ("eflx",   wksp.eflx,  ("Nx", "Ny")),
                    ],
                )

                recorders[rec_key] = recorder
               
            end
        end

        # check if key "Qflux" is set in configs
        if haskey(configs, "Qflux")
            if ! (typeof(configs["Qflux"]) <: Bool)
                throw(ErrorException("The value of `Qflux` in configs must be of type Bool."))
            end
        else
            throw(ErrorException("The key `Qflux` is not found in variable `configs`."))
        end

        return SOM_DATA(
            casename,
            map,
            occ,
            x2o,
            o2x,
            configs,
            output_vars,
            wksp,
            recorders,
        )

    end

    function run(
        MD            :: SOM_DATA;
        t             :: AbstractArray{Integer},
        t_cnt         :: Integer,
        t_flags       :: Dict,
        Δt            :: Float64,
        substep       :: Integer,
        write_restart :: Bool,
    )

        if MD.configs["enable_short_term_archive"]

            if MD.configs["daily_record"]
                
                if t_flags["new_month"] && substep == 1
                        filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
                        RecordTool.setNewNCFile!(
                            MD.recorders["daily_record"],
                            joinpath(MD.configs["short_term_archive_dir"], filename)
                        )
                        appendLine(MD.configs["short_term_archive_list"], filename)
                end

                RecordTool.record!(
                    MD.recorders["daily_record"];
                    avg_and_output = ( t_flags["new_day"] && substep==1 && t_cnt != 1)
                )

            end
            
            if MD.configs["monthly_record"]

                if t_flags["new_year"] && substep == 1
                        filename = format("{}.ocn.h.monthly.{:04d}.nc", MD.casename, t[1])
                        RecordTool.setNewNCFile!(
                            MD.recorders["monthly_record"],
                            joinpath(MD.configs["short_term_archive_dir"], filename)
                        )
                        appendLine(MD.configs["short_term_archive_list"], filename)
                end

                RecordTool.record!(
                    MD.recorders["monthly_record"];
                    avg_and_output = ( t_flags["new_month"] && substep==1 && t_cnt != 1)
                )


            end

        end

        wksp = MD.wksp

        # Only process incoming data for the first time!! 
        if substep == 1

            wksp.nswflx .*= -1.0
            wksp.swflx .*= -1.0

            wksp.eflx[:, :] = wksp.nswflx[:, :]
            wksp.eflx .+= wksp.swflx

            if MD.configs["Qflux"]
                wksp.eflx .+= wksp.tfdiv
            end

        end


        SOM.stepOceanColumnCollection!(
            MD.occ;
            eflx   = wksp.eflx,
            Δt     = Δt,
        )
        
        if write_restart
            restart_file = format("restart.ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
            SOM.takeSnapshot(MD.occ, restart_file)
             
            open(MD.configs["rpointer_file"], "w") do file
                write(file, restart_file)
            end

            println("(Over)write restart pointer file: ", MD.configs["rpointer_file"])
            println("Output restart file: ", restart_file)
        end

    end

    function final(MD::SOM_DATA)
        
    end


end
