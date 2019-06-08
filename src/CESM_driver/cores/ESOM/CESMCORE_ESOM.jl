include(joinpath(@__DIR__, "..", "..", "..", "ESOM", "ESOM.jl"))
module CESMCORE_ESOM

    include(joinpath(@__DIR__, "Workspace.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "RecordTool.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "AppendLine.jl"))

    using Formatting
    using ..NetCDFIO
    using ..ESOM
    using NCDatasets
    using .RecordTool

    name = "ESOM"
    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct ESOM_DATA

        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        occ         :: ESOM.OceanColumnCollection

        x2o         :: Dict
        o2x         :: Dict

        configs       :: Dict

        output_vars :: Dict
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

        init_file = (haskey(configs, "init_file")) ? configs["init_file"] : nothing

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
            occ = ESOM.loadSnapshot(init_file; gridinfo_file=configs["domain_file"])
        
        else

            throw(ErrorException("Variable `init_file` is absent in `configs`."))

        end

        wksp = Workspace(occ.Nx, occ.Ny)

        x2o = Dict(
            "SWFLX"  => wksp.swflx,
            "NSWFLX" => wksp.nswflx,
            "TAUX"   => wksp.τx,
            "TAUY"   => wksp.τy,
            "FRWFLX" => wksp.frwflx,
            "IFRAC"  => wksp.ifrac,
        )

        o2x = Dict(
            "SST"      => view(occ.Ts, :, :, 1),
            "QFLX2ATM" => occ.qflx2atm,
        )

        output_vars = Dict(
            #=
            "rain"      => wksp.frwflx,
            "mld"       => occ.h_ML,
            "sst"       => occ.T_ML,
            "qflx2atm"  => occ.qflx2atm,
            "sumflx"    => wksp.sumflx,
            "fric_u"    => wksp.fric_u,
            "ifrac"     => wksp.ifrac,
            =#
        )
       
        
        recorders = Dict()
 
        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]

                recorder = RecordTool.Recorder(
                    Dict(
                        "Nx" => occ.Nx,
                        "Ny" => occ.Ny,
                        "Nz" => 2,
                    ), [
                        ("T",   occ.Ts, ("Nx", "Ny", "Nz")),
                        ("S",   occ.Ss, ("Nx", "Ny", "Nz")),
                        ("M1x", occ.wksp.M1x, ("Nx", "Ny")),
                        ("M1y", occ.wksp.M1y, ("Nx", "Ny")),
                        ("sumflx", wksp.sumflx, ("Nx", "Ny")),
                    ],

                )

                recorders[rec_key] = recorder
            end
        end

        return ESOM_DATA(
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
        MD            :: ESOM_DATA;
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

        wksp.nswflx .*= -1.0
        wksp.swflx  .*= -1.0



        for i=1:MD.occ.Nx, j=1:MD.occ.Ny
            wksp.sumflx[i, j] = wksp.nswflx[i, j] + wksp.swflx[i, j]
            wksp.τx[i, j] *= (1.0 - wksp.ifrac[i, j])
            wksp.τy[i, j] *= (1.0 - wksp.ifrac[i, j])
        end

        ESOM.stepOceanColumnCollection!(
            MD.occ;
            τx     = wksp.τx,
            τy     = wksp.τy,
            swflx  = wksp.swflx,
            nswflx = wksp.nswflx,
            frwflx = wksp.frwflx,
            Δt     = Δt,
        )

        if write_restart
            restart_file = format("restart.ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
            ESOM.takeSnapshot(MD.occ, restart_file)
             
            open(MD.configs["rpointer_file"], "w") do file
                write(file, restart_file)
            end

            println("(Over)write restart pointer file: ", MD.configs["rpointer_file"])
            println("Output restart file: ", restart_file)

        end

    end

    function final(MD::ESOM_DATA)
        
    end

end
