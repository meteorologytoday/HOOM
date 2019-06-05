include(joinpath(@__DIR__, "..", "..", "..", "ESOM", "ESOM.jl"))
module CESMCORE_ESOM

    using Formatting
    using ..NetCDFIO
    using ..ESOM
    using NCDatasets

    name = "ESOM"

    include(joinpath(@__DIR__, "Workspace.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "StatObj.jl"))
    include(joinpath(@__DIR__, "..", "..", "..", "share", "AppendLine.jl"))

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

        sobjs       :: Dict
        sobj_dict   :: Dict
        rec_cnts    :: Dict

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
        
        sobj_dict = Dict(
            "M1x_T1"  => occ.wksp.M1x_T1,
            "M1y_T1"  => occ.wksp.M1y_T1,
            "M1x_T2"  => occ.wksp.M1x_T2,
            "M1y_T2"  => occ.wksp.M1y_T2,
            "T"      => occ.Ts,
            "S"      => occ.Ss,
        )

        sobjs    = Dict()
        rec_cnts = Dict()

        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]
                sobjs[rec_key] = StatObj(sobj_dict)
                rec_cnts[rec_key] = 0
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
            sobjs,
            sobj_dict,
            rec_cnts,
        )

    end

    function run(
        MD            :: ESOM_DATA;
        t             :: AbstractArray{Integer},
        t_cnt         :: Integer,
        t_flags       :: Dict,
        Δt            :: Float64,
        write_restart :: Bool,
    )

        if MD.configs["enable_short_term_archive"]

            if t_cnt == 1 
                for (k, sobj) in MD.sobjs
                    zeroStatObj!(sobj)
                end
            end

            if MD.configs["daily_record"]

                daily_file = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
                addStatObj!(MD.sobjs["daily_record"], MD.sobj_dict)

                # Create new file every month
                if t_flags["new_month"]
                    MD.rec_cnts["daily_record"] = 0
                    ESOM._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], daily_file), MD.map.missing_value)
                    appendLine(MD.configs["short_term_archive_list"], daily_file)
                end

                if t_flags["new_day"]
                    normStatObj!(MD.sobjs["daily_record"])
                    Dataset(daily_file, "a") do ds
                        for (k, v) in MD.sobj_dict

                            if length(size(v)) == 3
                                dim = ("Nx", "Ny", "Nz")
                            elseif length(size(v)) == 2
                                dim = ("Nx", "Ny")
                            end
                            
                            ESOM._write2NCFile_time(ds, k, dim, MD.rec_cnts["daily_record"] + 1, MD.sobjs["daily_record"].vars[k]; missing_value = MD.map.missing_value)
                        end
                    end
                    zeroStatObj!(MD.sobjs["daily_record"])

                    MD.rec_cnts["daily_record"] += 1

                end
 
            end

            if MD.configs["monthly_record"]

                monthly_file = format("{}.ocn.h.monthly.{:04d}.nc", MD.casename, t[1])
                addStatObj!(MD.sobjs["monthly_record"], MD.sobj_dict)

                # Create new file every year 
                if t_flags["new_year"]
                    MD.rec_cnts["monthly_record"] = 0
                    ESOM._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], monthly_file), MD.map.missing_value)
                    appendLine(MD.configs["short_term_archive_list"], monthly_file)
                end

                if t_flags["new_month"] 
                    normStatObj!(MD.sobjs["monthly_record"])
                    Dataset(monthly_file, "a") do ds
                        for (k, v) in MD.sobj_dict

                            if length(size(v)) == 3
                                dim = ("Nx", "Ny", "Nz")
                            elseif length(size(v)) == 2
                                dim = ("Nx", "Ny")
                            end
                            
                            ESOM._write2NCFile_time(ds, k, dim, MD.rec_cnts["monthly_record"] + 1, MD.sobjs["monthly_record"].vars[k]; missing_value = MD.map.missing_value)
                        end
                    end
                    zeroStatObj!(MD.sobjs["monthly_record"])

                    MD.rec_cnts["monthly_record"] += 1

                end
 
            end

        end
 
        wksp = MD.wksp

        wksp.nswflx .*= -1.0
        wksp.swflx  .*= -1.0

        for i=1:MD.occ.Nx, j=1:MD.occ.Ny
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