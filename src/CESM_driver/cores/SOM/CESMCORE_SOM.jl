include("../../../SOM/SOM.jl")
module CESMCORE_SOM

    include("../../../share/StatObj.jl")
    include("../../../share/AppendLine.jl")
    include("Workspace_SOM.jl")

    using Formatting
    using ..NetCDFIO
    using ..SOM
    using NCDatasets

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

        sobjs       :: Dict
        sobj_dict   :: Dict
        rec_cnts    :: Dict
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
            "SST"      => occ.T_ML,
            "QFLX2ATM" => occ.qflx2atm,
        )

        output_vars = Dict(#=
            "T_ML"      => occ.T_ML,
            "h_ML"      => occ.h_ML,
            "eflx"      => wksp.eflx,
            "tfdiv"     => wksp.tfdiv,
            =#
        )
        
        sobj_dict = Dict(
            "T_ML"     => occ.T_ML,
            "swflx"    => wksp.swflx,
            "nswflx"   => wksp.nswflx,
           # "tfdiv"    => wksp.tfdiv,
            "eflx"     => wksp.eflx,
        )


        sobjs    = Dict()
        rec_cnts = Dict()

        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]
                sobjs[rec_key] = StatObj(sobj_dict)
                rec_cnts[rec_key] = 0
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
            sobjs,
            sobj_dict,
            rec_cnts,
        )

    end

    function run(
        MD            :: SOM_DATA;
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

                daily_file = format("{}.xttocn_SOM.h.{:04d}.nc", MD.casename, t[1])
                addStatObj!(MD.sobjs["daily_record"], MD.sobj_dict)

                if t_flags["new_year"]
                    MD.rec_cnts["daily_record"] = 0
                    SOM._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], daily_file), MD.map.missing_value)
                    appendLine(MD.configs["short_term_archive_list"], daily_file)
                end

                if t_flags["new_day"]
                    normStatObj!(MD.sobjs["daily_record"])
                    Dataset(daily_file, "a") do ds
                        for v in keys(MD.sobj_dict)
                            SOM._write2NCFile_time(ds, v, ("Nx", "Ny",), MD.rec_cnts["daily_record"] + 1, MD.sobjs["daily_record"].vars[v]; missing_value = MD.map.missing_value)
                        end
                    end

                    MD.rec_cnts["daily_record"] += 1

                end
 
            end


#=
            if MD.configs["monthly_record"]
                # ===== monthly statistics begin =====

                addStatObj!(MD.sobj, MD.sobj_dict)
                
                # Do monthly average and output it by the end of month
                if days_of_mon[t[2]] == t[3] && t[4] == 0

                    avg_file = format("avg_{:04d}{:02d}.nc", t[1], t[2])
                    
                    normStatObj!(MD.sobj)

                    SOM._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], avg_file), MD.map.missing_value)
                    Dataset(avg_file, "a") do ds

                        for v in keys(MD.sobj_dict)
                            SOM._write2NCFile(ds, v, ("Nx", "Ny",), MD.sobj.vars[v], MD.map.missing_value)
                        end

                    end
                    println("Output monthly average: ", avg_file)
                    appendLine(MD.configs["short_term_archive_list"], avg_file)
                    
                    zeroStatObj!(MD.sobj)
                end
            end
=#
            if MD.configs["yearly_snapshot"]
                # Take snapshot every first day of the year.
                if t[2] == 1 && t[3] == 1 && t[4] == 0
                    snapshot_file = format("Snapshot_{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])

                    SOM.takeSnapshot(MD.occ, joinpath(MD.configs["short_term_archive_dir"], snapshot_file))
                    println("Output snapshot: ", snapshot_file)
                    appendLine(MD.configs["short_term_archive_list"], snapshot_file)
                end
            end
        end

        wksp = MD.wksp
        wksp.nswflx .*= -1.0
        wksp.swflx .*= -1.0

        wksp.eflx[:, :] = wksp.nswflx[:, :]
        wksp.eflx .+= wksp.swflx

        if MD.configs["Qflux"]
            wksp.eflx .+= wksp.tfdiv
        end
        
        SOM.stepOceanColumnCollection!(
            MD.occ;
            eflx   = wksp.eflx,
            Δt     = Δt,
        )
        
        if write_restart
            restart_file = format("restart.xtt_ocn.{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])
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
