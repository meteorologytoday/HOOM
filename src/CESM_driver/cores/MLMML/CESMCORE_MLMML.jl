include("../../../MLMML/MLMML.jl")
module CESMCORE_MLMML

    using Formatting
    using ..NetCDFIO
    using ..MLMML
    using NCDatasets

    zs = collect(Float64, range(0, -500, step=-5))
    name = "MLMML"

    include("Workspace_MLMML.jl")
    include("../../../share/StatObj.jl")
    include("../../../share/AppendLine.jl")

    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct MLMML_DATA

        casename    :: AbstractString
        map         :: NetCDFIO.MapInfo
        occ         :: MLMML.OceanColumnCollection

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
        init_file    :: Union{Nothing, AbstractString},
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
 
        else

            must_need_keys = [
                "z_coord",
            ]

            optional_keys = [
                "climatology_temperature_file",
                "climatology_temperature_relaxation_time",
                "climatology_salinity_file",
                "climatology_salinity_relaxation_time",
                "topography_file",
            ]


            for key in must_need_keys
                if !haskey(configs, key)
                    throw(ErrorException(format("Variable config must contain key: {}", key)))
                end
            end

            println("Optional keys are: ", optional_keys)

            if haskey(configs, "climatology_temperature_file")
                println("Using climatology temperature file: ", configs["climatology_temperature_file"])
            end

            if haskey(configs, "climatology_salinity_file")
                println("Using climatology salinity file: ", configs["climatology_salinity_file"])
            end

            if haskey(configs, "topography_file")
                println("Using topography file: ", configs["topography_file"])
            end

        end


        if typeof(init_file) <: AbstractString

            println("Initial ocean with profile: ", init_file)
            occ = MLMML.loadSnapshot(init_file)

        else
            println("No initial ocean profile. Using the naive one.")
            occ = let

                zs = collect(Float64, range(0, -500, step=-5))
                K_T = 1e-5
                K_S = 1e-5
                h_ML_min = 10.0
                h_ML_max = 499.0
                we_max   = 1e-2

                init_h_ML     = h_ML_min
                
                init_T_ML     = 288.0
                init_T_slope  = 2.0 / 4000.0
                init_ΔT       = 5.0

                init_S_ML     = MLMML.S_ref 
                init_S_slope  = 0.0
                init_ΔS       = 0.0

                MLMML.makeBasicOceanColumnCollection(
                    map.nx, map.ny, zs;
                    T_ML     = init_T_ML,
                    ΔT       = init_ΔT,
                    T_slope  = init_T_slope,
                    S_ML     = init_S_ML,
                    ΔS       = init_ΔS,
                    S_slope  = init_S_slope,
                    h_ML     = h_ML_min,
                    K_T      = K_T,
                    K_S      = K_S,
                    h_ML_min = h_ML_min,
                    h_ML_max = h_ML_max,
                    we_max   = we_max,
                    mask     = map.mask,
                )
            end

            snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
            println("Output snapshot: ", snapshot_file)
            MLMML.takeSnapshot(occ, joinpath(configs["short_term_archive_dir"], snapshot_file))
            appendLine(configs["short_term_archive_list"], snapshot_file)
        end

        wksp = Workspace(occ.Nx, occ.Ny, occ.Nz)

        x2o = Dict(
            "SWFLX"  => wksp.swflx,
            "NSWFLX" => wksp.nswflx,
            "TAUX"  => wksp.taux,
            "TAUY"  => wksp.tauy,
            "IFRAC" => wksp.ifrac,
            "FRWFLX" => wksp.frwflx,
        )

        o2x = Dict(
            "SST"      => occ.T_ML,
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
            "mld"    => occ.h_ML,
            "T"      => occ.Ts,
            "S"      => occ.Ss,
            "sumflx" => wksp.sumflx,
            "fric_u" => wksp.fric_u,
            "frwflx" => wksp.frwflx,
        )

        sobjs    = Dict()
        rec_cnts = Dict()

        for rec_key in ["daily_record", "monthly_record"]
            if configs[rec_key]
                sobjs[rec_key] = StatObj(sobj_dict)
                rec_cnts[rec_key] = 0
            end
        end

        return MLMML_DATA(
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
        MD            :: MLMML_DATA;
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
                    MLMML._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], daily_file), MD.map.missing_value)
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
                            
                            MLMML._write2NCFile_time(ds, k, dim, MD.rec_cnts["daily_record"] + 1, MD.sobjs["daily_record"].vars[k]; missing_value = MD.map.missing_value)
                        end
                    end

                    MD.rec_cnts["daily_record"] += 1

                end
 
            end

#=
            if MD.configs["monthly_record"]
                # ===== monthly statistics begin =====
                if t_cnt == 1 
                    zeroStatObj!(MD.sobj)
                end

                addStatObj!(MD.sobj, MD.sobj_dict)
                
                # Do monthly average and output it by the end of month
                if days_of_mon[t[2]] == t[3] && t[4] == 0

                    avg_file = format("avg_{:04d}{:02d}.nc", t[1], t[2])
                    
                    normStatObj!(MD.sobj)

                    MLMML._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], avg_file), MD.map.missing_value)
                    Dataset(avg_file, "a") do ds

                        for v in ["mld", "sumflx", "fric_u", "frwflx"]
                            MLMML._write2NCFile(ds, v, ("Nx", "Ny",), MD.sobj.vars[v], MD.map.missing_value)
                        end

                        for v in ["T", "S"]
                            MLMML._write2NCFile(ds, v, ("Nx", "Ny", "Nz"), MD.sobj.vars[v], MD.map.missing_value)
                        end

                    end
                    println("Output monthly average: ", avg_file)
                    appendLine(MD.configs["short_term_archive_list"], avg_file)
                    
                    zeroStatObj!(MD.sobj)
                end
            end

            if MD.configs["yearly_snapshot"]
                # Take snapshot every first day of the year.
                if t[2] == 1 && t[3] == 1 && t[4] == 0
                    snapshot_file = format("Snapshot_{:04d}{:02d}{:02d}_{:05d}.nc", t[1], t[2], t[3], t[4])

                    MLMML.takeSnapshot(MD.occ, joinpath(MD.configs["short_term_archive_dir"], snapshot_file))
                    println("Output snapshot: ", snapshot_file)
                    appendLine(MD.configs["short_term_archive_list"], snapshot_file)
                end
            end
=#
        end
 
        wksp = MD.wksp

        wksp.nswflx .*= -1.0
        wksp.swflx  .*= -1.0

        #wksp.sumflx[:, :]  = wksp.nswflx
        #wksp.sumflx      .+= wksp.swflx
        
        wksp.fric_u .= sqrt.(sqrt.((wksp.taux).^2.0 .+ (wksp.tauy).^2.0) / MLMML.ρ)
        wksp.weighted_fric_u .*= (1.0 .- wksp.ifrac)

        MLMML.stepOceanColumnCollection!(
            MD.occ;
            fric_u = wksp.weighted_fric_u,
            swflx  = wksp.swflx,
            nswflx = wksp.nswflx,
            frwflx = wksp.frwflx,
            Δt     = Δt,
        )



    end

    function final(MD::MLMML_DATA)
        
    end

end
