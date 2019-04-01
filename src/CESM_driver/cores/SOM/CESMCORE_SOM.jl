include("../../../SOM/SOM.jl")
module CESMCORE_SOM

    include("../../../share/StatObj.jl")
    include("Workspace_SOM.jl")

    using Formatting
    using ..NetCDFIO
    using ..SOM
    using NCDatasets

    name = "SOM"
    days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    mutable struct SOM_DATA
        map         :: NetCDFIO.MapInfo
        occ         :: SOM.OceanColumnCollection

        x2o         :: Dict
        o2x         :: Dict

        configs       :: Dict

        output_vars :: Dict
        wksp        :: Workspace

        sobj        :: StatObj
        sobj_dict   :: Dict
    end


    function init(;
        map       :: NetCDFIO.MapInfo,
        init_file :: Union{Nothing, AbstractString},
        t         :: AbstractArray{Integer},
        configs   :: Dict,
    )

        if init_file == nothing
            println("No initial ocean profile. Using the naive one.")
            occ = SOM.makeBlankOceanColumnCollection(map.nx, map.ny, map.mask)
            occ.Kh_T = 0.0
            occ.h_ML = 30.0
            occ.T_ML .= 288.0

            snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
            println("Output snapshot: ", snapshot_file)
            SOM.takeSnapshot(occ, joinpath(configs["short_term_archive_dir"], snapshot_file))
            appendLine(configs["short_term_archive_list"], snapshot_file)
        else
            println("Initial ocean with profile: ", init_file)
            occ = SOM.loadSnapshot(init_file)
        end

        wksp = Workspace(occ.Nx, occ.Ny)

        x2o = Dict(
            "SWFLX"  => wksp.swflx,
            "NSWFLX" => wksp.nswflx,
        )

        o2x = Dict(
            "SST"      => occ.T_ML,
            "QFLX2ATM" => occ.qflx2atm,
        )

        output_vars = Dict(
            "sst"       => occ.T_ML,
            "eflx"    => wksp.eflx,
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
            "T_ML"   => occ.T_ML,
            "eflx"   => wksp.eflx,
        )

        return SOM_DATA(
            map,
            occ,
            x2o,
            o2x,
            configs,
            output_vars,
            wksp,
            StatObj(sobj_dict),
            sobj_dict,
        )

    end

    function run(
        MD    :: SOM_DATA;
        t     :: AbstractArray{Integer},
        t_cnt :: Integer,
        Δt    :: Float64,
    )

        if MD.configs["enable_short_term_archive"]
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

                    SOM._createNCFile(MD.occ, joinpath(MD.configs["short_term_archive_dir"], avg_file), MD.map.missing_value)
                    Dataset(avg_file, "a") do ds

                        for v in ["eflx", "T_ML"]
                            SOM._write2NCFile(ds, v, ("Nx", "Ny",), MD.sobj.vars[v], MD.map.missing_value)
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

        SOM.stepOceanColumnCollection!(
            MD.occ;
            eflx   = wksp.eflx,
            Δt     = Δt,
        )

    end

    function final(MD::SOM_DATA)
        
    end

    function appendLine(filename, content)
        open(filename, "a") do io
            write(io, content)
            write(io, "\n")
        end
    end

end
