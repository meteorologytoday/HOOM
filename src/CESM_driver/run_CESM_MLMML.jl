include("default_config.jl")

include("../MLMML/MLMML.jl")

using ArgParse
using Printf

module CESM_CORE_MLMML

using Formatting
using ..NetCDFIO
using ..MLMML
using NCDatasets

name = "MLMML"

include("Workspace_MLMML.jl")
include("../share/StatObj.jl")

days_of_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

mutable struct MLMML_DATA
    map         :: NetCDFIO.MapInfo
    occ         :: MLMML.OceanColumnCollection

    x2o         :: Dict
    o2x         :: Dict

    output_vars :: Dict
    wksp        :: Workspace

    sobj        :: StatObj
end


function init(;
    map       :: NetCDFIO.MapInfo,
    init_file :: Union{Nothing, AbstractString},
    t         :: AbstractArray{Integer},
)

    if init_file == nothing
        println("No initial ocean profile. Using the naive one.")
        occ = let

            zs = collect(Float64, range(0, -500, step=-5))
            K = 1e-5

            init_b_ML     = (280.0 - MLMML.T_ref) * MLMML.g * MLMML.α
            init_h_ML     = MLMML.h_ML_min
            init_b_slope  = 30.0 / 5000.0 * MLMML.g * MLMML.α
            init_Δb       = 1.0 * MLMML.g * MLMML.α

            MLMML.makeBasicOceanColumnCollection(
                map.nx, map.ny, zs;
                b_slope  = init_b_slope,
                b_ML     = init_b_ML,
                h_ML     = MLMML.h_ML_min,
                Δb       = init_Δb,
                K        = K,
                mask     = map.mask,
            )
        end

        snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
        println("Output snapshot: ", snapshot_file)
        MLMML.takeSnapshot(occ, snapshot_file)

    else
        println("Initial ocean with profile: ", init_file)
        occ = MLMML.loadSnapshot(init_file)
    end

    wksp = Workspace(occ.Nx, occ.Ny, occ.Nz)

    x2o = Dict(
        "SWFLX" => wksp.swflx,
        "HFLX"  => wksp.hflx,
        "TAUX"  => wksp.taux,
        "TAUY"  => wksp.tauy,
        "IFRAC" => wksp.ifrac,
    )

    o2x = Dict(
        "SST"      => occ.sst,
        "QFLX2ATM" => occ.qflx2atm,
    )

    output_vars = Dict(
        "mld"       => occ.h_ML,
        "sst"       => occ.sst,
        "qflx2atm"  => occ.qflx2atm,
        "sumflx"    => wksp.sumflx,
        "fric_u"    => wksp.fric_u,
        "ifrac"     => wksp.ifrac,
    )
    
    return MLMML_DATA(
        map,
        occ,
        x2o,
        o2x,
        output_vars,
        wksp,
        StatObj(Dict(
            "mld" => occ.h_ML,
            "sst" => occ.sst,
            "sumflx" => wksp.sumflx,
            "fric_u" => wksp.fric_u,
        )),
    )

end

function run(
    MD    :: MLMML_DATA;
    t     :: AbstractArray{Integer},
    t_cnt :: Integer,
    Δt    :: Float64,
)


    if t_cnt == 1 
        zeroStatObj!(MD.sobj)
    end

    addStatObj!(MD.sobj, Dict(
        "mld" => MD.occ.h_ML,
        "sst" => MD.occ.sst,
        "sumflx" => MD.wksp.sumflx,
        "fric_u" => MD.wksp.fric_u,
    ))
    
    # Do monthly average and output it by the end of month
    if days_of_mon[t[2]] == t[3] && t[4] == 0
        avg_file = format("avg_{:04d}{:02d}.nc", t[1], t[2])
        
        normStatObj!(MD.sobj)

        MLMML._createNCFile(MD.occ, avg_file, MD.map.missing_value)

        Dataset(avg_file, "a") do ds
            for v in ("mld", "sst", "sumflx", "fric_u")
                MLMML._write2NCFile(ds, v,    ("Nx", "Ny",), MD.sobj.vars[v], MD.map.missing_value)
            end
        end
        println("Output monthly average: ", avg_file)
        
       zeroStatObj!(MD.sobj)
    end


    # Take snapshot every first day of the year.
    if t[2] == 1 && t[3] == 1 && t[4] == 0
        snapshot_file = format("Snapshot_{:04d}{:02d}{:02d}_00000.nc", t[1], t[2], t[3])
        MLMML.takeSnapshot(MD.occ, snapshot_file)
        println("Output snapshot: ", snapshot_file)
    end
    
    wksp = MD.wksp

    wksp.hflx   .*= -1
    wksp.swflx  .*= -1

    wksp.sumflx .= wksp.hflx + wksp.swflx
    
    wksp.fric_u .= sqrt.(sqrt.((wksp.taux).^2.0 .+ (wksp.tauy).^2.0) / MLMML.ρ)
    wksp.weighted_fric_u .*= (1.0 .- wksp.ifrac)

    wksp.hflx   .*= MLMML.αgρc
    wksp.swflx  .*= MLMML.αgρc
    
    MLMML.stepOceanColumnCollection!(
        MD.occ;
        fric_u = wksp.weighted_fric_u,
        B0     = wksp.hflx,
        J0     = wksp.swflx,
        Δt     = Δt,
    )



end

function final(MD::MLMML_DATA)
    
end

end

# ===== READ ALL CONFIG BEGIN =====
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "Configuration file."
            arg_type = String

        "--init-file"
            help = "Initial profile of ocean. If not provided then a default profile will be used."
            arg_type = String

    end

    return parse_args(s)
end

parsed_args = parse_commandline()
#for (arg,val) in parsed_args
#    println("$arg  =>  ", repr(val))
#end

config_file = parsed_args["config"]
if config_file != nothing
    config_file = normpath( (isabspath(config_file)) ? config_file : joinpath(pwd(), config_file) )
    println("Load config file: ", config_file)
    include(config_file)
end

init_file = parsed_args["init-file"]
if init_file != nothing
    init_file = normpath( (isabspath(init_file)) ? init_file : joinpath(pwd(), init_file) )
    println("Ocean init file : ", init_file)

    if !isfile(init_file)
        throw(ErrorException("File missing: ", init_file))
    end
end

# ===== READ ALL CONFIG END =====

OMMODULE = Main.CESM_CORE_MLMML



include("driver.jl")
