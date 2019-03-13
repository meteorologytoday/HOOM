include("default_config.jl")

include("../MLMML/SSM.jl")

using ArgParse
using Printf

#using .NetCDFIO
#using .SSM

module CESM_CORE_MLMML

using Formatting
using ..NetCDFIO
using ..SSM
using ..MLMML

name = "MLMML"

mutable struct MLMML_DATA
    map         :: NetCDFIO.MapInfo
    occ         :: SSM.OceanColumnCollection

    CESM_containers :: Dict
    output_vars     :: Dict

    sst :: AbstractArray{Float64}
    mld :: AbstractArray{Float64}
    qflx2atm :: AbstractArray{Float64} 
    sumflx :: AbstractArray{Float64} 
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

            tmp_oc = MLMML.makeSimpleOceanColumn(;
                zs       = zs,
                b_slope  = init_b_slope,
                b_ML     = init_b_ML,
                h_ML     = MLMML.h_ML_min,
                Δb       = init_Δb,
            )



            SSM.OceanColumnCollection(
                Nx    = map.nx,
                Ny    = map.ny,
                N     = length(zs)-1,
                zs    = zs,
                bs    = tmp_oc.bs,
                K     = K,
                b_ML  = tmp_oc.b_ML,
                h_ML  = tmp_oc.h_ML,
                FLDO  = tmp_oc.FLDO,
                mask  = map.mask
            )
        end
        snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
        println("Output snapshot: ", snapshot_file)
        SSM.takeSnapshot(occ, snapshot_file)

    else
        println("Initial ocean with profile: ", init_file)
        occ = SSM.loadSnapshot(init_file)
    end

    containers = Dict(
        "SWFLX" => occ.wksp.swflx,
        "HFLX"  => occ.wksp.hflx,
        "TAUX"  => occ.wksp.taux,
        "TAUY"  => occ.wksp.tauy,
        "IFRAC"  => occ.wksp.ifrac,
    )

    sst       = zeros(Float64, map.nx, map.ny)
    mld       = copy(sst)
    qflx2atm  = copy(sst)
    sumflx    = copy(sst)

    # Mask data
    SSM.maskData!(occ, sst)
    SSM.maskData!(occ, mld)
    SSM.maskData!(occ, qflx2atm)
    SSM.maskData!(occ, sumflx)

    output_vars = Dict(
        "mld"       => mld,
        "sst"       => sst,
        "sumflx"    => sumflx,
        "qflx2atm"  => qflx2atm,
        "fric_u"    => occ.wksp.fric_u,
        "ifrac"     => occ.wksp.ifrac,
    )
    
    SSM.getInfo!(occ=occ, sst=sst, mld=mld, qflx2atm=qflx2atm)

    return MLMML_DATA(
        map,
        occ,
        containers,
        output_vars,
        sst,
        mld,
        qflx2atm,
        sumflx,
    )

end

function run(
    MD    :: MLMML_DATA;
    t     :: AbstractArray{Integer},
    t_cnt :: Integer,
    Δt    :: Float64,
)

    # Take snapshot every first day of the year.
    if t[2] == 1 && t[3] == 1 && t[4] == 0
        snapshot_file = format("Snapshot_{:04d}0101_00000.nc", t[1])
        println("Output snapshot: ", snapshot_file)
        SSM.takeSnapshot(MD.occ, snapshot_file)
    end
    
    wksp = MD.occ.wksp
    wksp.hflx   .*= -1
    wksp.swflx  .*= -1

    MD.sumflx .= wksp.hflx + wksp.swflx

    SSM.stepOceanColumnCollection!(;
        occ = MD.occ,
        Δt  = Δt,
    )

    SSM.getInfo!(occ=MD.occ; sst=MD.sst, mld=MD.mld, qflx2atm=MD.qflx2atm)

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
