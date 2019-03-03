include("../MLMML/SSM.jl")
include("julia_lib/MailboxPipe.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")

using ArgParse
using Printf
using Formatting
using JSON

using .BinaryIO
using .MailboxPipe
using .NetCDFIO
using .SSM

module CESM_CORE_MLMML

mutable struct MLMML_DATA
    map         :: NetCDFIO.MapInfo
    occ         :: SSM.OceanColumnCollection
end


function init(map::NetCDFIO.MapInfo)

    init_b_ML     = 280.0 * MLMML.g * MLMML.α
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


    occ = SSM.OceanColumnCollection(
        N_ocs = map.lsize,
        N     = length(zs)-1,
        zs    = zs,
        bs    = tmp_oc.bs,
        K     = K,
        b_ML  = tmp_oc.b_ML,
        h_ML  = tmp_oc.h_ML,
        FLDO  = tmp_oc.FLDO,
        mask  = map.mask
    )
    
    return MLMML_DATA(map, occ)

end

function run(
    MD    :: MLMML_DATA;
    t     :: AbstractArray{Float64},
    t_cnt :: Integer,
    Δt    :: Float64,
)

    wksp = MD.occ.wksp
    wksp.fric_u .= sqrt.(sqrt.((wksp.taux).^2.0 + (wksp.tauy).^2.0) / MLMML.ρ)
    wksp.hflx   *= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)
    wksp.swflx  *= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)

    SSM.stepOceanColumnCollection!(;
        occ = MD.occ,
        Δt  = Δt,
    )
end

function final(MD::MLMML_DATA)
    
end

end

include("default_config.jl")

# ===== READ ALL CONFIG BEGIN =====
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "Configuration file."
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
# ===== READ ALL CONFIG END =====


include("driver.jl")
