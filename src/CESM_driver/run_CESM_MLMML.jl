include("default_config.jl")

include("../MLMML/SSM.jl")

using ArgParse
using Printf

#using .NetCDFIO
#using .SSM

module CESM_CORE_MLMML

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
    map :: NetCDFIO.MapInfo,
    t   :: AbstractArray{Integer},
)

    init_b_ML     = (280.0 - MLMML.T_ref) * MLMML.g * MLMML.α
    init_h_ML     = MLMML.h_ML_min
    init_b_slope  = 10.0 / 5000.0 * MLMML.g * MLMML.α
    init_Δb       = 1.0 * MLMML.g * MLMML.α

    zs = collect(Float64, range(0, -500, step=-5))
    K = 1e-5

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
 
    containers = Dict(
        "SWFLX" => occ.wksp.swflx,
        "HFLX"  => occ.wksp.hflx,
        "TAUX"  => occ.wksp.taux,
        "TAUY"  => occ.wksp.tauy,
        "IFRAC"  => occ.wksp.ifrac,
    )

    sst       = zeros(Float64, map.lsize)
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

OMMODULE = Main.CESM_CORE_MLMML



include("driver.jl")
