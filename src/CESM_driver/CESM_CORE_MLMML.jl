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

    sst :: AbstractArray{Float64}
    mld :: AbstractArray{Float64}
    qflux2atm :: AbstractArray{Float64} 
    sumflx :: AbstractArray{Float64} 
end


function init(map::NetCDFIO.MapInfo)

    init_b_ML     = 280.0 * MLMML.g * MLMML.α
    init_h_ML     = MLMML.h_ML_min
    init_b_slope  = 30.0 / 5000.0 * MLMML.g * MLMML.α
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
        "SWFLX" => zeros(Float64, map.lsize),
        "HFLX"  => zeros(Float64, map.lsize),
        "TAUX"  => zeros(Float64, map.lsize),
        "TAUY"  => zeros(Float64, map.lsize),
    )


    sst       = zeros(Float64, map.lsize)
    mld       = copy(sst)
    qflux2atm = copy(sst)
    sumflx      = copy(sst)

    # Mask data
    SSM.maskData!(occ, sst)
    SSM.maskData!(occ, mld)
    SSM.maskData!(occ, qflux2atm)
    SSM.maskData!(occ, sumflx)

    output_vars = Dict(
        "mld"       => reshape(mld,       map.nx, map.ny),
        "sst"       => reshape(sst,       map.nx, map.ny),
        "sumflx"    => reshape(sumflx,    map.nx, map.ny),
        "qflux2atm" => reshape(qflux2atm, map.nx, map.ny),
        "fric_u"    => reshape(occ.wksp.fric_u, map.nx, map.ny),
    )

       
    return MLMML_DATA(
        map,
        occ,
        containers,
        sst,
        mld,
        qflux2atm,
        sumflx
    )

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
