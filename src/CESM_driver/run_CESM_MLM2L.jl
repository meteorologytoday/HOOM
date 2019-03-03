include("default_config.jl")

include("../MLM2L/MLM2L.jl")

using ArgParse
using Printf

module CESM_DRIVER_MLM2L

using ..NetCDFIO
using ..MLM2L

name = "MLM2L"

mutable struct MLM2L_DATA
    map         :: NetCDFIO.MapInfo
    occ         :: MLM2L.OceanColumnCollection
    
    CESM_containers :: Dict
    output_vars     :: Dict

    sst :: AbstractArray{Float64}
    mld :: AbstractArray{Float64}
    qflux2atm :: AbstractArray{Float64} 
    sumflx :: AbstractArray{Float64}

end


function init(
    map::NetCDFIO.MapInfo,
    t     :: AbstractArray{Integer},
)

    daycnt = date2daycnt(t[2], t[3])

    init_b_ML     = 280.0 * MLM2L.g * MLM2L.α
    init_h_ML     = MLM2L.h_ML_min

    b_DO          = (273.15 + 2.0) * MLM2L.g * MLM2L.α

    occ = MLM2L.makeBlankOceanColumnCollection(
        N_ocs    = map.lsize,
        period_n = 365,
        period   = 86400.0 * 365,
        t        = collect(Float64, 0:364) * 86400.0,
        mask  = map.mask,
    )
 
    sst       = zeros(Float64, map.lsize)
    mld       = copy(sst)
    qflux2atm = copy(sst)
    sumflx    = copy(sst)

    CESM_containers = Dict(
        "SWFLX" => copy(sst),
        "HFLX"  => copy(sst),
    )

    output_vars = Dict(
        "mld"       => mld,
        "sst"       => sst,
        "sumflx"    => sumflx,
        "qflux2atm" => qflux2atm,
    )
    
    MLM2L.getInfo!(occ=occ, sst=sst, mld=mld, idx=daycnt)

    return MLM2L_DATA(
        map,
        occ,
        CESM_containers,
        output_vars,
        sst,
        mld,
        qflux2atm,
        sumflx,
    )

end

function run(
    MD    :: MLM2L_DATA;
    t     :: AbstractArray{Integer},
    t_cnt :: Integer,
    Δt    :: Float64,
)
    daycnt = date2daycnt(t[2], t[3]),
   
    MD.sumflx .=  MD.CESM_containers["HFLX"]
    MD.sumflx .+= MD.CESM_containers["SWFLX"]

    MLM2L.stepOceanColumnCollection!(;
        occ = MD.occ,
        F   = MD.sumflx,
        idx = daycnt,
        Δt  = Δt,
    )

    SSM.getInfo!(occ=MD.occ; sst=MD.sst, mld=MD.mld, idx=daycnt)

end

function final(MD::MLM2L_DATA)
    
end

d_of_m = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
beg_of_m = zeros(Integer, length(d_of_m))
beg_of_m[1] = 1
for i = 2:length(d_of_m)
    beg_of_m[i] = beg_of_m[i-1] + d_of_m[i-1]
end

function date2daycnt(
    m :: Integer,
    d :: Integer,
)
    if m < 1 || m >length(d_of_m) || d < 1 || d > d_of_m[m]
        throw(ErrorException("[date2daycnt] invalid input."))
    end

    return beg_of_m[m] + d - 1
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

OMMODULE = Main.CESM_DRIVER_MLM2L



include("driver.jl")
