include("../../../share/constants.jl")
include("../../../share/MapInfo.jl")
include("../../../share/RecordTool.jl")

include("ShallowWater.jl")

using .ShallowWater

using NCDatasets
using Formatting
using ArgParse
using JSON




function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--run-days"
            help = "Simulation days."
            arg_type = Int64
            default = 10
 
        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
    end

    return parse_args(ARGS, s)
end



parsed = parse_commandline()
print(json(parsed, 4))





z_bnd_f = collect(Float64, range(0, -200, length=21))
height_level_counts = [2, 2, 6, 10]
Dh = [1e3, 1e3]
Dv = [1e-3, 1e-3] * 0.0


mi = ModelMap.MapInfo{Float64}("/seley/tienyiah/CESM_domains/domain.lnd.fv4x5_gx3v7.091218.nc")
cutoff = 5
gi = DisplacedPoleCoordinate.GridInfo(
    Re,
    mi.nx,
    mi.ny,
    mi.xc,
    mi.yc,
    mi.xv,
    mi.yv,
    mi.area;
    angle_unit=:deg,
    sub_yrng = cutoff+1:mi.ny-cutoff,
)

println(rad2deg.(gi.c_lat[1, :]))

model = ShallowWater.Model(
    gi = gi,
    z_bnd_f = z_bnd_f,
    height_level_counts = height_level_counts,
    Dh = Dh,
    Dv = Dv,
    f = nothing,
    ϵ = nothing,
)

# Setting up recorder
complete_varlist = ShallowWater.getCompleteVariableList(model)
varlist = []
for varname in keys(complete_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(varlist, ( varname, complete_varlist[varname]... ) )
end


model.state.T[1, :, :] .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat)
model.state.T[2, :, :] .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat * 2)





recorder = RecordTool.Recorder(
    Dict(
        "NX"      => model.env.NX,
        "Nyp1"    => model.env.Ny+1,
        "Nz_fp1"  => model.env.Nz_f+1,
        "Nx"      => model.env.Nx,
        "Ny"      => model.env.Ny,
        "Nz_f"    => model.env.Nz_f,
        "Nz_c"    => model.env.Nz_c,
        "z_bnd_f" => length(model.env.z_bnd_f),
    ), varlist, Dict()
)
 
RecordTool.setNewNCFile!(recorder, parsed["output-file"])

sum_dσ=sum(gi.dσ)
function integrateT()
    s = 0.0
    for k=1:model.env.Nz_f
        s += sum(gi.dσ .* view(model.state.T, k, :, :))
    end
    return s / sum_dσ / model.env.Nz_f
end

# output initial state
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)
println(integrateT())


model.state.u_f[1, :, :] .= 0.5 * sin.((model.env.gi.c_lon .- deg2rad(45))) .* cos.(model.env.gi.c_lat)
for step=1:parsed["run-days"]
    println("Run day ", step)

    for substep = 1:1
        ShallowWater.stepModel!(model, 86400.0/1)
        RecordTool.record!(recorder)
    end

    RecordTool.avgAndOutput!(recorder)
    println(integrateT())
end

