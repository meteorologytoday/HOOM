
include("../../../share/constants.jl")
include("../../../share/MapInfo.jl")
include("../../../share/RecordTool.jl")

include("ShallowWater.jl")

using .ShallowWater

using NCDatasets
using Formatting
using ArgParse
using JSON



if !isdefined(Main, :REPL)

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

else

    run_days = 5
    output_file = "output.nc"

end



z_bnd_f = collect(Float64, range(0, -200, length=11))
height_level_counts = [2, 8]
Dh = [1e3, 1e3]
Dv = [1e-3, 1e-3] * 0.0


mi = ModelMap.MapInfo{Float64}("/seley/tienyiah/CESM_domains/domain.lnd.fv4x5_gx3v7.091218.nc")
cutoff = 5
gi = PolelikeCoordinate.CurvilinearSphericalGridInfo(;
    R=Re,
    Ω=Ωe,
    Nx=mi.nx,
    Ny=mi.ny,
    c_lon=mi.xc,
    c_lat=mi.yc,
    vs_lon=mi.xv,
    vs_lat=mi.yv,
    area=mi.area,
    angle_unit=:deg,
    sub_yrng = cutoff+1:mi.ny-cutoff,
)

println(rad2deg.(gi.c_lat[1, :]))

model = ShallowWater.DynModel(
    gi = gi,
    Δt = 86400.0,
    z_bnd_f = z_bnd_f,
    height_level_counts = height_level_counts,
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
 
RecordTool.setNewNCFile!(recorder, output_file)

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
for step=1:run_days
    println("Run day ", step)

    for substep = 1:1
        ShallowWater.stepModel!(model)
        RecordTool.record!(recorder)
    end

    RecordTool.avgAndOutput!(recorder)
    println(integrateT())
end

