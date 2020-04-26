
include("../../share/PolelikeCoordinate.jl")
include("../../share/constants.jl")
include("../../share/MapInfo.jl")
include("../../share/RecordTool.jl")


include("Tmd.jl")

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

    run_days = 1
    output_file = "output.nc"

end



z_bnd = collect(Float64, range(0, -100, length=11))

println("Create Gridinfo");
gi = PolelikeCoordinate.RegularCylindricalGridInfo(;
    R = 5000e3,
    Ω = Ωe,
    Nx = 360,
    Ny = 180,
    Ly = 100e3 * 100,
    lat0 = 0.0 |> deg2rad,
    β    = Ωe / Re,
);

Δt = 3600.0

model = Tmd.TmdModel(
    gi = gi,
    Δt = Δt,
    z_bnd = z_bnd,
)



# Setting up recorder
complete_varlist = ShallowWater.getCompleteVariableList(model)
varlist = []
for varname in keys(complete_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(varlist, ( varname, complete_varlist[varname]... ) )
end

coord_varlist = ShallowWater.getCoordinateVariable(model)
cvarlist = []
for varname in keys(coord_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(cvarlist, ( varname, coord_varlist[varname]... ) )
end



#model.state.T[1, :, :] .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat)
#model.state.T[2, :, :] .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat * 2)

#model.state.Φ .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat * 2)
σ = 100e3 * 10.0
#model.state.u_c[1, :, :] .= 1.0

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
    ), varlist, Dict(),
    other_vars = cvarlist
)
 
RecordTool.setNewNCFile!(recorder, output_file)

sum_dσ=sum(gi.dσ)
function integrateT()
    s = 0.0
    for k=1:model.env.Nz_f
        s += sum(gi.dσ .* view(model.state.T, :, :, k))
    end
    return s / sum_dσ / model.env.Nz_f
end

#a = 0.1 * exp.(- (gi.c_y.^2 + gi.c_x.^2) / (σ^2.0) / 2) .* sin.(gi.c_lon*3)
#b = 0.1 * exp.(- (gi.c_y.^2 + gi.c_x.^2) / (σ^2.0) / 2) .* cos.(gi.c_lon*3)
a=b=0
run_days=10

#model.state.v_c[:, 2:end, 1] .= 1.0 * exp.(- (gi.c_y.^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2) .* cos.(gi.c_lon*3)
#model.state.v_c[:, 2:end, 1] .= 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.u_c[:, :, 1] .= 1.0 * exp.(- (gi.c_y.^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.u_c[:, :, 1] .= 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
model.state.Φ[:, :] .= 0.1 * exp.(- (gi.c_y.^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.Φ[:, :] .= g * 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

# output initial state
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)
#println(integrateT())


# 0.5 * sin.((model.env.gi.c_lon .- deg2rad(45))) .* cos.(model.env.gi.c_lat)

@time for step=1:run_days
    println("Run day ", step)

    if step <= 10
        model.state.τx[:, :] .= a * exp(- (step-.5)/ 10.0)
        model.state.τy[:, :] .= b * exp(- (step-.5)/ 10.0)
    else
        model.state.τx[:, :] .= 0
        model.state.τy[:, :] .= 0
    end

    for substep = 1:Int64(86400/Δt)
        ShallowWater.stepModel!(model)
        RecordTool.record!(recorder)
        RecordTool.avgAndOutput!(recorder)
    end

end

