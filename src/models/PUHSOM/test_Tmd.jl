
include("../../share/PolelikeCoordinate.jl")
include("../../share/constants.jl")
include("../../share/MapInfo.jl")
include("../../share/RecordTool.jl")

include("Tmd/Tmd.jl")

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

    run_days = 2
    output_file = "output_tmd.nc"

end



z_bnd = collect(Float64, range(0, -100, length=11))

println("Create Gridinfo");
gi = PolelikeCoordinate.RegularCylindricalGridInfo(;
    R = 5000e3,
    Ω = Ωe,
    Nx = 120,
    Ny = 100,
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
complete_varlist = Tmd.getCompleteVariableList(model)
varlist = []
for varname in keys(complete_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(varlist, ( varname, complete_varlist[varname]... ) )
end


coord_varlist = Tmd.getCoordinateVariable(model)
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
        "Nz_p1"  => model.env.Nz+1,
        "Nx"      => model.env.Nx,
        "Ny"      => model.env.Ny,
        "Nz"    => model.env.Nz,
        "z_bnd" => length(model.env.z_bnd),
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

# output initial state
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)


@time for step=1:run_days
    
    println("Run day ", step)

    for substep = 1:Int64(86400/Δt)
        Tmd.stepModel!(model)
        RecordTool.record!(recorder)
        RecordTool.avgAndOutput!(recorder)
    end

end

