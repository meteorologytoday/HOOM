
include("../../share/GridFiles.jl")
include("../../share/PolelikeCoordinate.jl")
include("../../share/constants.jl")

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

    run_days = 10
    output_file = "output_tmd.nc"

end



z_bnd = collect(Float64, range(0, -100, length=11))

println("Create Gridinfo");

#=
gi = PolelikeCoordinate.RegularCylindricalGridInfo(;
    R = 5000e3,
    Ω = Ωe,
    Nx = 60,
    Ny = 30,
    Ly = 100e3 * 100,
    lat0 = 0.0 |> deg2rad,
    β    = Ωe / Re,
);
=#

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx1v6.nc"


hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv0.9x1.25.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv4x5_gx3v7.091218.nc"


gf = GridFiles.CurvilinearSphericalGridFile(
        hrgrid_file;
        R   = Re,
        Ω   = Ωe,
)


xcutoff = 1
ycutoff = 5

gf.mask                       .= 1
gf.mask[:, 1:ycutoff]         .= 0 
gf.mask[:, end-ycutoff+1:end] .= 0 

gf.mask[1:xcutoff, :]         .= 0 
gf.mask[end-xcutoff+1:end, :] .= 0 



Dataset(topo_file, "r") do ds
    global mask_idx = (ds["depth"][:] |> nomissing) .< 1000.0
end

gi = PolelikeCoordinate.genGridInfo(gf);

Δt = 3600.0 * 24 

model = Tmd.TmdModel(
    gi                    = gi,
    Δt                    = Δt,
    substeps              = 8,
    z_bnd                 = z_bnd,
    mask2                 = gf.mask,
    Kh_X                  = [0.0, 0.0],
    Kv_X                  = [0.0, 0.0],
    MLT_rng               = [3, 1e5],
    radiation_scheme      = :exponential_decay,
    MLT_scheme            = :prescribe_MLT,
    convective_adjustment = true,
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


model.state.h_ML .= 7.0

model.state.S_ML .= 35.0
model.state.S .= 35.0

model.state.T .= 10.0
model.state.T_ML .= 10.0
model.state.T_ML .= 20 .+ 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat)
#model.state.T[1, :, :] .= model.state.T_ML


#model.forcing.swflx .= -1000.0
model.forcing.h_ML  .= model.state.h_ML
#model.state.Φ .= 1.0 * 10.0 * sin.((model.env.gi.c_lon )) .* sin.(model.env.gi.c_lat * 2)
σ = 100e3 * 10.0
#model.state.u_c[1, :, :] .= 1.0

model.forcing.u_U[1, :, :] .= 1 .* cos.(2*model.env.gi.c_lat)#* sin.((model.env.gi.c_lon )) .* cos.(model.env.gi.c_lat*2) #.* exp.( - (model.env.gi.c_y / σ).^2.0 / 2.0 )

Tmd.initialization!(model)

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
    for k=1:model.env.Nz
        s += sum(gi.dσ .* view(model.state.T, k, :, :))
    end
    return s / sum_dσ / model.env.Nz
end

# output initial state
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)

println(integrateT())

@time for step=1:run_days
    
    println("Run day ", step)

    for substep = 1:model.env.substeps
        Tmd.stepModel!(model)
    end
    RecordTool.record!(recorder)
    RecordTool.avgAndOutput!(recorder)

end

