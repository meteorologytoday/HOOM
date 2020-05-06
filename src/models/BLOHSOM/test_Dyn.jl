
include("../../share/constants.jl")
include("../../share/GridFiles.jl")
include("../../share/RecordTool.jl")

include("Dyn/Dyn.jl")

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



z_bnd_f = collect(Float64, range(0, -4000, length=11))
height_level_counts = [1, 3, 6]
height_level_counts = [10]

println("Create Gridinfo");
hrgrid_file = "/seley/tienyiah/CESM_domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc"

gf_ref = GridFiles.CurvilinearSphericalGridFile(
        hrgrid_file;
        R   = Re,
        Ω   = Ωe,
)

Ly = 100e3 * 150.0
gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = gf_ref.Nx,
        Ny   = gf_ref.Ny,
        Ly   = Ly,
        lat0 = 0.0 |> deg2rad,
        β    = Ωe / Re,
)

#xcutoff = 10
ycutoff = 2

gf.mask .= 1.0 .- gf_ref.mask

#gf.mask                       .= 1
gf.mask[:, 1:ycutoff]         .= 0 
gf.mask[:, end-ycutoff+1:end] .= 0 

#gf.mask[1:xcutoff, :]         .= 0 
#gf.mask[end-xcutoff+1:end, :] .= 0 



#=
hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx1v6.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/domain.ocn.gx1v6.090206.nc"
topo_file = "/seley/tienyiah/CESM_domains/ocean_topog_gx1v6.nc"


hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv0.9x1.25.nc"


gf = GridFiles.CurvilinearSphericalGridFile(
        hrgrid_file;
        R   = Re,
        Ω   = Ωe,
)

Dataset(topo_file, "r") do ds
    global mask_idx = (ds["depth"][:] |> nomissing) .< 1000.0
end

gf.mask = 1 .- gf.mask
gf.mask[mask_idx] .= 0.0

=#


gi = PolelikeCoordinate.genGridInfo(gf);


#=
mi = ModelMap.MapInfo{Float64}("/seley/tienyiah/CESM_domains/domain.lnd.fv4x5_gx3v7.091218.nc")
cutoff = 5
gi = PolelikeCoordinate.CurvilinearSphericalGridInfo(;
    R=Re,
    Ω=Ωe,
    Nx=gf.nx,
    Ny=gf.ny,
    c_lon=gf.xc,
    c_lat=gf.yc,
    vs_lon=gf.xv,
    vs_lat=gf.yv,
    area=gf.area,
    angle_unit=:deg,
    sub_yrng = cutoff+1:gf.ny-cutoff,
)
=#

#println(rad2deg.(gi.c_lat[1, :]))


Δt = 3600.0

model = Dyn.DynModel(
    gi                  = gi,
    Δt                  = Δt,
    Dh                  = 10000.0,
    z_bnd_f             = z_bnd_f,
    height_level_counts = height_level_counts,
    mask                = gf.mask,
)

#=
U = similar(model.state.U)
V = similar(model.state.V)
T = similar(model.state.Φ)

U.=1000
V.=1000
T.=1000

U_filtered = copy(U)
V_filtered = copy(V)
T_filtered = copy(T)


Dyn.mul2!(U_filtered, model.core.s_ops.filter_U, U)
Dyn.mul2!(V_filtered, model.core.s_ops.filter_V, V)
Dyn.mul2!(T_filtered, model.core.s_ops.filter_T, T)



Dataset("dissect.nc", "c") do ds

    defDim(ds, "Nx", gf.nx)
    defDim(ds, "Ny", gf.ny)
    defDim(ds, "Nyp1", gf.ny+1)

    for (varname, vardata, vardim, attrib) in [
        ("mask", gf.mask, ("Nx", "Ny"), Dict()),
        ("U_filtered", U_filtered, ("Nx", "Ny"), Dict()),
        ("V_filtered", V_filtered, ("Nx", "Nyp1"), Dict()),
        ("T_filtered", T_filtered, ("Nx", "Ny"), Dict()),
    ]

        println("Doing var: ", varname)

        var = defVar(ds, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end
=# 
#println("Pause")
#readline()
# Setting up recorder
complete_varlist = Dyn.getCompleteVariableList(model)
varlist = []
for varname in keys(complete_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(varlist, ( varname, complete_varlist[varname]... ) )
end

coord_varlist = Dyn.getCoordinateVariable(model)
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
function integrateKE()
    s = 0.0
    KE = (model.state.U.^2 + model.state.V.^2)
    return avg(KE)
end

function avg(x)
    return sum(gi.dσ .* x) / sum_dσ
end


#a = 0.1 * exp.(- (gi.c_y.^2 + gi.c_x.^2) / (σ^2.0) / 2) .* sin.(gi.c_lon*3)
#b = 0.1 * exp.(- (gi.c_y.^2 + gi.c_x.^2) / (σ^2.0) / 2) .* cos.(gi.c_lon*3)
a=b=0
run_days=100

#model.state.v_c[:, 2:end, 1] .= 1.0 * exp.(- (gi.c_y.^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2) .* cos.(gi.c_lon*3)
#model.state.v_c[:, 2:end, 1] .= 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.u_c[:, :, 1] .= 1.0 * exp.(- (gi.c_y.^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.u_c[:, :, 1] .= 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.Φ[:, :] .= 0.01 * exp.(- ((gi.c_lat * gi.R).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
model.state.Φ[:, :] .= 0.01 * exp.(- ((gi.c_y).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#model.state.Φ[:, :] .= g * 1.0 * exp.(- ((gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

# output initial state
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)
#println(integrateT())


# 0.5 * sin.((model.env.gi.c_lon .- deg2rad(45))) .* cos.(model.env.gi.c_lat)

for step=1:run_days
    println("Run day ", step)

    #=
    if step <= 10
        model.state.τx[:, :] .= a * exp(- (step-.5)/ 10.0)
        model.state.τy[:, :] .= b * exp(- (step-.5)/ 10.0)
    else
        model.state.τx[:, :] .= 0
        model.state.τy[:, :] .= 0
    end
    =#


    @time for substep = 1:Int64(86400/Δt)
        Dyn.stepModel!(model)
        RecordTool.record!(recorder)
        if mod(substep,2)==0
            RecordTool.avgAndOutput!(recorder)
        end

    end
    println("Avg ϕ: ", avg(model.state.Φ.^2.0))
#    println("Avg KE: ", integrateKE())
end

