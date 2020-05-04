using ArgParse
using JSON

include("BLOHSOM.jl")
include("../../share/constants.jl")

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

    run_days = 300
    output_file = "output.nc"

end

Δt = 86400.0 / 24

Nz_f = 20

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx3v7.120323.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx3v7.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv4x5_gx3v7.091218.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv4x5.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc"

mi = ModelMap.MapInfo{Float64}(hrgrid_file)
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
)

cutoff = 20

mi.mask .= 1.0
mi.mask[:, 1:cutoff] .= 0 
mi.mask[:, end-cutoff+1:end] .= 0 

topo = similar(mi.mask)
topo .= -2000

ocn_env = BLOHSOM.OcnEnv(
    hrgrid_file,
    topo_file,
    "topo",
    Δt,
    12,
    8,
    collect(Float64, range(0.0, -100.0, length=21)),
    [1, 1, 1, 1, 1, 2, 2, 2, 9],
    0,
    1000.0,
    1e-3,
    1e-3,
    [1e-3, 1e-3], 
    [1e-3, 1e-3], 
    0.48,
    23.0,
    1e-2,
    [10.0, 1000.0],
    ["T.nc", "S.nc"],
    ["T", "S"],
    [NaN, NaN],
    :dynamic,
    :exponential_decay,
    true,
    false,
    false;
    mask2 = mi.mask,
    topo  = topo,
)

println("##### Initialize model #####")
model = BLOHSOM.init!(ocn_env)

du = model.shared_data.data_units
#du[:Φ].data .= 0.01 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

σ = 750e3
#du[:SWFLX].data .= -1000.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
#du[:X].odata[30:40, 15:25, 1:5, 1] .+= 10.0
#du[:X_ML].odata[30:40, 15:25, 1] .+= 10.0

@sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
    @spawnat pid let
        
        BLOHSOM.tmd_slave.model.state.h_ML[:, :] .= 10.0
        BLOHSOM.tmd_slave.model.state.X_ML[:, :, :, 1] .= 10.0
        BLOHSOM.tmd_slave.model.state.X[:, :, :, 1] .= 10.0

        BLOHSOM.tmd_slave.model.state.X[1, :, :, 1] .= 10 .+ 20.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)
        BLOHSOM.tmd_slave.model.state.X_ML[:, :, 1]   .= BLOHSOM.tmd_slave.model.state.X[1, :, :, 1]
        BLOHSOM.Tmd.initialization!(BLOHSOM.tmd_slave.model)
    end
end

@sync BLOHSOM.touchTmd!(model, :END_TMD2MAS, :S2M)

recorder = BLOHSOM.getBasicRecorder(model)
RecordTool.setNewNCFile!(recorder, output_file)
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)


@time for step=1:run_days
    println("##### Run day ", step)

    if step >= 20 || true
        du[:SWFLX].data .= 0.0
    end
    @time BLOHSOM.stepModel!(model, false)
    RecordTool.record!(recorder)
    RecordTool.avgAndOutput!(recorder)
end
