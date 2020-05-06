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

    run_days = 365
    output_file = "output.nc"

end


coupling = 4
Δt_day = 86400.0

Nz_f = 20

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx3v7.120323.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx3v7.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv4x5_gx3v7.091218.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv4x5.nc"

hrgrid_file = "/seley/tienyiah/CESM_domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc"
#=

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
=#
Ly = 100e3 * 60.0
gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = 240,
        Ny   = 120,
        Ly   = Ly,
        lat0 = 0.0 |> deg2rad,
        β    = 2*Ωe / Re,
)

gi = PolelikeCoordinate.genGridInfo(gf);

xcutoff = 25
ycutoff = 25

gf.mask                       .= 1
gf.mask[:, 1:ycutoff]         .= 0 
gf.mask[:, end-ycutoff+1:end] .= 0 

gf.mask[1:xcutoff, :]         .= 0 
gf.mask[end-xcutoff+1:end, :] .= 0 


topo = similar(gf.mask)
topo .= -4000
z_bnd_f = collect(Float64, range(0.0, -500.0, length=21))
#push!(z_bnd_f, -4000)
ocn_env = BLOHSOM.OcnEnv(
    hrgrid                = gf,
    topo_file             = topo_file,
    topo_varname          = "topo",
    Δt                    = Δt_day / coupling,
    substeps_dyn          = 3,
    substeps_tmd          = 3,
    z_bnd_f               = z_bnd_f,
    #height_level_counts   = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5],
    height_level_counts   = [2, 2, 2, 2, 2, 2, 8],
    NX_passive            = 0,
    deep_threshold        = 1000.0,
    Kh_m                  = 30000.0,
    Kv_m                  = 1e-3,
    Kh_X                  = [1000.0, 1e-5], 
    Kv_X                  = [1e-5, 1e-5], 
    R                     = 0.48,
    ζ                     = 23.0,
    we_max                = 1e-2,
    MLT_rng               = [10.0, 1000.0],
    X_wr_file             = ["T.nc", "S.nc"],
    X_wr_varname          = ["T", "S"],
    t_X_wr                = [NaN, NaN],
    MLT_scheme            = :prescribe,
    radiation_scheme      = :exponential_decay,
    convective_adjustment = true,
    use_Qflux             = false,
    finding_Qflux         = false;
    mask2                 = gf.mask,
    topo                  = topo,
)

println("##### Initialize model #####")
model = BLOHSOM.init!(ocn_env)

du = model.shared_data.data_units
#du[:Φ].data .= 0.01 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

σz = 50.0
σ = 750e3
σ_warmpool = 750e3

#du[:X].odata[30:40, 15:25, 1:5, 1] .+= 10.0
#du[:X_ML].odata[30:40, 15:25, 1] .+= 10.0



z_mid = (z_bnd_f[1:end-1] + z_bnd_f[2:end]) / 2
z_mid = repeat(reshape(z_mid, :, 1, 1), outer=(1, gf.Nx, gf.Ny))
mask3 = repeat(reshape(gf.mask, 1, gf.Nx, gf.Ny), outer=(length(z_bnd_f)-1, 1, 1))

basic_T   = z_mid * 0
anomaly_T = z_mid * 0
    
@. basic_T = 10 #+ 15 * exp(z_mid / σz)

#warmpool_bnd_z = - 300.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
#warmpool_bnd_z = - 300.0 * exp.(- ( (gi.c_y .- ( Ly/4.0)).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
warmpool_bnd_z = - 130.0 * exp.(- ( (gi.c_y ).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
warmpool_bnd_z = repeat(reshape(warmpool_bnd_z, 1, size(warmpool_bnd_z)...), outer=(size(z_mid)[1], 1, 1))
anomaly_T[z_mid .> warmpool_bnd_z] .= 20.0

total_T = basic_T + anomaly_T
total_T[1:4, :, :] .= 30.0

h_ML = gf.mask * 0 .+ 10.0

total_T[mask3 .== 0.0] .= 0

@sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
    @spawnat pid let
        m = BLOHSOM.tmd_slave.model

        m.forcing.h_ML .= h_ML
        
        m.state.X_ML[:, :, 1] .= total_T[1, :, :]
        m.state.X[:, :, :, 1] .= total_T

        BLOHSOM.Tmd.initialization!(BLOHSOM.tmd_slave.model)

    end
end

@sync @spawnat model.job_dist_info.dyn_slave_pid let

#    BLOHSOM.dyn_slave.model.state.Φ .= 0.01 * exp.(- ( (gi.c_y ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

end
BLOHSOM.touchTmd!(model, :END_TMD2MAS, :S2M)
BLOHSOM.touchDyn!(model, :END_DYN2MAS, :S2M)

recorder = BLOHSOM.getBasicRecorder(model)
RecordTool.setNewNCFile!(recorder, output_file)
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)


@time for step=1:run_days
    println("##### Run day ", step)

#    du[:NSWFLX].data .= (du[:X_ML].odata[:, :, 1] .- 30.0) * 1026*3996*25 / (86400*10)

    @time for c = 1:coupling    
        BLOHSOM.stepModel!(model, false)
    end

    RecordTool.record!(recorder)
    RecordTool.avgAndOutput!(recorder)
end
