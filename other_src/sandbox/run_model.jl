include("./HOOM/src/models/HOOM/HOOM.jl")
include("./HOOM/RecordTool.jl")

using .HOOM
using .RecordTool

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
     
            "--init-file"
                help = "Output file."
                arg_type = String
                required = true
     
            "--domain-file"
                help = "Domain file."
                arg_type = String
                required = true
     
            "--output-file"
                help = "Output file."
                arg_type = String
                required = true
     
        end

        return parse_args(ARGS, s)
    end

    parsed = parse_commandline()

else
    parsed = Dict(
        "run-days"    => 50,
        "init-file"   => "gen_data/init_aqp.nc",
        "output-file" => "gen_data/record_revise.nc",
        "domain-file" => "raw_data/fv45_remove_lnd.nc",
    )
end
print(json(parsed, 4))


print("Loading initial file: ", parsed["init-file"], " ...")
ocn = HOOM.loadSnapshot(
    parsed["init-file"];
    gridinfo_file = parsed["domain-file"] 
)


ocn.Dh_T = ocn.Dv_T
#ocn.Ts[1, :, :] .= 10.0
ocn.Ts[1, :, :] .+= 10.0 * cos.(ocn.gi.c_lat) .* cos.(ocn.gi.c_lon)

#ocn.u[1:2,  :, :]   .= repeat( reshape( 1.0 * cos.(ocn.gi.c_lat) .* sin.(ocn.gi.c_lon), 1, ocn.Nx, ocn.Ny), outer=( 2, 1, 1))
#ocn.u[3:12,  :, :]  .= repeat( reshape(-1.0 * cos.(ocn.gi.c_lat) .* sin.(ocn.gi.c_lon), 1, ocn.Nx, ocn.Ny), outer=(10, 1, 1))

ocn.u[1:2,  :, :]   .= 1.0
ocn.u .= 1.0


ocn.T_ML[:, :]  .= ocn.Ts[1, :, :]

ocn.Ss   .= 35.0
ocn.S_ML .= 35.0


HOOM.init(ocn)
println("done.")


# Setting up recorder
complete_varlist = HOOM.getCompleteVariableList(ocn)
varlist = []

for varname in keys(complete_varlist)
    println(format("Using varaible: {:s}", varname))
    push!(varlist, ( varname, complete_varlist[varname]... ) )
end

recorder = RecordTool.Recorder(
    Dict(
        "Nx_ext"  => ocn.Nx+1,
        "Ny_ext"  => ocn.Ny+1,
        "Nx"      => ocn.Nx,
        "Ny"      => ocn.Ny,
        "Nz_bone" => ocn.Nz_bone,
        "zs_bone" => length(ocn.zs_bone),
        "NP_zs_bone" => ocn.Nz_bone+1,
    ), varlist, Dict()
)
 
RecordTool.setNewNCFile!(recorder, parsed["output-file"])





SECS_PER_DAY = 86400.0
DAYS_PER_MON = 30
MONS_PER_YEAR= 12
DAYS_PER_YEAR = DAYS_PER_MON * MONS_PER_YEAR
SECS_PER_YEAR = DAYS_PER_YEAR * SECS_PER_DAY

TOTAL_YEARS = 1

TOTAL_DAYS = parsed["run-days"]
TOTAL_SECS = TOTAL_DAYS * SECS_PER_DAY

t_sim = collect(Float64, range(0.0, step=SECS_PER_DAY, stop=TOTAL_SECS))[1:end-1]
total_steps = length(t_sim)
Δt = t_sim[2] - t_sim[1]

Δt == SECS_PER_DAY || throw(ErrorException("Δt does not equal to SECS_PER_DAY"))

for step = 1:total_steps-1

    print(format("\rStep : {:d} / {:d}", step, length(t_sim)))

    RecordTool.record!(recorder)
    RecordTool.avgAndOutput!(recorder)

    #=
    ocn.in_flds.swflx  .=  153.0
    ocn.in_flds.nswflx .=  25.0
    #ocn.in_flds.frwflx .=  -7e-5
    ocn.in_flds.taux   .=  0.0
    ocn.in_flds.tauy   .=  0.0
    ocn.in_flds.taux .= .1 * sin.((ocn.gi.c_lon .+ step * π/180.0)) .* sin.(ocn.gi.c_lat)
    ocn.in_flds.tauy .= .1 * cos.((ocn.gi.c_lon .+ step * π/180.0) * 2) .* sin.(ocn.gi.c_lat)

    ocn.in_flds.Tclim .= 20
    ocn.in_flds.Sclim .= 35

    =#
    ocn.in_flds.h_ML   .=  50.0#=rand(ocn.Nx, ocn.Ny) * 20 .+ =# #step .+50.0;


    #if step < 20
    #    ocn.in_flds.swflx .= -100 * sin.(ocn.gi.c_lon * 5.0) .* cos.(ocn.gi.c_lat * 3.0)
    #    ocn.in_flds.frwflx .= -0.01/3600.0
    #end

    # HOOM
    #= 
    HOOM.run!(
        ocn;
        substeps      = 12,
        qflx_scheme   = :energy_flux,
        use_h_ML      = true,
        Δt            = Δt,
        do_vert_diff  = false,
        do_horz_diff  = false,
        do_relaxation = true,
        do_convadjust = false,
        rad_scheme    = :exponential,
        adv_scheme    = :static,
    )
    =#
    # EntOM / SOM
    HOOM.run!(
        ocn;
        substeps      = 12,
        do_qflx       = false,
        do_qflx_finding  = false,
        use_h_ML      = true,
        Δt            = Δt,
        do_vert_diff  = false,
        do_horz_diff  = false,
        do_relaxation = false,
        do_convadjust = false,
        do_seaice_nudging = false,
        rad_scheme    = :step,
#        adv_scheme    = :static,
#        adv_scheme    = :ekman_codron2012_partition,
        adv_scheme    = :test,
    )


end

println()
println("Simulation done.")
