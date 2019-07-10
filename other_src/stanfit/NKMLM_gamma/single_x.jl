program_beg_time = Base.time()

using Printf
using Formatting
using NCDatasets
import Statistics: mean, std

using ArgParse
using JSON

include("config.jl")

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--selected-index"
            help = "Selected fitted index of x-direction"
            arg_type = Int64
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))
print(json(config, 4))
    
ρ    = 1026.0  # kg / m^3
c_p  = 3996.0  # J / kg / K

Dataset(config["SST-file"], "r") do ds
    global nlon = ds.dim["nlon"]
    global nlat = ds.dim["nlat"]
end

target_i = parsed["selected-index"]


mkpath(config["main-dir"])
# construct tmp folder
tmp_dir = joinpath(config["main-dir"], "stan_tmp", format("{:03d}", target_i))
mkpath(tmp_dir)

total_sub_output = ceil(Integer, nlat / config["sub-output-size"])
output_filenames = [
    normpath(joinpath(config["main-dir"], format("{:03d}_{:03d}.jld", target_i, i))) for i = 1:total_sub_output
]

file_exists = [
    isfile(output_filenames[i]) for i = 1:length(output_filenames)
]

println(format("This program is going to fit {}/{} ", target_i, nlon))

if all(file_exists)
    println("Files are all present, so nothing more to do. End program now.")
    exit()
else
    println("Some files are missing, means we have work to do!")
end

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan
@printf("done\n")

script_path = normpath(joinpath(dirname(@__FILE__), "NKMLM_gamma.stan"))
model_script = read(script_path, String)

stanmodel = Stanmodel(
    name="STAN",
    nchains     = config["stan-chains"],
    num_samples = config["stan-samples"],
    num_warmup  = config["stan-warmup"],
    model       = model_script,
    pdir        = tmp_dir,
)

println("Model built.")

#display(stanmodel)

h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q.{}", i) for i = 1:12 ]

data = Dict()
time_stat = Dict()

total_time = 0.0

β_mean = zeros(Float64, config["sub-output-size"], 25)
β_std = zeros(Float64, config["sub-output-size"], 25)


Dataset(config["SHF-file"], "r") do ds
    global F = convert(Array{Float64}, nomissing(ds["SHF"][target_i, :, :], NaN))
end

Dataset(config["SST-file"], "r") do ds
    global SST = convert(Array{Float64}, nomissing(ds["SST"][target_i, :, 1, :], NaN))
end

Dataset(config["MLD-file"], "r") do ds
    global MLD = convert(Array{Float64}, nomissing(ds["HMXL"][target_i, :, :], NaN)) / 100.0  # raw data is in cm
end

N  = size(F)[end]
dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
steps = [1 for _ in 1:12]
for i = 1:total_sub_output

    if file_exists[i]

        println("File ", output_filenames[i], " already exists. Going to skip this.")
        
        continue
    end

    β_mean .= NaN
    β_std  .= NaN

    lat_rng = 1 + config["sub-output-size"] * (i-1) : min(config["sub-output-size"] * i, nlat)
    rng_offset = config["sub-output-size"] * (i-1)

    filename = output_filenames[i]

    for j = lat_rng

        println(format("Doing (lat, lon) = ({:d} / {:d}, {:d} / {:d})", j, nlat, target_i, nlon))

        if !isfinite(SST[j, 1])

            println("This grid point has NaN. Skip.")

        else
            
            mean_seasonal_h = mean( reshape(MLD[j, :], 12, :), dims=(2,) )[:, 1]
            mean_h = mean(mean_seasonal_h)

            init_partial_h = [mean_h for _ in 1:11]

            F_avg = mean( reshape(F[j, :], 12, :), dims=(2,) )[:, 1]

       
            println("mean_seasonal_h: ", mean_seasonal_h)
            println("mean_h: ", mean_h)
            println("F_avg: ", F_avg)
     
            beg_time = Base.time()
            data = Dict(
                "raw_N"  => N, 
                "dom"    => dom,
                "steps"  => steps,
                "h_mean" => mean_h,
                "raw_T"  => SST[j, :], 
                "raw_F"  => F[j, :],
                "T_std"  => config["T-sigma"],
                "c_sw"   => c_p,
                "rho_sw" => ρ,
            )

            init = Dict(
                "partial_h" => init_partial_h,
            )
            
            rc, sim1 = stan(
                stanmodel,
                [data];
                init = [init],
                CmdStanDir=ENV["CMDSTAN_HOME"]
            )

            if rc != 0
                println("There are errors!!")
                continue
            end
            
            println("Extracting result...")

            h_mean = zeros(12)
            h_std  = zeros(12)
            Q_mean = zeros(12)
            Q_std  = zeros(12)

            data_h = sim1[:, h_key,  :].value
            data_Q = sim1[:, Q_key,  :].value
            data_γ = sim1[:, ["gamma"], :].value

            for i = 1:12
                h_mean[i] = mean(data_h[:, i, :])
                h_std[i]  = std(data_h[:, i, :])
                Q_mean[i] = mean(data_Q[:, i, :])
                Q_std[i]  = std(data_Q[:, i, :])
            end

            println("Annual mean of MLD (old vs new):", mean_h, " vs ", mean(h_mean))

            γ_mean = mean(data_γ)
            γ_std  = std(data_γ)

            jj = j - rng_offset

            β_mean[jj,  1:12] = h_mean
            β_mean[jj, 13:24] = Q_mean
            β_mean[jj,    25] = γ_mean

            β_std[jj,  1:12] = h_std
            β_std[jj, 13:24] = Q_std
            β_std[jj,    25] = γ_std

            println("h_mean", h_mean)
            println("Q_mean", Q_mean)
            println("γ_mean", γ_mean)


            time_stat = Base.time() - beg_time

            global total_time += time_stat
            println("##########")
            println(format("Done (lat, lon) = ({:d} / {:d}, {:d} / {:d})", j, nlat, target_i, nlon))
            println(format("Stan fit: {:.2f} min. Total time: {:.2f} min. ", time_stat / 60.0, total_time / 60.0 ))
            println("##########")

        end


    end


    println("# Output filename: ", filename)
    using JLD
    save(
        filename,
        Dict(
            "β_mean" => β_mean[1:length(lat_rng), :],
            "β_std"  =>  β_std[1:length(lat_rng), :],
        )
    )

end

program_end_time = Base.time()

println("##########")
@printf("Total time used: %.2f min for %d points, with nchains = %d, num_samples = %d\n",
    (program_end_time-program_beg_time)/ 60.0,
    nlat,
    config["stan-chains"],
    config["stan-samples"],
)
println("##########")
