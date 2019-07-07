program_beg_time = Base.time()

include("config.jl")
using Printf
using Formatting
using NCDatasets
import Statistics: mean, std


if length(ARGS) != 1 
    throw(ErrorException("Length of ARGS must be 1. That is the longitude index."))
end

target_i = parse(Int, ARGS[1])

# construct tmp folder
tmp_dir = joinpath(main_dir, "stan_tmp", format("{:03d}", target_i))
mkpath(tmp_dir)


sub_output_N = ceil(Integer, nlat / sub_output_size)
output_filenames = [
    normpath(joinpath(main_dir, format("{:03d}_{:03d}.jld", target_i, i))) for i = 1:sub_output_N
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

script_path = normpath(joinpath(dirname(@__FILE__), "..", "..", "STAN_code", "NKMLM.stan"))
model_script = read(script_path, String)

stanmodel = Stanmodel(
    name="STAN",
    nchains=nchains,
    num_samples=num_samples,
    num_warmup=num_warmup,
    model=model_script,
    pdir=tmp_dir,
)

#display(stanmodel)

h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q.{}", i) for i = 1:12 ]

data = Dict()
time_stat = Dict()

total_time = 0.0

β_mean = zeros(Float64, sub_output_size, 24)
β_std = zeros(Float64, sub_output_size, 24)


Dataset(F_filename, "r") do ds
    global F = convert(Array{Float64}, nomissing(ds["SHF"][target_i, :, :], NaN))
end

Dataset(SST_filename, "r") do ds
    global SST = convert(Array{Float64}, nomissing(ds["SST"][target_i, :, :], NaN))
end

Dataset(MLD_filename, "r") do ds
    global MLD = convert(Array{Float64}, nomissing(ds["omlmax"][target_i, :, :], NaN))
end


init_h = zeros(Float64, 12) .+ 30.0

N  = size(F)[end]
dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

for i = 1:sub_output_N

    if file_exists[i]
        continue
    end

    β_mean .= NaN
    β_std  .= NaN

    β_mean .= i
    β_std  .= i

    lat_rng = 1 + sub_output_size * (i-1) : min(sub_output_size * i, nlat)

    filename = output_filenames[i]

    for j = lat_rng

        continue 
        
        println(format("Doing (lat, lon) = ({:d} / {:d}, {:d} / {:d})", j, nlat, target_i, nlon))

        if isfinite(SST[j, 1])

            beg_time = Base.time()
            data = Dict(
                "raw_N"  => N, 
                "dom"    => dom,
                "steps"  => steps,
                "h_max"  => maximum(MLD[j, :]) + 100.0,
                "raw_T"  => SST[j, :], 
                "raw_F"  => F[j, :],
                "T_std"  => σ_T,
                "c_sw"   => c_p,
                "rho_sw" => ρ,
            )

            init = Dict(
                "h" => init_h
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
            data_Td = sim1[:, ["Td"], :].value

            for i = 1:12
                h_mean[i] = mean(data_h[:, i, :])
                h_std[i]  = std(data_h[:, i, :])

                Q_mean[i] = mean(data_Q[:, i, :])
                Q_std[i]  = std(data_Q[:, i, :])
            end

            Td_mean = mean(data_Td)
            Td_std  = std(data_Td)

            β_mean[j,  1:12] = h_mean
            β_mean[j, 13:24] = Q_mean
            β_mean[j,    25] = Td_mean

            β_std[j,  1:12] = h_std
            β_std[j, 13:24] = Q_std
            β_std[j,    25] = Td_std

            println("h_mean", h_mean)
            println("Q_mean", Q_mean)
            println("Td_mean", Td_mean)


            time_stat = Base.time() - beg_time

            global total_time += time_stat
            println("##########")
            println(format("Done (lat, lon) = ({:d} / {:d}, {:d} / {:d})", j, nlat, target_i, nlon))
            println(format("Stan fit: {:.2f} min. Total time: {:.2f} min. ", time_stat / 60.0, total_time / 60.0 ))
            println("##########")

        end


    end

    using JLD
    println("# Output filename: ", filename)
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
    nchains,
    num_samples
)
println("##########")
