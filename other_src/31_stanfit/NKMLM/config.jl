using Formatting

# Check:
#
# model_name
# nchains, num_samples, num_warmup
# exp_name [contains initial condition status]
# 
# sbatch --job-name, --output, --cpus-per-task, --time, --array

include("setup_src_data.jl")

data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data")

nchains     = 2
num_samples = 1000
num_warmup  = 50


ρ    = 1026.0  # kg / m^3
c_p  = 3996.0  # J / kg / K
σ_T  = 1.0

sub_output_size = 10


model_name = "NCAR_LENS_g37"
exp_name = format("NKMLM_{}_c{:d}_s{:d}_w{:d}", model_name, nchains, num_samples, num_warmup)

main_dir = joinpath(data_path, exp_name)
mkpath(main_dir)
