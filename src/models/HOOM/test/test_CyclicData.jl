include("../CyclicData.jl")

using .CyclicData
using Plots
using CFTime

varnames = ["vice_target"]

cdm = CyclicData.CyclicDataManager(
    filename     = "/seley/tienyiah/projects/paper2021/paper2021/simulation_shared_files/vice_target_file/forcing.vice.gx1v6.RCP85.nc" ,
    varname_time = "time",
    varnames     = varnames,
    beg_time     = 0.0,
    cyc_time     = 365.0,
)



data = Dict()
for varname in varnames
    data[varname] = zeros(Float64, size(cdm.data[varname])[1])
end

data_shaped = reshape(data["vice_target"], 320, 384)
raw_data_shaped = reshape(cdm.data["vice_target"], 320, 384, :)

t_needed = range(0, 365, length=366) |> collect
data_needed = t_needed |> copy

first_day_of_month = [ timeencode( DateTimeNoLeap(0, m, 1), "days since 0000-01-01 00:00:00", "noleap") for m=1:12]

for (i, t) in enumerate(t_needed)

    CyclicData.getData!(cdm, t, varnames, data)
    data_needed[i] = data_shaped[217, 15]

end

days_of_month = 
p = scatter(t_needed, data_needed, ticks=first_day_of_month)
scatter!(p, cdm.t_vec, raw_data_shaped[217, 15, :])

