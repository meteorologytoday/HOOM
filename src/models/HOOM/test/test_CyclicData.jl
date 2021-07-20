include("../CyclicData.jl")

using .CyclicData
using Plots

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

t_needed = range(0, 365, length=366) |> collect
data_needed = t_needed |> copy

for (i, t) in enumerate(t_needed)

    CyclicData.getData!(cdm, t, varnames, data)
    data_needed[i] = data_shaped[217, 15]

end


scatter(t_needed, data_needed)
