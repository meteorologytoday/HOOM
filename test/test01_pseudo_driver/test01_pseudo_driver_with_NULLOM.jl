include("HOOM/src/driver/pseudo_driver.jl")

include("CORE_NULLOM.jl")

using CFTime

t_start = DateTimeNoLeap(1, 1, 1)
Δt = Dates.Second(1800)
steps = Int64(2*86400 / Δt.value)
read_restart = false

configs = Dict(
    :casename => "mini_model_NULLOM"
)

runModel(
    CORE_NULLOM, 
    t_start,
    Δt,
    steps,
    read_restart,
    configs,
)
