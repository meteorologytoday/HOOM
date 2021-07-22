include("HOOM/src/driver/pseudo_driver.jl")

using CFTime

t_start = DateTimeNoLeap(1, 1, 1)
Δt = Dates.Second(1800)
steps = Int64(2*86400 / Δt.value)

runOfflineModel(t_start, Δt, steps, nothing, nothing)
