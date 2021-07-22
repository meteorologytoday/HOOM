include("HOOM/src/driver/pseudo_driver.jl")
include("HOOM/src/models/HOOM/CESMCORE_HOOM.jl")

using CFTime

t_start = DateTimeNoLeap(1, 1, 1)
Δt = Dates.Second(1800)
steps = Int64(2*86400 / Δt.value)

configs = Dict(
    :substeps    => 8, # This controls how many steps will occur for each CESM coupling. Example: ocean couple to atmosphere every 24 hours but itself steps every 3 hours. This means we would expect `Δt` = 86400, and we set `substeps` = 8.
    :daily_record              => [],
    :monthly_record            => :ESSENTIAL,
    :enable_archive            => true,
    :archive_list              => "archive_list.txt",
    :rpointer_file             => "rpointer.hoom",

    :casename     => "Sandbox",
    :caseroot     => @__DIR__,
    :caserun      => @__DIR__,
    :domain_file  => "domain.lnd.fv4x5_gx3v7.091218.nc",
    :archive_root => joinpath(@__DIR__, "hist"),
    :enable_archive             => true,
    :daily_record               => [],
    :monthly_record             => :ESSENTIAL,
    :yearly_snapshot            => true,
    :substeps                   => 8,
    :init_file                  => nothing,
    
    :MLD_scheme                   => :datastream,
    :Qflux_scheme                 => :on,
    :Qflux_finding                => :off,
    :seaice_nudging               => :off,
    :vertical_diffusion_scheme    => :off,
    :horizontal_diffusion_scheme  => :off,
    :relaxation_scheme            => :on,
    :convective_adjustment_scheme => :on,
    :radiation_scheme             => :step,
    :advection_scheme             => :ekman_codron2012_partition,
)

read_restart = false

runModel(
    CESMCORE_HOOM, 
    t_start,
    Δt,
    steps,
    read_restart,
    configs,
)
