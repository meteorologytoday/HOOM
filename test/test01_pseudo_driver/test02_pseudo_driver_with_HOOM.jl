include("HOOM/src/driver/pseudo_driver.jl")
include("HOOM/src/models/HOOM/CESMCORE_HOOM.jl")

using CFTime
using Dates

t_start = DateTimeNoLeap(1, 1, 1)
Δt = Second(86400)
steps = 59

configs = Dict(
    :substeps    => 8, # This controls how many steps will occur for each CESM coupling. Example: ocean couple to atmosphere every 24 hours but itself steps every 3 hours. This means we would expect `Δt` = 86400, and we set `substeps` = 8.
    :daily_record              => [],
    :monthly_record            => :ESSENTIAL,
    :enable_archive            => true,
    :archive_list              => "archive_list.txt",
    :rpointer_file             => "rpointer.hoom",

    :casename     => "Sandbox",
    :caseroot     => joinpath(@__DIR__, "Sandbox", "caseroot"),
    :caserun      => joinpath(@__DIR__, "Sandbox", "caserun"),
    :domain_file  => "domain.lnd.fv4x5_gx3v7.091218.nc",
    :archive_root => joinpath(@__DIR__, "Sandbox", "hist"),
    :enable_archive             => true,
    :daily_record               => :ESSENTIAL,
    :monthly_record             => :ESSENTIAL,
    :yearly_snapshot            => true,
    :substeps                   => 8,
    :init_file                  => joinpath(@__DIR__, "ocn_init.nc"),
    
    :MLD_scheme                   => :datastream,
    :Qflux_scheme                 => :off,
    :Qflux_finding                => :off,
    :seaice_nudging               => :off,
    :vertical_diffusion_scheme    => :off,
    :horizontal_diffusion_scheme  => :off,
    :relaxation_scheme            => :off,
    :convective_adjustment_scheme => :off,
    :radiation_scheme             => :step,
    :advection_scheme             => :static,#ekman_codron2012_partition,
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
