include("HOOM/src/driver/pseudo_driver_MPI.jl")
include("HOOM/src/models/HOOM_beta/CESMCORE_HOOM.jl")
include("HOOM/src/share/Log.jl")
include("HOOM/src/share/PolelikeCoordinate.jl")

using MPI
using CFTime
using Dates
using .PolelikeCoordinate



MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0

    t_start = DateTimeNoLeap(1, 1, 1)
    Δt = Second(86400)
#    steps = 31
    t_end = DateTimeNoLeap(1, 2, 1)
    #t_end = t_start + Δt * steps
    read_restart = false

    configs = Dict(
        :substeps    => 8, # This controls how many steps will occur for each CESM coupling. Example: ocean couple to atmosphere every 24 hours but itself steps every 3 hours. This means we would expect `Δt` = 86400, and we set `substeps` = 8.
        :enable_archive     => true,
        :archive_list       => "archive_list.txt",
        :rpointer_file      => "rpointer.hoom",

        :casename                     => "Sandbox",
        :caseroot                     => joinpath(@__DIR__, "Sandbox", "caseroot"),
        :caserun                      => joinpath(@__DIR__, "Sandbox", "caserun"),
        :domain_file                  => "domain.ocn_aqua.fv4x5_gx3v7.091218.nc",
        #:domain_file                  => "domain.ocn.gx1v6.090206.nc",
        :archive_root                 => joinpath(@__DIR__, "Sandbox", "hist"),
        :enable_archive               => true,
        :daily_record                 => :ALL,
        :monthly_record               => :ALL,
        :yearly_snapshot              => true,
        :substeps                     => 8,
        :init_file                    => nothing,#joinpath(@__DIR__, "ocn_init.nc"),
        
        :MLD_scheme                   => :datastream,
        :Qflux_scheme                 => :off,
        :Qflux_finding                => :off,
        :seaice_nudging               => :off,
        :vertical_diffusion_scheme    => :off,
        :horizontal_diffusion_scheme  => :off,
        :relaxation_scheme            => :off,
        :convective_adjustment_scheme => :off,
        :radiation_scheme             => :exponential,
        :advection_scheme             => :static,#ekman_codron2012_partition,
    )

    gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
        configs[:domain_file],
        R  = 6371229.0,
        Ω  = 2π / (86400 / (1 + 365/365)),
    )




end


coupler_funcs = (

    before_model_init! = function()

        global comm, rank
        global t_start, read_restart, configs
         
        is_master = (rank == 0)
        println("My rank: ", rank) 
        writeLog("[Coupler] before model init.")
        writeLog("Reading configuration.")


        if ! is_master
            t_start = nothing 
            read_restart = nothing
            configs = nothing
        end
        
        t_start = MPI.bcast(t_start, 0, comm) 
        read_restart = MPI.bcast(read_restart, 0, comm) 
        configs = MPI.bcast(configs, 0, comm) 

        return t_start, read_restart, configs
    end,

    after_model_init! = function(OMMODULE, OMDATA)
        writeLog("[Coupler] After model init")
    end,

    before_model_run! = function(OMMODULE, OMDATA)
        writeLog("[Coupler] Before model run")
        writeLog("[Coupler] This is where flux exchange happens.")

        global comm, rank
        
        is_master = rank == 0
        
        # Test if t_end is reached
        return_values = nothing
        if is_master
            t_end_reached = OMDATA.clock.time >= t_end
           
            if t_end_reached
                return_values = ( :END, Δt, t_end_reached )
            else

                return_values = ( :RUN, Δt, t_end_reached )
            end
        end

        return_values = MPI.bcast(return_values, 0, comm)

        # Deal with coupling
        if is_master

            writeLog("[Coupler] Need to broadcast forcing fields.")
            # compute flux
            lat = gf.yc
            lon = gf.xc

            #OMDATA.x2o["SWFLX"][1, :, :] .= + (cos.(deg2rad.(lat)) .+ 1) / 2 .* (sin.(deg2rad.(lon)) .+ 1)/2 * 100.0
            #OMDATA.x2o["SWFLX"][1, :, :] .= + 200.0
            OMDATA.x2o["TAUX_east"][1, :, :]   .= 1.0
            OMDATA.x2o["TAUY_north"][1, :, :]  .= 0.0
#cos.(deg2rad.(lat))
            
        end


        return return_values

    end,
    after_model_run! = function(OMMODULE, OMDATA)
        writeLog("[Coupler] After model run")
    end,
    finalize! = function(OMMODULE, OMDATA)
        writeLog("[Coupler] Finalize")
    end 
)



runModel(
    CESMCORE_HOOM, 
    coupler_funcs,
)
