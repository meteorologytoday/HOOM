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
    t_end   = DateTimeNoLeap(1, 2, 1)
    read_restart = false
    Δt = Second(86400)

    config = Dict{Any, Any}(

        :DRIVER => Dict(
            :casename           => "Sandbox",
            :caseroot           => joinpath(@__DIR__, "Sandbox", "caseroot"),
            :caserun            => joinpath(@__DIR__, "Sandbox", "caserun"),
            :archive_root       => joinpath(@__DIR__, "Sandbox", "archive"),
        ),

        :MODEL_MISC => Dict(
            :timetype               => "DateTimeNoLeap",
            :init_file              => "",#nothing,#joinpath(@__DIR__, "ocn_init.nc"),
            :rpointer_file          => "rpointer.hoom",
            :daily_record           => [:ALL,],
            :monthly_record         => [:ALL,],
            :enable_archive         => true,
        ),

        :MODEL_CORE => Dict(
            #:domain_file                  => "domain.ocn.gx1v6.090206.nc",
            :domain_file                  => joinpath(@__DIR__, "domain.ocn_aqua.fv4x5_gx3v7.091218.nc"),
            :cdata_file                   => joinpath(@__DIR__, "forcing.nc"),

            :cdata_beg_time               => DateTimeNoLeap(1, 1, 1, 0, 0, 0),
            :cdata_end_time               => DateTimeNoLeap(2, 1, 1, 0, 0, 0),
            :cdata_align_time             => DateTimeNoLeap(1, 1, 1, 0, 0, 0),

            :z_w               => collect(Float64, 0:-10:-350),

            :substeps           => 8,
            :MLD_scheme                   => :static,
            :Qflx                         => :off,
            :Qflx_finding                 => :off,
            :weak_restoring               => :on,
            :convective_adjustment        => :off,
            :advection_scheme             => :ekman_codron2012_partition,

            :τwk_TEMP => 86400.0 * 30.0,
            :τwk_SALT => 86400.0 * 30.0,


            :Ekman_layers      => 5,
            :Returnflow_layers => 25,
        ),

    )

    gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
        config[:MODEL_CORE][:domain_file],
        R  = 6371229.0,
        Ω  = 2π / (86400 / (1 + 365/365)),
    )




end


coupler_funcs = (

    before_model_init! = function()

        global comm, rank
        global t_start, read_restart, config
         
        is_master = (rank == 0)
        println("My rank: ", rank) 
        writeLog("[Coupler] before model init.")
        writeLog("Reading configuration.")


        if ! is_master
            t_start = nothing 
            read_restart = nothing
            config = nothing
        end
        
        t_start = MPI.bcast(t_start, 0, comm) 
        read_restart = MPI.bcast(read_restart, 0, comm) 
        config = MPI.bcast(config, 0, comm) 

        return t_start, read_restart, config
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

            #OMDATA.x2o["SWFLX"][1, :, :] .= - (cos.(deg2rad.(lat).+0.1) .+ 1) / 2 .* (sin.(deg2rad.(lon)) .+ 1)/2 * 100.0
            #OMDATA.x2o["SWFLX"][1, :, :] .= + 200.0
            #OMDATA.x2o["TAUX_east"][1, :, :]   .= 0.2 * (cos.(deg2rad.(lat).+0.1) .+ 1) / 2
            #OMDATA.x2o["TAUY_north"][1, :, :]  .= 0.1 * (sin.(deg2rad.(lon)) .+ 1) / 2
#cos.(deg2rad.(lat))
            
        end


        return return_values

    end,
    after_model_run! = function(OMMODULE, OMDATA)
        writeLog("[Coupler] After model run")
    end,
    final = function(OMMODULE, OMDATA)
        writeLog("[Coupler] Finalize")
    end 
)



runModel(
    CESMCORE_HOOM, 
    coupler_funcs,
)
