println("===== Universal Driver Initialization BEGIN =====")

if isdir(configs[:caserun])
    cd(configs[:caserun])
else
    throw(ErrorException("Caserun directory [ " * configs[:caserun] * " ] does not exist."))
end


function parseCESMTIME!(ts::AbstractString, timeinfo::AbstractArray{Integer})

    timeinfo[1] = parse(Int, ts[1:4])
    timeinfo[2] = parse(Int, ts[5:6])
    timeinfo[3] = parse(Int, ts[7:8])
    timeinfo[4] = parse(Int, ts[10:17])
end

function copy_list_to!(
    a,
    b,
)
    for i = 1:length(a)
        b[i][:] = a[i]
    end
end

output_vars = Dict()

stage = :INIT

# Need to find a way to avoid receiving the msg
# from last run

PTI = ProgramTunnelInfo(
    reverse_role  = true,
    recv_channels = 2,
)

#=
PTI = ProgramTunnelInfo(
    path             = configs[:tmp_folder],
    timeout          = configs[:timeout],
    buffer           = 0.1,
    recv_first_sleep = 0.1,
    reverseRole      = true,
)
=#

map = NetCDFIO.MapInfo{Float64}(configs[:domain_file])

t_cnt = 1
output_filename = ""
nullbin  = [zeros(Float64, 1)]

timeinfo = zeros(Integer, 4) 
timeinfo_old = copy(timeinfo)
timeinfo_old .= -1
t_flags = Dict()

ocn_run_time = 0.0
ocn_run_N    = 0

println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
while true

    global OMDATA, stage, t_cnt, output_filename

    end_time = Base.time()

    println(format("Execution time          : {:d} s", floor(end_time - beg_time)))
    println(format("# Time Counter for RUN  : {:d}", t_cnt))
    println(format("# Stage                 : {}", String(stage)))

    msg = parseMsg( recvData!(PTI, nullbin, which=1) )

    println("==== MESSAGE RECEIVED ====")
    print(json(msg, 4))
    println("==========================")

    if msg["MSG"] in ["INIT", "RUN"]
        parseCESMTIME!(msg["CESMTIME"], timeinfo)
    end

    if stage == :INIT && msg["MSG"] == "INIT"

        println("===== INITIALIZING MODEL: ", OMMODULE.name , " =====")

        OMDATA = OMMODULE.init(
            casename     = configs[:casename],
            map          = map,
            t            = timeinfo,
            configs      = configs,
            read_restart = (msg["READ_RESTART"] == "TRUE") ? true : false,
        )
 
        # Must declare SharedArray to avoid implicit copying!!!!
        global send_data_list_shared = SharedArray{Float64}[OMDATA.o2x["SST"], OMDATA.o2x["QFLX2ATM"]]
        global recv_data_list_shared = SharedArray{Float64}[]
 
        global send_data_list = Array{Float64}[OMDATA.o2x["SST"], OMDATA.o2x["QFLX2ATM"]]
        global recv_data_list = Array{Float64}[]
     
        rm(configs[:short_term_archive_list], force=true)

        global x2o_available_varnames  = split(msg["VAR2D"], ",")
        global x2o_wanted_varnames = keys(OMDATA.x2o)
        global x2o_wanted_flag     = [(x2o_available_varnames[i] in x2o_wanted_varnames) for i = 1:length(x2o_available_varnames)]

        println("List of available x2o variables:")
        for (i, varname) in enumerate(x2o_available_varnames)
            push!(recv_data_list_shared, ( x2o_wanted_flag[i] ) ? OMDATA.x2o[varname] :  SharedArray{Float64}(map.lsize))
            push!(recv_data_list       , ( x2o_wanted_flag[i] ) ? OMDATA.x2o[varname] :  zeros(Float64, map.lsize))
            println(format(" ({:d}) {:s} => {:s}", i, varname, ( x2o_wanted_flag[i] ) ? "Wanted" : "Abandoned" ))
        end

        
        copy_list_to!(send_data_list_shared, send_data_list)

        sendData(PTI, "OK", send_data_list)
        sendData(PTI, "MASKDATA", [map.mask])

        stage = :RUN
        
    elseif stage == :RUN && msg["MSG"] == "RUN"

        t_flags[:new_year]  = (timeinfo[1] != timeinfo_old[1])
        t_flags[:new_month] = (timeinfo[2] != timeinfo_old[2])
        t_flags[:new_day]   = (timeinfo[3] != timeinfo_old[3])

        timeinfo_old[:] = timeinfo

        recvData!(
            PTI,
            recv_data_list,
            which=2
        )
      
        copy_list_to!(recv_data_list, recv_data_list_shared)
        
        println("Calling ", OMMODULE.name, " to do MAGICAL calculations")

        Δt = parse(Float64, msg["DT"])

        cost = @elapsed let

            OMMODULE.run(
                OMDATA;
                t             = timeinfo,
                t_cnt         = t_cnt,
                t_flags       = t_flags,
                Δt            = Δt
                write_restart = msg["WRITE_RESTART"] == "TRUE",
            )

        end

        global ocn_run_time += cost
        global ocn_run_N += 1

        println(format("*** It takes {:.2f} secs. (Avg: {:.2f} secs) ***", cost, ocn_run_time / ocn_run_N))
       
        copy_list_to!(send_data_list_shared, send_data_list)
        sendData(PTI, "OK", send_data_list)

        t_cnt += 1

    elseif stage == :RUN && msg["MSG"] == "END"

        # move short_term_archive_files to long term archive directory
        if configs[:enable_long_term_archive]
            println("===== Long term archiving files BEGIN =====")
            sdir = configs[:short_term_archive_dir]
            ldir = configs[:long_term_archive_dir]
            mkpath(ldir)
            for fname in eachline(configs[:short_term_archive_list])
                src = joinpath(sdir, fname)
                dst = joinpath(ldir, fname)
                
                if isfile(src)
                    mv(src, dst, force=true)
                    println("Long term archiving file: ", fname)
                else
                    println("File does not exist: ", fname)
                end
            end

            println("===== Long term archiving files END =====")
        end

        OMMODULE.final(OMDATA) 
        
        # Print some report... ?
        
        println("Simulation ends peacefully.")
        break
    else
        #OMMODULE.crash(OMDATA) 
        sendData(PTI, "CRASH", send_data_list)
        throw(ErrorException("Unknown status: stage " * string(stage) * ", MSG: " * string(msg["MSG"])))
    end

    flush(stdout)
end


