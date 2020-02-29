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

PTI = ProgramTunnelInfo(
    reverse_role  = true,
    recv_channels = 2,
)

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

overall_run_time = 0.0
overall_run_N    = 0

println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
end_time = Base.time()
while true

    global OMDATA, stage, t_cnt, output_filename, end_time

    new_end_time = Base.time()
    exe_time = new_end_time - end_time
    end_time = new_end_time

    println(format("# Current real time:      : {:s}", Dates.format(now(), "yyyy/mm/dd HH:MM:SS")))
    println(format("# Time Counter for RUN    : {:d}", t_cnt))
    println(format("# Total execution time    : {:d} s", floor(end_time - beg_time)))
    println(format("# Average execution time  : {:.2f} s", floor(end_time - beg_time) / t_cnt))
    println(format("# Last run execution time : {:.2f} s", exe_time))
    println(format("# Stage                   : {}", String(stage)))

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
     
        rm(configs[:archive_list], force=true)

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
                Δt            = Δt,
                write_restart = msg["WRITE_RESTART"] == "TRUE",
            )

        end

        global ocn_run_time += cost
        global ocn_run_N += 1

        println(format("*** Ocean takes {:.2f} secs. (Avg: {:.2f} secs) ***", cost, ocn_run_time / ocn_run_N))
       
        copy_list_to!(send_data_list_shared, send_data_list)
        sendData(PTI, "OK", send_data_list)

        t_cnt += 1

    elseif stage == :RUN && msg["MSG"] == "END"

        # move short_term_archive_files to long term archive directory
        println("===== Archiving files BEGIN =====")
        sdir = configs[:caseroot]
        ldir = configs[:archive_root]
        mkpath(ldir)
        for line in eachline(configs[:archive_list])

            args = split(line, ",")

            if length(args) == 0
                continue
            end
          
            action = args[1]
            args = args[2:end]

            if action in ["mv", "cp"]

                fname, src_dir, dst_dir = args

                if ! isdir(dst_dir)
                    mkpath(dst_dir)
                end
     
                src_file = joinpath(src_dir, fname)
                dst_file = joinpath(dst_dir, fname)

                if isfile(src_file)

                    if action == "mv"
                        mv(src_file, dst_file, force=true)
                        println(format("Moving file: {:s} ( {:s} => {:s} )", fname, src_dir, dst_dir))
                    elseif action == "cp"
                        cp(src_file, dst_file, force=true)
                        println(format("Copying file: {:s} ( {:s} => {:s} )", fname, src_dir, dst_dir))
                    end

                else
                    println("File does not exist: ", src_file)
                end

            elseif action == "rm"
                fname, fdir = args
                rm(joinpath(fdir, fname), force=true)
                println(format("Removing file: {:s} in {:s}", fname, fdir))
            else
                throw(ErrorException(format("Unknown action in archive list: {:s}", action)))
            end

        end

        println("===== Archiving files END =====")

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


