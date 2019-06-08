println("===== Universal Driver Initialization BEGIN =====")

if isdir(configs["caserun"])
    cd(configs["caserun"])
else
    throw(ErrorException("Caserun directory [ " * configs["caserun"] * " ] does not exist."))
end


function parseCESMTIME!(ts::AbstractString, timeinfo::AbstractArray{Integer})

    timeinfo[1] = parse(Int, ts[1:4])
    timeinfo[2] = parse(Int, ts[5:6])
    timeinfo[3] = parse(Int, ts[7:8])
    timeinfo[4] = parse(Int, ts[10:17])
end


output_vars = Dict()

stage = :INIT

mkpath(configs["tmp_folder"])

PTI = ProgramTunnelInfo(path=configs["tmp_folder"], timeout=configs["timeout"])
reverseRole!(PTI)

map = NetCDFIO.MapInfo{Float64}(configs["domain_file"])

t_cnt = 1
output_filename = ""
buffer2d = zeros(UInt8, map.lsize * 8)
null2d   = zeros(Float64, map.lsize)

timeinfo = zeros(Integer, 4) 
timeinfo_old = copy(timeinfo)
timeinfo_old .= -1
t_flags = Dict()


println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
while true

    global OMDATA, stage, t_cnt, output_filename

    end_time = Base.time()

    println(format("Execution time          : {:d} s", floor(end_time - beg_time)))
    println(format("# Time Counter for RUN  : {:d}", t_cnt))
    println(format("# Stage                 : {}", String(stage)))

    msg = parseMsg(recvText(PTI))
    println("==== MESSAGE RECEIVED ====")
    print(json(msg, 4))
    println("==========================")

    if msg["MSG"] in ["INIT", "RUN"]
        parseCESMTIME!(msg["CESMTIME"], timeinfo)
    end

    if stage == :INIT && msg["MSG"] == "INIT"

        println("===== INITIALIZING MODEL: ", OMMODULE.name , " =====")

        OMDATA = OMMODULE.init(
            casename     = configs["casename"],
            map          = map,
            t            = timeinfo,
            configs      = configs,
            read_restart = (msg["READ_RESTART"] == "TRUE") ? true : false,
        )
        
        rm(configs["short_term_archive_list"], force=true)

        global x2o_available_varnames  = split(msg["VAR2D"], ",")
        global x2o_wanted_varnames = keys(OMDATA.x2o)
        global x2o_wanted_flag     = [(x2o_available_varnames[i] in x2o_wanted_varnames) for i = 1:length(x2o_available_varnames)]

        writeBinary!(joinpath(configs["tmp_folder"], "SST.bin"), OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(joinpath(configs["tmp_folder"], "QFLX2ATM.bin"), OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)
        writeBinary!(joinpath(configs["tmp_folder"], "MASK.bin"), map.mask, buffer2d; endianess=:little_endian)

        sendText(PTI, "OK")

        stage = :RUN
        
    elseif stage == :RUN && msg["MSG"] == "RUN"

        t_flags["new_year"]  = (timeinfo[1] != timeinfo_old[1])
        t_flags["new_month"] = (timeinfo[2] != timeinfo_old[2])
        t_flags["new_day"]   = (timeinfo[3] != timeinfo_old[3])

        timeinfo_old[:] = timeinfo

        for i = 1:length(x2o_available_varnames)
            varname = x2o_available_varnames[i]

            readBinary!(
                joinpath(configs["tmp_folder"], varname * ".bin"),
                (x2o_wanted_flag[i]) ? OMDATA.x2o[varname] : null2d,
                buffer2d;
                endianess=:little_endian, delete=false
            )

        end
       
        println("Calling ", OMMODULE.name, " to do MAGICAL calculations")

        Δt = parse(Float64, msg["DT"])
        Δt_substeps = Δt / configs["substeps"]


        cost = @elapsed for substep = 1:configs["substeps"]

            print(format("Substep: {:d}/{:d}\r", substep, configs["substeps"]))

            OMMODULE.run(OMDATA;
                t             = timeinfo,
                t_cnt         = t_cnt,
                t_flags       = t_flags,
                Δt            = Δt_substeps,
                substep       = substep,
                write_restart = ( msg["WRITE_RESTART"] == "TRUE" ) ? true : false,
            )

        end

        println(format("*** It takes {:.2f} secs. ***", cost))

        writeBinary!(joinpath(configs["tmp_folder"], "SST.bin"), OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(joinpath(configs["tmp_folder"], "QFLX2ATM.bin"), OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)

        sendText(PTI, "OK")

        t_cnt += 1

    elseif stage == :RUN && msg["MSG"] == "END"

        # move short_term_archive_files to long term archive directory
        if configs["enable_long_term_archive"]
            println("===== Long term archiving files BEGIN =====")
            sdir = configs["short_term_archive_dir"]
            ldir = configs["long_term_archive_dir"]
            mkpath(ldir)
            for fname in eachline(configs["short_term_archive_list"])
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
        OMMODULE.crash(OMDATA) 
        sendText(PTI, "CRASH")
        throw(ErrorException("Unknown status: stage " * stage * ", MSG: " * String(msg["MSG"])))
    end

    flush(stdout)
end


