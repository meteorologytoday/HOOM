println("===== Universal Driver Initialization BEGIN =====")

if isdir(wdir)
    cd(wdir)
else
    throw(ErrorException("Working directory [ " * wdir * " ] does not exist."))
end


function parseCESMTIME!(ts::AbstractString, timeinfo::AbstractArray{Integer})

    timeinfo[1] = parse(Int, ts[1:4])
    timeinfo[2] = parse(Int, ts[5:6])
    timeinfo[3] = parse(Int, ts[7:8])
    timeinfo[4] = parse(Int, ts[10:17])
end


output_vars = Dict()

vars_x2o = [
    "SWFLX",
    "HFLX",
    "TAUX",
    "TAUY",
    "IFRAC",
]

vars_o2x = [
    "SST",
    "QFLX",
]

stage = :INIT
mail = MailboxInfo()
#mkPipe(mail)

map = NetCDFIO.MapInfo{Float64}(domain_file)

loop_i = 1
nc_cnt = 1 
output_filename = ""
buffer2d = zeros(UInt8, map.lsize * 8)
timeinfo = zeros(Integer, 4) 

println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
while true

    global OMDATA, stage, loop_i, output_filename, nc_cnt

    end_time = Base.time()

    println(format("Execution time: {:d}", floor(end_time - beg_time)))
    println(format("# Time counter : {:d}", loop_i))
    println(format("# Stage        : {}", String(stage)))

    msg = parseMsg(recv(mail))
    println("==== MESSAGE RECEIVED ====")
    print(json(msg, 4))
    println("==========================")

    if msg["MSG"] in ["INIT", "RUN"]
        parseCESMTIME!(msg["CESMTIME"], timeinfo)
        if loop_i == 1 || (timeinfo[2] == 1 && timeinfo[3] == 1 && timeinfo[4] == 0)
            output_filename = format("SSM_output_{:04d}.nc", Int(timeinfo[1]))
            NetCDFIO.createNCFile(map, output_filename)
            nc_cnt = 1
        end
    end

    if stage == :INIT && msg["MSG"] == "INIT"

        println("===== INITIALIZING MODEL: ", OMMODULE.name , " =====")
        OMDATA = OMMODULE.init(map=map, init_file=init_file, t=timeinfo)
        for (varname, var) in OMDATA.output_vars
            output_vars[varname] = reshape(var, map.nx, map.ny)
        end

        writeBinary!(msg["SST"], OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(msg["QFLX"], OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)
        send(mail, "OK")

        NetCDFIO.write2NCFile(
            map,
            output_filename,
            output_vars;
            time=nc_cnt,
            missing_value=map.missing_value
        )
        nc_cnt += 1


        stage = :RUN
        
    elseif stage == :RUN && msg["MSG"] == "RUN"

        for (varname, container) in OMDATA.x2o
            readBinary!(
                msg[varname],
                container,
                buffer2d;
                endianess=:little_endian,
                delete=false
            )
        end
       
        println("Calling ", OMMODULE.name, " to do MAGICAL calculations")
        OMMODULE.run(OMDATA;
            t     = timeinfo,
            t_cnt = loop_i,
            Δt    = parse(Float64, msg["DT"])
        ) 

        NetCDFIO.write2NCFile(
            map,
            output_filename,
            output_vars;
            time=nc_cnt,
            missing_value=map.missing_value
        )
        nc_cnt += 1

        writeBinary!(msg["SST"], OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(msg["QFLX"], OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)
        send(mail, "OK")

    elseif stage == :RUN && msg["MSG"] == "END"
        OMMODULE.final(OMDATA) 
        
        # Print some report... ?
        
        println("Simulation ends peacefully.")
        break
    else
        OMMODULE.crash(OMDATA) 
        send(mail, "CRASH")
        throw(ErrorException("Unknown status: stage " * stage * ", MSG: " * String(msg["MSG"])))
    end

    loop_i += 1
    flush(stdout)
end


