println("===== Universal Driver Initialization BEGIN =====")

vars_from_CESM = [
    "SWFLX",
    "HFLX",
    "TAUX",
    "TAUY",
]

vars_to_CESM = [
    "SST",
    "QFLX",
]

if isdir(wdir)
    cd(wdir)
else
    throw(ErrorException("Working directory [ " * wdir * " ] does not exist."))
end

stage = :INIT
mail = MailboxInfo()
map = NetCDFIO.MapInfo{Float64}(domain_file)

time_i = 1 
output_filename = ""
buffer2d  = zeros(UInt8, map.lsize * 8)

println("===== INITIALIZING MODEL: ", OMMODULE.name , " =====")
OMDATA = OMMODULE.init(map)
println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
while true

    global OMDATA, stage, time_i, output_filename

    end_time = Base.time()

    println(format("Execution time: {:d}", floor(end_time - beg_time)))
    println(format("# Time counter : {:d}", time_i))
    println(format("# Stage        : {}", String(stage)))

    msg = parseMsg(recv(mail))
    println("==== MESSAGE RECEIVED ====")
    print(json(msg, 4))
    println("==========================")

    # need to parse time

    if stage == :INIT && msg["MSG"] == "INIT"

        writeBinary!(msg["SST"], OMDATA.sst, buffer2d; endianess=:little_endian)
        send(mail, msg["SST"])
        time_i += 1

        stage = :RUN
        
    elseif stage == :RUN && msg["MSG"] == "RUN"

        for (varname, container) in OMDATA.CESM_containers
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
            t     = [parse(Float64, msg["CESMTIME"])],
            t_cnt = time_i,
            Δt    = parse(Float64, msg["DT"])
        ) 

        writeBinary!(msg["SST_NEW"], OMDATA.sst, buffer2d; endianess=:little_endian)
        send(mail, msg["SST_NEW"])

        time_i += 1

    elseif stage == :RUN && msg["MSG"] == "END"
        OMMODULE.final(OMDATA) 

        println("Simulation ends peacefully.")
        break
    else
        OMMODULE.crash(OMDATA) 
        send(mail, "CRASH")
        throw(ErrorException("Unknown status: stage " * stage * ", MSG: " * String(msg["MSG"])))
    end

    flush(stdout)
end


