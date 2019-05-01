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

#=
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
=#
stage = :INIT

tmp_folder = joinpath(configs["caserun"], "x_tmp")
mkpath(tmp_folder)

PTI = ProgramTunnelInfo(path=tmp_folder)
reverseRole!(PTI)

map = NetCDFIO.MapInfo{Float64}(configs["domain_file"])

loop_i = 1
nc_cnt = 1 
output_filename = ""
buffer2d = zeros(UInt8, map.lsize * 8)
null2d   = zeros(Float64, map.lsize)
timeinfo = zeros(Integer, 4) 

println("===== ", OMMODULE.name, " IS READY =====")

beg_time = Base.time()
while true

    global OMDATA, stage, loop_i, output_filename, nc_cnt

    end_time = Base.time()

    println(format("Execution time       : {:d} s", floor(end_time - beg_time)))
    println(format("# Time Index counter : {:d}", loop_i))
    println(format("# Stage              : {}", String(stage)))

    msg = parseMsg(recvText(PTI))
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
        OMDATA = OMMODULE.init(
            map          = map,
            init_file    = init_file,
            t            = timeinfo,
            configs      = configs,
            read_restart = (msg["READ_RESTART"] == "TRUE") ? true : false,
        )
        
        rm(configs["short_term_archive_list"], force=true)

        for (varname, var) in OMDATA.output_vars
            output_vars[varname] = reshape(var, map.nx, map.ny)
        end

        global x2o_available_varnames  = split(msg["VAR2D"], ",")
        global x2o_wanted_varnames = keys(OMDATA.x2o)
        global x2o_wanted_flag     = [(x2o_available_varnames[i] in x2o_wanted_varnames) for i = 1:length(x2o_available_varnames)]

        sendText(PTI, "OK")

        writeBinary!(joinpath(tmp_folder, "SST.bin"), OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(joinpath(tmp_folder, "QFLX2ATM.bin"), OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)

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

        for i = 1:length(x2o_available_varnames)
            varname = x2o_available_varnames[i]

            readBinary!(
                joinpath(tmp_folder, varname * ".bin"),
                (x2o_wanted_flag[i]) ? OMDATA.x2o[varname] : null2d,
                buffer2d;
                endianess=:little_endian, delete=false
            )

        end
       
        println("Calling ", OMMODULE.name, " to do MAGICAL calculations")

        cost = @elapsed OMMODULE.run(OMDATA;
            t             = timeinfo,
            t_cnt         = loop_i,
            Δt            = parse(Float64, msg["DT"]),
            write_restart = ( msg["WRITE_RESTART"] == "TRUE" ) ? true : false,
        )

        println(format("*** It takes {:.2f} secs. ***", cost))

        NetCDFIO.write2NCFile(
            map,
            output_filename,
            output_vars;
            time=nc_cnt,
            missing_value=map.missing_value,
        )
        nc_cnt += 1

        println("Gonna send \"OK\"")
        sendText(PTI, "OK")
 
        writeBinary!(joinpath(tmp_folder, "SST.bin"), OMDATA.o2x["SST"], buffer2d; endianess=:little_endian)
        writeBinary!(joinpath(tmp_folder, "QFLX2ATM.bin"), OMDATA.o2x["QFLX2ATM"], buffer2d; endianess=:little_endian)


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

    loop_i += 1
    flush(stdout)
end


