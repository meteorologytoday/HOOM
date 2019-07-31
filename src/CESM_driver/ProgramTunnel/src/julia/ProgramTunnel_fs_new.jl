
include("TBIO.jl")

module ProgramTunnel_fs

    using Formatting
    using ..TBIO

    export ProgramTunnelInfo, hello, recvText, sendText, reverseRole!

    mutable struct ProgramTunnelInfo

        nchars       :: Integer
        
        recv_fn      :: AbstractString
        send_fn      :: AbstractString

        chk_freq          :: AbstractFloat
        timeout           :: AbstractFloat
        timeout_limit_cnt :: Integer
        buffer_cnt        :: Integer

        recv_first_sleep_max :: AbstractFloat
        recv_first_sleep :: AbstractFloat
        recv_first_cnt   :: Integer

        rotate       :: Integer
        recv_trackno :: Integer
        send_trackno :: Integer

        path         :: AbstractArray
        error_sleep  :: Float64

        history_len  :: Integer

        recv_fns     :: AbstractArray
        send_fns     :: AbstractArray


        function ProgramTunnelInfo(;
            nchars        :: Integer            = 512,
            recv_fn       :: AbstractString     = "Y2X",
            send_fn       :: AbstractString     = "X2Y",
            chk_freq      :: AbstractFloat                  = 0.05,
            path          :: Union{AbstractString, Nothing} = "x_tmp",
            timeout       :: AbstractFloat                  = 10.0,
            buffer        :: AbstractFloat                  = 0.1,
            recv_first_sleep_max :: AbstractFloat = 5.00,
            recv_first_sleep :: AbstractFloat = 0.0,
            reverseRole   :: Bool = false,
            rotate        :: Integer = 100,
            error_sleep   :: Float64 = 1.0,
            history_len   :: Integer = 20,
        )

            if chk_freq <= 0.0
                ErrorException("chk_freq must be positive.") |> throw
            end
            

            for i = 1:rotate
                push!(paths, joinpath(MPI.path, format("{:03d}")))
            end

            PTI = new(
                nchars,
                recv_fn,
                send_fn,
                chk_freq,
                timeout,
                ceil(timeout / chk_freq),
                ceil(buffer / chk_freq),
                recv_first_sleep_max,
                recv_first_sleep,
                ceil(recv_first_sleep / chk_freq),
                rotate, 1, 1,
                path,
                1.0, history_len,
                [], [],
            )

            if path != nothing
                appendPath(PTI, path)
            end

            if reverseRole
                reverseRole!(PTI)
            end
            
            makepath(path)
            updateFiles!(PTI)        


            return PTI
        end
    end

    function updateFiles!(PTI::ProgramTunnelInfo)
        recv_fns = []
        send_fns = []
        for i = 1:PTI.rotate
            push!(recv_fns, joinpath(PTI.path, format("{:s}_{:03d}.tb", PTI.recv_fn)))
            push!(send_fns, joinpath(PTI.path, format("{:s}_{:03d}.tb", PTI.send_fn)))
        end

        PTI.recv_fns = recv_fns
        PTI.send_fns = send_fns
    end

    function cleanHistory(PTI::ProgramTunnelInfo)
        cleanDir(PTI, ((PTI.send_trackno-1 - PTI.history_length -1) % PTI.rotate) + 1)
    end


    function reverseRole!(PTI::ProgramTunnelInfo)
        PTI.recv_fn, PTI.send_fn = PTI.send_fn, PTI.recv_fn
    end

    function incTrackno(PTI::ProgramTunnelInfo, which::Symbol)
        if which == :recv
            rm( PTI.recv_fns[mod( PTI.recv_trackno - 1 - PTI.history_len, PTI.rotate) + 1], force=true)
            PTI.recv_trackno = mod(PTI.recv_trackno, PTI.rotate) + 1

        elseif which == :send
            rm( PTI.send_fns[mod( PTI.send_trackno - 1 - PTI.history_len, PTI.rotate) + 1], force=true)
            PTI.send_trackno = mod(PTI.send_trackno, PTI.rotate) + 1
        end
    end

    function sendData(
        PTI  :: ProgramTunnelInfo,
        msg  :: AbstractString,
        arrs :: AbstractArray,
    )

        while true
            try
                writeTB(
                    PTI.send_fns[PTI.send_trackno],
                    format("{:d}#", PTI.send_trackno),
                    PTI.nchars,
                    arrs,
                )
            catch ex
                println(string(ex))
                println("keep sending msg.")
                sleep(PTI.error_sleep)
            end
            break
        end

        incTrackno(PTI, :send)

    end


    function recvData!(
        PTI  :: ProgramTunnelInfo,
        arrs :: AbstractArray,
    )
        recv_fn = PTI.recv_fns[PTI.recv_trackno]

        get_through = false
        sleep(PTI.recv_first_sleep)

        msg_get = false

        
        if isfile(recv_fn)
            PTI.recv_first_sleep -= PTI.chk_freq
            PTI.recv_first_sleep = max(0.0, PTI.recv_first_sleep)
            get_through = true
            println("[recvText] Message is already there. Adjust recv_first_sleep to : ", PTI.recv_first_sleep)
        else
            for cnt in 1:(PTI.timeout_limit_cnt - PTI.recv_first_cnt)

                sleep(PTI.chk_freq)

                if isfile(recv_fn)
                    get_through = true

                    if cnt <= PTI.buffer_cnt

                        println("[recvText] Good guess of the recv_first_sleep : ", PTI.recv_first_sleep)

                    elseif PTI.recv_first_sleep < PTI.recv_first_sleep_max

                        # Out of buffer, need to adjust: increase PTI.recv_first_sleep
                        PTI.recv_first_sleep += PTI.chk_freq 
                        PTI.recv_first_sleep = min(PTI.recv_first_sleep_max, PTI.chk_freq)
                        println("[recvText] Out of buffer. Adjust recv_first_sleep to : ", PTI.recv_first_sleep)

                    else
                        println("[recvText] Out of buffer. But reach to recv_first_sleep_max : ", PTI.recv_first_sleep)
                    end
                        
                    break

                end

            end
        end

        if ! get_through
            ErrorException("[recvText] No further incoming message within timeout.") |> throw
        end

        local msg

        while true
            try
                msg = readTB!(
                    recv_fn,
                    PTI.nchars,
                    arrs,
                )

                if msg == nothing

                    println("File does not exist or not the expected size.")
                    println("Keep receiving...")
                    sleep(PTI.error_sleep)
                    continue
                end

            catch ex
                println(string(ex))
                println("Keep receiving...")
                sleep(PTI.error_sleep)
                continue
            end

            recv_no, msg = split(msg, "#")
            if parse(Int, recv_no) != PTI.recv_trackno
                println(format("Recive file trackno does not match. Expect {:d} but got {:s}", PTI.recv_trackno, recv_no))
                println("Keep receiving...")
                sleep(PTI.error_sleep)
                continue
            end

            break
        end

        incTrackno(PTI, :recv)
        
        return msg 
    end



end
