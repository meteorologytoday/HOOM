
include("TBIO.jl")

module ProgramTunnel_fs

    using Formatting
    using ..TBIO

    export ProgramTunnelInfo, sendData, recvData!

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

        path         :: AbstractString
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
            timeout       :: AbstractFloat                  = 60.0 * 30,
            buffer        :: AbstractFloat                  = 0.1,
            recv_first_sleep_max :: AbstractFloat = 5.00,
            recv_first_sleep :: AbstractFloat = 0.0,
            reverse_role  :: Bool = false,
            rotate        :: Integer = 100,
            error_sleep   :: Float64 = 0.05,
            history_len   :: Integer = 5,
        )

            if chk_freq <= 0.0
                ErrorException("chk_freq must be positive.") |> throw
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

            if reverse_role
                reverseRole!(PTI)
            end
            
            mkpath(path)
            updateFiles!(PTI)        


            return PTI
        end
    end

    function updateFiles!(PTI::ProgramTunnelInfo)
        recv_fns = []
        send_fns = []
        for i = 1:PTI.rotate
            push!(recv_fns, joinpath(PTI.path, format("{:s}_{:03d}.tb", PTI.recv_fn, i)))
            push!(send_fns, joinpath(PTI.path, format("{:s}_{:03d}.tb", PTI.send_fn, i)))
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
            PTI.recv_trackno += 1

        elseif which == :send
            rm( PTI.send_fns[mod( PTI.send_trackno - 1 - PTI.history_len, PTI.rotate) + 1], force=true)
            PTI.send_trackno += 1

        end
    end

    function sendData(
        PTI  :: ProgramTunnelInfo,
        msg  :: AbstractString,
        arrs :: AbstractArray,
    )

        send_fn = PTI.send_fns[mod(PTI.send_trackno - 1, PTI.rotate) + 1]

        while true
            try
            
                writeTB(
                    send_fn,
                    format("{:d}#{:s}", PTI.send_trackno, msg),
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
        recv_fn = PTI.recv_fns[mod(PTI.recv_trackno - 1, PTI.rotate) + 1]

        println(format("[recvData!] Expecting file: {:s}", recv_fn))
        get_through = false
        sleep(PTI.recv_first_sleep)

        msg_get = false

        
        if isfile(recv_fn)
            PTI.recv_first_sleep -= PTI.chk_freq
            PTI.recv_first_sleep = max(0.0, PTI.recv_first_sleep)
            get_through = true
            println("[recvData] Message is already there. Adjust recv_first_sleep to : ", PTI.recv_first_sleep)
        else
            for cnt in 1:(PTI.timeout_limit_cnt - PTI.recv_first_cnt)

                sleep(PTI.chk_freq)

                if isfile(recv_fn)
                    get_through = true

                    if cnt <= PTI.buffer_cnt

                        println("[recvData] Good guess of the recv_first_sleep : ", PTI.recv_first_sleep)

                    elseif PTI.recv_first_sleep < PTI.recv_first_sleep_max

                        # Out of buffer, need to adjust: increase PTI.recv_first_sleep
                        PTI.recv_first_sleep += PTI.chk_freq 
                        PTI.recv_first_sleep = min(PTI.recv_first_sleep_max, PTI.chk_freq)
                        println("[recvData] Out of buffer. Adjust recv_first_sleep to : ", PTI.recv_first_sleep)

                    else
                        println("[recvData] Out of buffer. But reach to recv_first_sleep_max : ", PTI.recv_first_sleep)
                    end
                        
                    break

                end

            end
        end

        if ! get_through
            ErrorException("[recvData] No further incoming message within timeout.") |> throw
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

                recv_no, msg = split(msg, "#")
                if parse(Int, recv_no) != PTI.recv_trackno
                    println(format("Recive file trackno does not match. Expect {:d} but got {:s}", PTI.recv_trackno, recv_no))
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


            break
        end

        incTrackno(PTI, :recv)
        
        return msg 
    end



end
