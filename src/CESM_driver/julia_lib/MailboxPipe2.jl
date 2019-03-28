
module MailboxPipe2

using Formatting
using ..BinaryIO

export MailboxInfo, hello, recv, send, mkPipe, rmPipe

default_max_try = 100

mutable struct MailboxInfo

    recv_txt_fn :: AbstractString
    send_txt_fn :: AbstractString
    recv_bin_fn :: AbstractString
    send_bin_fn :: AbstractString


    function MailboxInfo(;
        recv_txt :: AbstractString,
        send_txt :: AbstractString,
        recv_bin :: AbstractString,
        send_bin :: AbstractString,
    )
        return new(recv_txt, send_txt, recv_bin, send_bin)
    end

    function MailboxInfo(path::AbstractString)
        MI = new(
            "_cesm2mymodel_txt.fifo",
            "_mymodel2cesm_txt.fifo",
            "_cesm2mymodel_bin.fifo",
            "_mymodel2cesm_bin.fifo"
        )

        appendPath(MI, path)
        return MI
    end
end

function mkPipe(MI::MailboxInfo)
    
    for fn in [MI.send_txt_fn, MI.recv_txt_fn, MI.send_bin_fn, MI.recv_bin_fn]
        if !isfifo(fn)
            println(fn, " is not a fifo or does not exist. Remove it and create a new one.")
            rm(fn, force=true)
            run(`mkfifo $fn`)
        end
    end

end


function reverseRole!(MI::MailboxInfo)
    MI.send_txt_fn, MI.recv_txt_fn, MI.send_bin_fn, MI.recv_bin_fn = MI.recv_txt_fn, MI.send_txt_fn, MI.recv_bin_fn, MI.send_bin_fn
end

function appendPath(MI::MailboxInfo, path::AbstractString)
    MI.recv_txt_fn = joinpath(path, MI.recv_txt_fn)
    MI.send_txt_fn = joinpath(path, MI.send_txt_fn)
    MI.recv_bin_fn = joinpath(path, MI.recv_bin_fn)
    MI.send_bin_fn = joinpath(path, MI.send_bin_fn)

end

function recvText(MI::MailboxInfo)
    local result

    open(MI.recv_txt_fn, "r") do io
        result = strip(read(io, String))
    end

    return result
end

function sendText(MI::MailboxInfo, msg::AbstractString)

    open(MI.send_txt_fn, "w") do io
        write(io, msg)
    end
end

function recvBinary!(
    MI        :: MailboxInfo,
    arr       :: AbstractArray{Float64},
    buffer    :: AbstractArray{UInt8};
    endianess :: Symbol=:little_endian,
)
    BinaryIO.readBinary!(MI.recv_bin_fn, arr, buffer; endianess=endianess)
end

function sendBinary!(
    MI        ::MailboxInfo, 
    arr       :: AbstractArray{Float64},
    buffer    :: AbstractArray{UInt8};
    endianess :: Symbol=:little_endian,
)
    BinaryIO.writeBinary!(MI.send_bin_fn, arr, buffer, endianess=endianess)
end




function hello(MI::MailboxInfo)
    send(MI, "<<TEST>>")
    recv_msg = recv_txt(MI) 
    if recv_msg != "<<TEST>>"
        throw(ErrorException("Weird message: " * recv_msg))
    end
end


end
