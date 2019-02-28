
module MailboxPipe
using Formatting

export MailboxInfo, hello, recv, send

default_max_try = 100

mutable struct MailboxInfo

    recv_fn :: AbstractString
    send_fn :: AbstractString

    function MailboxInfo(;
        recv :: AbstractString,
        send :: AbstractString,
    )
        return new(recv, send)
    end

    function MailboxInfo()
        return new("cesm2mymodel.pipe", "mymodel2cesm.pipe")
    end
end

function verify(MI::MailboxInfo)
    if !isfifo(MI.recv_fn) || !isfifo(MI.send_fn)
        throw(ErrorException("Either one or both of the pipe file is not a fifo file."))
    end
end




function appendPath(MI::MailboxInfo, path::AbstractString)
    MI.recv_fn = joinpath(path, MI.recv_fn)
    MI.send_fn = joinpath(path, MI.send_fn)
end

function recv(MI::MailboxInfo)
    local result

    open(MI.recv_fn, "r") do io
        result = strip(read(io, String))
    end

    return result
end

function send(MI::MailboxInfo, msg::AbstractString)

    open(MI.send_fn, "w") do io
        write(io, msg)
    end
end



function hello(MI::MailboxInfo)
    send(MI, "<<TEST>>")
    recv_msg = recv(MI) 
    if recv_msg != "<<TEST>>"
        throw(ErrorException("Weird message: " * recv_msg))
    end
end


end
