include("../../julia_lib/BinaryIO.jl")
include("../../julia_lib/MailboxPipe2.jl")


using Formatting
using Printf
using .MailboxPipe2

MI = MailboxPipe2.MailboxInfo(".")
MailboxPipe2.mkPipe(MI)

#hello(MI)

# Mimic
data = rand(Float64, parse(Int, MailboxPipe2.recvText(MI)))
buffer = zeros(UInt8, length(data) * 8)


while true

    println("Try to recv new msg")
    msg = MailboxPipe2.recvText(MI)
    println("Msg recv: [", msg, "]")

    if msg == "<<END>>"
        println("Simulation ends!")
        break
    end

    MailboxPipe2.recvBinary!(MI, "FLX", data, buffer)
    println("Received binary data: ", data)

    println("Now I am doing some magical SMARTSLAB computation.")
    sst_fn = format("SST_{:03d}.nc", convert(Integer, floor(rand() * 1000)))

    data .+= 100.0
    sleep(2)

    println("Gonna send SST file back : ", sst_fn)
    MailboxPipe2.sendText(MI, sst_fn)

    
    MailboxPipe2.sendBinary!(MI, "SST", data, buffer)
    println("Received binary data: ", data)


end



