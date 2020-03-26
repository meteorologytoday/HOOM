include("../../julia_lib/BinaryIO.jl")
include("../../julia_lib/ProgramTunnel.jl")


using Formatting
using Printf
using .ProgramTunnel

TS = defaultTunnelSet(path=".")
reverseRole!(TS)
mkTunnel(TS)

# Mimic
data = rand(Float64, parse(Int, recvText(TS)))
buffer = zeros(UInt8, length(data) * 8)

while true

    println("Try to recv new msg")
    msg = recvText(TS)
    println("Msg recv: [", msg, "]")

    if msg == "<<END>>"
        println("Simulation ends!")
        break
    end

    recvBinary!(TS, data, buffer)
    println("Received binary data: ", data)

    println("Now I am doing some magical SMARTSLAB computation.")
    sst_fn = format("SST_{:03d}.nc", convert(Integer, floor(rand() * 1000)))

    data .+= 100.0
    #sleep(1)

    println("Gonna send SST file back : ", sst_fn)
    sendText(TS, sst_fn)

    
    sendBinary!(TS, data, buffer)
    println("Received binary data: ", data)


end



