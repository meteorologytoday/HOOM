include("BinaryIO.jl")
include("MailboxPipe2.jl")



MI = MailboxPipe2.MailboxInfo(".")

n = parse(Int, MailboxPipe2.recvText(MI))

data = zeros(Float64, n)
buffer = zeros(UInt8, length(data) * 8)


MailboxPipe2.recvBinary!(MI, data, buffer)

println("What I got here?")
println(data)


