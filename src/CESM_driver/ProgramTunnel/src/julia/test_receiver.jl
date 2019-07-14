include("ProgramTunnel_fs.jl")

using .ProgramTunnel_fs
using Formatting

PTI = ProgramTunnelInfo(tag_and_init=[ (:default, 1.8) ])

println(PTI)

for i=1:100
    println(format("[{:3d}]", i))
    msg = recvText(PTI)
end

