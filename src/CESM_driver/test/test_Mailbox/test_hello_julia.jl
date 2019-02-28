include("../julia_lib/Mailbox.jl")

using Formatting
using Printf
using .Mailbox

recv_fifo = "cesm2mymodel.fifo"
send_fifo = "mymodel2cesm.fifo"


MI = MailboxInfo()

println("Program start.")

send(MI, "xxx", 5)
hello(MI; max_try=5)

println("Program ends.")

