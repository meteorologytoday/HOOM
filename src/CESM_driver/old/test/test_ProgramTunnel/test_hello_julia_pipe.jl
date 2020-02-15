include("../../julia_lib/MailboxPipe.jl")

using Formatting
using Printf
using .MailboxPipe

MI = MailboxInfo()

println("Program start.")

hello(MI)

println("Program ends.")

