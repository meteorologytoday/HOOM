using ArgParse
using Printf

include("julia_lib/MailboxPipe.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")
include("julia_lib/parseMsg.jl")

using Formatting
using Printf
using JSON
using .MailboxPipe
using .BinaryIO
using .NetCDFIO


include("config.jl")
include("driver.jl")
