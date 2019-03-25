using Distributed

using Formatting
using Printf
using ArgParse

include("julia_lib/MailboxPipe.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")
include("julia_lib/parseMsg.jl")

using JSON
using .MailboxPipe
using .BinaryIO
using .NetCDFIO


include("config.jl")
include("driver.jl")
