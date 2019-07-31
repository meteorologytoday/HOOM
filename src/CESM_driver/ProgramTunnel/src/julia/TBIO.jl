using Formatting

module TBIO

    export readTB!, writeTB

    function readTB!(
        filename   :: AbstractString,
        txt_nchars :: Integer,
        arrs       :: AbstractArray{T};
        endianess  :: Symbol = :LITTLE,   # Small: 0x04030201,  Big: 0x01020304
    ) where T <: AbstractArray{Float64}

        if endianess != :LITTLE && endianess != :BIG
            throw(ErrorException("Unknown symbol: " * string(endianess)))
        end

        if txt_nchars < 0
            throw(ErrorException("txt_nchars cannot be negative."))
        end

        msg = nothing

        if isfile(filename) &&
           filesize(filename) == (sum(length.(arrs)) * 8 + txt_nchars)

            open(filename, "r") do io
                msg = String(read(io, txt_nchars))
                for i = 1:length(arrs)
                    read!(io, arrs[i])
                end

                if     endianess == :LITTLE && Base.ENDIAN_BOM == 0x01020304
                    for i = 1:length(arrs)
                        arrs[i][:] = ltoh.(arrs[i])
                    end
                elseif endianess == :BIG && Base.ENDIAN_BOM == 0x04030201
                    for i = 1:length(arrs)
                        arrs[i][:] = ntoh.(arrs[i])
                    end
                end
            end
        end

        return msg
    end


    function writeTB(
        filename   :: AbstractString,
        msg        :: AbstractString,
        txt_nchars :: Integer,
        arrs       :: AbstractArray{T};
        endianess  :: Symbol = :LITTLE,
    ) where T <: AbstractArray{Float64}


        if endianess != :LITTLE && endianess != :BIG
            throw(ErrorException("Unknown symbol: " * string(endianess)))
        end

        if length(msg) > txt_nchars
            throw(ErrorException("Message length exceeds txt_nchars."))
        end

        if txt_nchars < 0
            throw(ErrorException("txt_nchars cannot be negative."))
        end

        open(filename, "w") do io
            write(io, msg)
            
            for i = 1:(txt_nchars - length(msg))
                write(io, " ")
            end

            if     endianess == :LITTLE && Base.ENDIAN_BOM == 0x01020304
                for i = 1:length(arrs)
                    write(io, htol.(arrs[i]))
                end
            elseif endianess == :BIG && Base.ENDIAN_BOM == 0x04030201
                for i = 1:length(arrs)
                    write(io, hton.(arrs[i]))
                end
            else
                for i = 1:length(arrs)
                    write(io, arrs[i])
                end
            end
        end
    end

end

using .TBIO

filename = "Test.tb"
msg = "HELLO THIS IS CESM"

obin = [collect(Float64, 1:10) / 2.0]
ibin = [zeros(Float64, 10)]

print(typeof(obin))


writeTB(filename, msg, 100, obin)
println(ibin)
println("==========")
println(readTB!(filename, 100, ibin))
println(ibin)



