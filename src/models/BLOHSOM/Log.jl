
module Log

    using Formatting

    mutable struct LogInfo

        file            :: String
        labels          :: AbstractArray{Symbol}
        msg             :: Dict
        print_to_screen :: Bool

        function LogInfo(
            file,
            stages = [:main],
            ;
            print_to_screen=true
        )
            msg = Dict()
            
            if length(labels) == 0
                throw(ErrorException("length of labels should be at least one."))
            end

            for label in labels
                msg[label] = []
            end

            return new(
                file,
                labels,
                msg,
                print_to_screen,
            )
        end
    end

    function log!(
        li    :: LogInfo,
        label :: Symbol,
        fmt   :: String,
        args...,
     )
        msg1 = Format(fmt, args...)
        msg2 = Format("[{:s}] : {}", string(label), msg1)
        
        if log.print_to_screen
            println(msg2)
        end
        
    end

end
