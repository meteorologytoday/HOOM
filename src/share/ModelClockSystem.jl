module ModelClockSystem

    export ModelAlarm, ModelClock, advanceClock!, setClock!, clock2str, dt2str, addAlarm!

    using CFTime, Dates, Formatting

    mutable struct ModelAlarm
        name      :: String
        time      :: AbstractCFDateTime
        callbacks :: Array{Function, 1}
        done      :: Bool
    end

    mutable struct ModelClock
        
        name         :: String
        time         :: AbstractCFDateTime
        alarms       :: Array{ModelAlarm, 1}
        alarms_dict  :: Dict
        alarm_ptr    :: Integer 
        
        function ModelClock(
            name :: String,
            time :: AbstractCFDateTime,
        )
            alarms = Array{ModelAlarm}(undef, 0)
            alarms_dict = Dict()
            return new(
                name,
                time + Dates.Second(0),
                alarms,
                alarms_dict,
                0,
            )
        end

    end

    function checkType(time1 :: AbstractCFDateTime, time2 :: AbstractCFDateTime)
        if typeof(time1) != typeof(time2)
            throw(ErrorException("Time type does not match!"))
        end
    end

    function setClock!(
        mc :: ModelClock,
        time :: AbstractCFDateTime,
    )

        checkType(mc.time, time)
        mc.time = time + Second(0)

    end

    function advanceClock!(
        mc :: ModelClock,
        t :: Union{Second, AbstractCFDateTime},
    )
        if typeof(t) <: Second  # inteprete as time interval
            mc.time += t
        elseif typeof(Δt) <: AbstractCFDateTime
            setClock!(mc, t)
        end

        # Only check alarms when there is any
        if length(mc.alarms) > 0
            if mc.alarm_ptr == 0 && mc.time >= mc.alarms[1].time
                mc.alarm_ptr = 1
            end

            while (0 < mc.alarm_ptr <= length(mc.alarms)) && (mc.time >= mc.alarms[mc.alarm_ptr].time)
                    ringAlarm!(mc, mc.alarms[mc.alarm_ptr])
                    mc.alarm_ptr += 1
            end
        end

    end

    function ringAlarm!(mc :: ModelClock, alarm :: ModelAlarm)
        println(format("Alarm '{:s}' rang!", alarm.name))
        alarm.done = true
        for callback in alarm.callbacks
            callback(mc, alarm)
        end
    end

    function clock2str(mc :: ModelClock)
        return dt2str(mc.time)
    end

    function dt2str(dt)
        return format("{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}", Dates.year(dt), Dates.month(dt), Dates.day(dt), Dates.hour(dt), Dates.minute(dt), Dates.second(dt))
    end

    function addAlarm!(
        mc   :: ModelClock,
        name :: String,
        time :: AbstractCFDateTime,
        callback :: Union{Array{Function, 1}, Function, Nothing} = nothing, 
    )

        checkType(time, mc.time)
        
        if callback == nothing
            callbacks = Array{Function, 1}(undef, 0)
        elseif typeof(callback) <: Function
            callbacks = [ callback ]
        end        


        alarm = ModelAlarm(name, time, callbacks, false)
        
        if mc.alarm_ptr == 0
            mc.alarm_ptr = 1
        elseif alarm.time <= mc.time #mc.alarms[mc.alarm_ptr].timestamp
            println("alarm.time = ", dt2str(alarm.time))
            println("mc.alarms[mc.alarm_ptr].time = ", dt2str(mc.alarms[mc.alarm_ptr].time))
            throw(ErrorException("Alarm needs to be set in the future."))
        end
        
        push!(mc.alarms, alarm)
        if ! haskey(mc.alarms_dict, name)
            mc.alarms_dict[name] = []
        end

        push!(mc.alarms_dict[name], alarm)
        sort!(mc.alarms, by = (x)->x.time)

    end

end
