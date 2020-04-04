
module ModelTime

    mutable struct TimeInfo

        t            :: AbstractArray{Integer, 1}
        t_old        :: AbstractArray{Integer, 1}
        t_flags      :: Dict

        function TimeInfo(
            init_t :: Union{Nothing, AbstractArray{Integer,1}} = nothing,
        )

            if init_t == nothing
                t = zeros(Integer, 4)
            else
                if length(init_t) == 4
                    t = copy(init_t)
                else
                    throw(ErrorException("Length of `t` array must be 4."))
                end
            end


            t_old = copy(t)

            t_flags = Dict(
                :new_year  => true,
                :new_month => true,
                :new_day   => true,
            )

            return new(
                t,
                t_old,
                t_flags,
            )
        end
    end

    function reset!(
        ti :: TimeInfo,
        t  :: AbstractArray{Integer, 1},
    )
        ti.t[:]     = t
        ti.t_old[:] = t

        merge!(ti.t_flags, Dict(
            :new_year  => true,
            :new_month => true,
            :new_day   => true,
        ))

    end

    function update!(
        ti :: TimeInfo,
        t  :: AbstractArray{Integer, 1},
    )
        if length(t) != 4
            throw(ErrorException("Length of `t` array must be 4."))
        end

        ti.t_old .= ti.t
        ti.t     .= t
        ti.t_flags[:new_year]  = (ti.t[1] != ti.t_old[1])
        ti.t_flags[:new_month] = (ti.t[2] != ti.t_old[2])
        ti.t_flags[:new_day]   = (ti.t[3] != ti.t_old[3])


    end

end
