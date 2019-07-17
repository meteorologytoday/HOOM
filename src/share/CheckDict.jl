
using Formatting

function checkDict!(
    d     :: AbstractArray,
    infos :: AbstractArray,
)

    for (key, required, valid_vts, default) in infos

        if ! haskey(d, key)
            msg = format(
                "Missing key: `{:s}`. Valid values: `{:s}`.",
                key,
                join(valid_vts, "` ,`")
            )

            if required 
                ErrorException(format("[Required] {:s}", msg)) |> throw
            else
                d[key] = default
                println(format("[Optional] {:s} Set to default: {:s}", d[key]))

            end

        else
            pass = false
            for valid_vt in valid_vts

                if typeof(valid_vt) == DataType 
                    if d[key] <: valid_vt
                        pass = true
                    end
                else
                    if d[key] == valid_vt
                        pass = true
                    end
                end

                if pass 
                    break
                end

            end
            
            if ! pass
                ErrorException(format(
                    "[Error] Invalid value of key `{:s}`: {:s}. Valid values/datatype: `{:s}`.",
                    key,
                    d[key],
                    join(valid_vts, "` ,`")
                )) |> throw
            end
        end
    end

end
