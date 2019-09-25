
using Formatting


"""

infos is a list of four elements:
    1. name of the key, 
    2. A boolean value indicating it is required or not.
    3. A list of valid values or data types. For example: (Bool, ) only allows true/false value, (:on, :off, String) means it can be :on, :off, or any subtype of String
    4. Default value if key is not required and does not exist in the dict yet.

"""
function checkDict!(
    d     :: Any,
    infos :: Any,
)

    for (key, required, valid_vts, default) in infos

        if ! haskey(d, key)
            msg = format(
                "Missing key: `{:s}`. Valid values: `{:s}`.",
                string(key),
                join(valid_vts, "` ,`")
            )

            if required 
                ErrorException(format("[Required] {:s}", msg)) |> throw
            else
                d[key] = default
                println(format("[Optional] {:s} Set to default: {:s}", key, d[key]))

            end

        else
            pass = false
            for valid_vt in valid_vts

                if typeof(valid_vt) <: Union{DataType, UnionAll}

                    if typeof(d[key]) <: valid_vt
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
