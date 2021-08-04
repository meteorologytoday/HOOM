
module ConfigCheck

    using Formatting

    export ConfigEntry, validateConfig

    mutable struct ConfigEntry
        name      :: Symbol
        required  :: Symbol
        valid_vts :: AbstractArray
        default   :: Any
        desc      :: String
        function ConfigEntry(
            name      :: Symbol,
            required  :: Symbol,
            valid_vts :: AbstractArray,
            default   :: Any = nothing;
            desc      :: String = "",
        )

            return new(
                name,
                required,
                valid_vts,
                default,
                desc,
            )
        end

    end

    """

    infos is a list of four elements:
        1. name of the key, 
        2. A boolean value indicating it is required or not.
        3. A list of valid values or data types. For example: (Bool, ) only allows true/false value, (:on, :off, String) means it can be :on, :off, or any subtype of String
        4. Default value if key is not required and does not exist in the dict yet.

    """
    function validateConfig(
        cfg         :: Dict,
        cfg_entries :: AbstractArray{ConfigEntry, 1},
    )

        new_cfg = Dict()
        
        for entry in cfg_entries

            name      = entry.name
            required  = entry.required
            valid_vts = entry.valid_vts
            default   = entry.default

            if ! ( required in [:optional, :required] )
                throw(ErrorException("The `required` only takes :optional or :required"))
            end

            if haskey(cfg, name)
                pass = false
                for valid_vt in valid_vts

                    if typeof(valid_vt) <: Union{DataType, UnionAll}
                        if typeof(cfg[name]) <: valid_vt
                            pass = true
                        end
                    else
                        if cfg[name] == valid_vt
                            pass = true
                        end
                    end

                    if pass 
                        break
                    end

                end
                
                if pass
                    
                    new_cfg[name] = cfg[name]
                    println(format("[Validation] Config `{:s}` is validated and the value is : {:s}", string(name), string(new_cfg[name])))
                else
                    throw(ErrorException(format(
                        "[Error] Invalid value of key `{:s}`: {:s}. Valid values/types: `{:s}`.",
                        string(name),
                        string(cfg[name]),
                        join(string.(valid_vts), "` ,`")
                    )))
                end


            else

                msg = format(
                    "Missing config: `{:s}`. Valid values/types: `{:s}`.",
                    string(name),
                    join(string.(valid_vts), "` ,`")
                )
                
                if required == :required
                    throw(ErrorException(format("[Required] {:s}", msg)))
                else
                    new_cfg[name] = default
                    println(format("[Optional] {:s} Set to default: {:s}", string(name), string(new_cfg[name])))

                end
            end
        end

        #dropped_names = filter(x -> !( haskey(new_cfg, x) ), keys(cfg))
        #for dropped_name in dropped_names
        #    println(format("The config `{:s}` is not used.", string(dropped_name)))
        #end

        return new_cfg

    end

end
