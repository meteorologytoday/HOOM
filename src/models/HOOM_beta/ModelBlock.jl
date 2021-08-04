mutable struct ModelBlock

    ev :: Env
    fi :: Field
    co :: Union{Core, Nothing}

    function ModelBlock(
        ev :: Env;
        init_core :: Bool = false,
        configs   :: Dict,
    ) 

        return new(
            ev,
            Field(ev),
            (init_core) ? Core(ev, configs) : nothing,
        )
    end
end


