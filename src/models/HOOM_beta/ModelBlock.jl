mutable struct ModelBlock

    ev :: Env
    fi :: Field
    co :: Union{Core, Nothing}

    function ModelBlock(
        ev :: Env;
        init_core :: Bool = false,
    ) 

        return new(
            ev,
            Field(ev),
            (init_core) ? Core(ev) : nothing,
        )
    end
end


