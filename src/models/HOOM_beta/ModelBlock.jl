mutable struct ModelBlock

    ev :: Env
    fi :: Field
    co :: Core

    function ModelBlock(ev :: Env)
        return new(
            ev,
            Field(ev),
            Core(ev),
        )
    end
end


