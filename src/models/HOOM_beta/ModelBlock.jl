mutable struct ModelBlock

    ev :: Env
    fi :: Field
    co :: Union{Core, Nothing}

    function ModelBlock(
        ev :: Env,
    ) 

        return new(
            ev,
            Field(ev, (ev.id == 0) ? :shared : :local),
            (ev.id == 0) ? nothing : Core(ev),
        )
    end
end


