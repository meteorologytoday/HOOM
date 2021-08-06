mutable struct ModelBlock

    ev :: Env
    fi :: Field
    co :: Union{Core, Nothing}
    dt :: Union{DataTable, Nothing}


    function ModelBlock(
        ev :: Env;
        init_core :: Bool = false,
    ) 

        fi = Field(ev)

        mb = new(
            ev,
            fi,
            nothing,
            nothing,
        )
        
        dt = DataTable(Nz = ev.Nz, Nx = ev.Nx, Ny = ev.Ny)
        for (k, (varref, grid_type)) in HOOM.getDynamicVariableList(mb; varsets=[:ALL,])
            regVariable!(dt, k, grid_type, varref) 
        end

        co = (init_core) ? Core(ev) : nothing

        mb.dt = dt
        mb.co = co

        return mb
    end
end


