mutable struct DataTable

    ev         :: Env
    data_units :: Dict{ Symbol, DataUnit }
    flags      :: Dict{ Symbol, Int64 }    # What are flags for? 

    function DataTable(
        ev::Env
    )

        data_units = Dict{Symbol, DataUnit}()
        flags      = Dict{Symbol, Int64}()

        return new(
            ev,
            data_units,
            flags,
        )

    end

end

function regVariable!(
    dt       :: DataTable,
    id       :: Symbol,
    grid     :: Symbol,
    dtype    :: DataType;
    data     :: Union{Nothing, SharedArray{T}} = nothing,
) where T

    ev = dt.ev
    Nx, Ny, Nz = ev.gd.Nx, ev.gd.Ny, ev.gd.Nz

    if haskey(dt.data_units, id)
        throw(ErrorException("Error: variable id " * String(id) *  " already exists."))
    end

    dim = Dict(
        :T  => [Nz,   Nx  , Ny  ],
        :U  => [Nz,   Nx  , Ny  ],
        :V  => [Nz,   Nx  , Ny+1],
        :W  => [Nz+1, Nx  , Ny  ],
        :sT => [1, Nx, Ny  ],
        :sU => [1, Nx, Ny  ],
        :sV => [1, Nx, Ny+1],
        :sW => [1, Nx, Ny  ],
    )[grid]

    if ! (dtype in (Float64, Int64))
        throw(ErrorException("Invalid data type. Only Float64 and Int64 are accepted"))
    end

    if data == nothing
        println("Create SharedArray of " * string(id))
        data = SharedArray{dtype}(dim...)
    else

        if Tuple(dim) != size(data)
            println("Expect ", dim)
            println("Get ", size(data))
            throw(ErrorException("Provided data does not have correct dimension."))
        end

        if dtype != T
            throw(ErrorException("dtype and provided data does not match."))
        end
    end

    dt.data_units[id] = DataUnit(
        id,
        grid,
        data,
    )

    dt.flags[id] = 0
end
