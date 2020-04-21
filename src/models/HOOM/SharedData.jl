
mutable struct DataMapping

    # dim can be either :xyz or :zxy
    fr_shape  :: Symbol
    to_shape  :: Symbol
    fr      :: Symbol
    to      :: Symbol
    
end

mutable struct DataUnit

    id          :: Symbol

    # could be (f, c, s) x (T, U, V, W)
    grid        :: Symbol

    # could be :xy, :xyz, :zxy
    shape         :: Symbol
 
    data        :: AbstractArray
    
end


mutable struct SharedData

    env        :: OcnEnv
    data_units :: Dict{ Symbol, DataUnit }

    function SharedData(
        env::OcnEnv
    )


        data_units = Dict{Symbol, DataUnit}()

        return new(
            env,
            data_units,
        )

    end

end

function regVariable!(
    sd       :: SharedData,
    id       :: Symbol,
    grid     :: Symbol,
    shape    :: Symbol,
    dtype    :: DataType,
)

    env = sd.env
    Nx, Ny, Nz_f, Nz_c = env.Nx, env.Ny, env.Nz_f, env.Nz_c

    dim = Dict(
        :fT => [Nx  , Ny  , Nz_f],
        :fU => [Nx+1, Ny  , Nz_f],
        :fV => [Nx  , Ny+1, Nz_f],
        :fW => [Nx  , Ny  , Nz_f+1],
        :cT => [Nx  , Ny  , Nz_c],
        :cU => [Nx+1, Ny  , Nz_c],
        :cV => [Nx  , Ny+1, Nz_c],
        :cW => [Nx  , Ny  , Nz_c+1],
        :sT => [Nx  , Ny],
        :sU => [Nx  , Ny],
        :sV => [Nx  , Ny],
    )[grid]

    if ! (shape in (:xy, :xyz, :zxy))
        throw(ErrorException("Error: only :xy, :xyz, :zxy are accepted in shape"))
    end

    if grid in (:sT, :sU, :sV)
        if shape in (:xyz, :zxy)
            throw(ErrorException("Error: grid and shape mismatch. Slab grid should only be with :xy shape."))
        end
    else
        if shape == :zxy
            circshift!(dim, 1)
        end
    end


    if ! (dtype in (Float64, Int64))
        throw(ErrorException("Invalid data type. Only Float64 and Int64 are accepted"))
    end
    sd.data_units[id] = DataUnit(
        id,
        grid,
        shape,
        SharedArray{dtype}(dim...),
    )
end
