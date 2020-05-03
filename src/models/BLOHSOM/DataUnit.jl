mutable struct DataUnit
    id          :: Symbol
    grid        :: Symbol    # could be (f, c, s) x (T, U, V, W)
    shape       :: Symbol    # could be :xy, :xyz, :zxy
    data        :: AbstractArray  # original data
    odata       :: AbstractArray  # oriented data with shape (x, y, z, X)
    has_Xdim    :: Bool
    function DataUnit(
        id          :: Symbol,
        grid        :: Symbol,   # could be (f, c, s) x (T, U, V, W)
        shape       :: Symbol,    # could be :xy, :xyz, :zxy
        data        :: AbstractArray,
        has_Xdim    :: Bool,
    )
     
        if ! (shape in (:xyz, :zxy, :xy))
            throw(ErrorException("Unknown shape: " * string(here.shape)))
        end

        if shape == :xy
            odata = data
        else 

            permute_xyz = [1, 2, 3]
            permute_zxy = [2, 3, 1]

            if has_Xdim
                push!(permute_xyz, 4)
                push!(permute_zxy, 4)
            end
            
            if shape == :zxy
                odata = PermutedDimsArray(data, permute_zxy)
            else
                odata = data
            end

        end

        return new(
            id,  
            grid, 
            shape,
            data,
            odata,
            has_Xdim,
        )
    end        
end
