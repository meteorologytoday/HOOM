mutable struct DataUnit
    id          :: Symbol
    grid        :: Symbol    # could be ('', s) x (T, U, V, W, UV)
    data        :: AbstractArray  # original data
    odata       :: AbstractArray  # oriented data with shape (x, y, z, X)
    
    function DataUnit(
        id          :: Symbol,
        grid        :: Symbol,   # could be (f, c, s) x (T, U, V, W)
        data        :: AbstractArray,
    )
    
        odata = PermutedDimsArray(data, [2,3,1])

        return new(
            id,  
            grid, 
            data,
            odata,
        )
    end        
end
