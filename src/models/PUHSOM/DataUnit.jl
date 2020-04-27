mutable struct DataUnit
    id          :: Symbol
    grid        :: Symbol    # could be (f, c, s) x (T, U, V, W)
    shape       :: Symbol    # could be :xy, :xyz, :zxy
    data        :: AbstractArray
    has_Xdim    :: Bool
end
