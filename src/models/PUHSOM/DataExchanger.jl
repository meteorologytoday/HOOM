struct DataBinding
    
    name    :: String

    data_here  :: AbstractArray
    data_there :: AbstractArray
    
    data_here_shape  :: Symbol   
    data_there_shape :: Symbol

    function DataBinding(
        name    :: String,
        data_here  :: AbstractArray{T_here},
        data_there :: AbstractArray{T_there},
        data_here_shape  :: Symbol,
        data_there_shape :: Symbol,
    ) where T_here where T_there

        if length(data_here) != length(data_there)
            println("size(data_here) = ", size(data_here), "; size(data_there) = ", size(data_there))
            throw(ErrorException("Data do not share same length."))
        end


        if data_here_shape == :xyz && data_there_shape == :zxy
            data_there = PermutedDimsArray(data_there, (2, 3, 1))
        elseif data_here_shape == :zxy && data_there_shape == :xyz
            data_there = PermutedDimsArray(data_there, (3, 1, 2))
        elseif data_here_shape != data_there_shape
            throw(ErrorException("Data do not share same shape"))
        elseif ! (data_here_shape in (:xyz, :zxy, :xy))
            throw(ErrorException("Unknown shape: " * string(data_here_shape)))
        end

        if T_here != T_there
            println(format("Warning : data binding {} do not share the same data type. Implicit data conversion will be expected.", name))
        end

        return new(
            name,
            data_here,
            data_there,
            data_here_shape,
            data_there_shape,
        )

    end
end


struct DataExchanger
    
    data_bindings :: AbstractArray{DataBinding}

    function DataExchanger()
        
        data_bindings = Array{DataBinding}(undef, 0)

        return new(data_bindings)
    end

end


function createBinding!(
    data_exchanger   :: DataExchanger,
    name             :: String,
    data_here        :: AbstractArray,
    data_here_shape  :: Symbol,
    data_there       :: AbstractArray,
    data_there_shape :: Symbol,
)
    push!(
        data_exchanger.data_bindings, 
        DataBinding(
            name,
            data_here,
            data_there,
            data_here_shape,
            data_there_shape,
        )
    )

end


function syncData!(
    data_exchanger :: DataExchanger,
    direction      :: Symbol,
)
    if direction == :push
        for data_binding in data_exchanger.data_bindings
            data_binding.data_there[:] = data_binding.data_here
        end 
    elseif direction == :pull
        for data_binding in data_exchanger.data_bindings
            data_binding.data_here[:] = data_binding.data_there
        end 
    else
        throw(ErrorException("Unknown direction keyword: " * string(direction)))
    end
end
