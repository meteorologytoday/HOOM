struct DataBinding
    
    here  :: DataUnit
    there :: DataUnit
   
    # views 
    here_view   :: Any
    there_view  :: Any


    function DataBinding(
        here        :: DataUnit,
        there       :: DataUnit,
        here_yrng        :: Any,
        there_yrng       :: Any,
    )

        if ! (here.shape in (:xyz, :zxy, :xy))
            throw(ErrorException("Unknown shape: " * string(here.shape)))
        end

        if ! (there.shape in (:xyz, :zxy, :xy))
            throw(ErrorException("Unknown shape: " * string(there.shape)))
        end


        if here.has_Xdim != there.has_Xdim
            throw(ErrorException("has_Xdim does not match"))
        end

        has_Xdim = here.has_Xdim

        if here.shape == :xy

            here_rng = [Colon(), here_yrng]
            there_rng = [Colon(), there_yrng]

            if has_Xdim
                push!(here_rng, Colon())
                push!(there_rng, Colon())
            end

            here_view  = view( here.data,  here_rng...)
            there_view = view(there.data, there_rng...)

        else 

            permute_xyz = [1, 2, 3]
            permute_zxy = [2, 3, 1]

            if has_Xdim
                push!(permute_xyz, 4)
                push!(permute_zxy, 4)
            end
            

            if here.shape == :zxy
                here_view = PermutedDimsArray(here.data, permute_zxy)
            else
                here_view = here.data
            end

            if there.shape == :zxy
                there_view = PermutedDimsArray(there.data, permute_zxy)
            else
                there_view = there.data
            end

            # at this point views are arranged in :xyz

            here_rng  = [Colon(), here_yrng, Colon()]
            there_rng = [Colon(), there_yrng, Colon()]

            if has_Xdim
                push!(here_rng, Colon())
                push!(there_rng, Colon())
            end

            here_view  = view( here_view,  here_rng...)
            there_view = view(there_view, there_rng...)
        end

        if length(here_view) != length(there_view)
            throw(ErrorException("Range mismatch."))
        end
        
        
        return new(
            here,
            there,
            here_view,
            there_view,
        )

    end

end


struct DataExchanger
    
    data_bindings :: Dict

    function DataExchanger(
        group_labels :: AbstractArray{Symbol} = Array{Symbol}(undef, 0),
    )

        if ! (:ALL in group_labels)
            push!(group_labels, :ALL)
        end


        data_bindings = Dict() 
        for label in group_labels
            data_bindings[label] = Array{DataBinding}(undef, 0)
        end

        return new(
            data_bindings,
        )
    end

end

function createBinding!(
    data_exchanger   :: DataExchanger,
    here             :: DataUnit,
    there            :: DataUnit,
    here_yrng        :: Any,
    there_yrng       :: Any;
    labels           :: Union{Symbol, AbstractArray{Symbol}} = [:ALL],
)

    db = DataBinding(
        here,
        there,
        here_yrng,
        there_yrng,
    )

    if typeof(labels) == Symbol
        labels = [labels]
    end
    
    if ! ( :ALL in labels )
        push!(labels, :ALL)
    end

    for label in labels
        push!(
            data_exchanger.data_bindings[label], 
            db,
        )
    end

end


function syncData!(
    data_exchanger :: DataExchanger,
    group          :: Symbol,
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
