struct DataBinding
    
    here  :: DataUnit
    there :: DataUnit
   
    # views 
    here_pull_view   :: Any
    there_pull_view  :: Any

    here_push_view   :: Any
    there_push_view  :: Any

    # pull from there
    # push from here

    function DataBinding(
        here        :: DataUnit,
        there       :: DataUnit,

        here_pull_yrng        :: Any,
        there_pull_yrng       :: Any,

        here_push_yrng        :: Any,
        there_push_yrng       :: Any,
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


        here_pull_view  = _helper_genView(here, here_pull_yrng)
        there_pull_view = _helper_genView(there, there_pull_yrng)

        here_push_view  = _helper_genView(here, here_push_yrng)
        there_push_view = _helper_genView(there, there_push_yrng)


        if length(here_pull_view) != length(there_pull_view)
            throw(ErrorException("Range mismatch in pull view."))
        end
 
        if length(here_push_view) != length(there_push_view)
            throw(ErrorException("Range mismatch in push view."))
        end
        
        
        return new(
            here,
            there,
            here_pull_view,
            there_pull_view,
            here_push_view,
            there_push_view,
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
    here_pull_yrng        :: Any,
    there_pull_yrng       :: Any,
    here_push_yrng        :: Any,
    there_push_yrng       :: Any;
    labels           :: Union{Symbol, AbstractArray{Symbol}} = [:ALL],
)

    db = DataBinding(
        here,
        there,
        here_pull_yrng,
        there_pull_yrng,
        here_push_yrng,
        there_push_yrng,
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
    group_label    :: Symbol,
    direction      :: Symbol,
)

    dbs = data_exchanger.data_bindings[group_label]
   

     
    if direction == :PUSH
        for db in dbs
    #        println("db: id: ", db.here.id, "; ", db.there.id)
            db.there_push_view[:] = db.here_push_view
        end 
    elseif direction == :PULL
        for db in dbs
            db.here_pull_view[:] = db.there_pull_view
        end 
    else
        throw(ErrorException("Unknown direction keyword: " * string(direction)))
    end
end

# input dataunit and yrange, output a view of corresponding yrng with shape :xy or :xyz
function _helper_genView(
    du,
    yrng,
)
     
    if du.shape == :xy

        rng = [Colon(), yrng]

        if du.has_Xdim
            push!(rng, Colon())
        end

        new_view = view( du.data, rng...)

    else 

        permute_xyz = [1, 2, 3]
        permute_zxy = [2, 3, 1]

        if du.has_Xdim
            push!(permute_xyz, 4)
            push!(permute_zxy, 4)
        end
        
        if du.shape == :zxy
            new_view = PermutedDimsArray(du.data, permute_zxy)
        else
            new_view = du.data
        end

        # at this point views are arranged in :xyz
        rng  = [Colon(), yrng, Colon()]

        if du.has_Xdim
            push!(rng, Colon())
        end

        new_view = view( new_view, rng...)
    end
   
end


