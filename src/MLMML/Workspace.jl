
mutable struct Workspace
    taux   :: Array{Float64, 1}
    tauy   :: Array{Float64, 1}
    fric_u :: Array{Float64, 1}
    hflx   :: Array{Float64, 1}
    swflx  :: Array{Float64, 1}
    ifrac  :: Array{Float64, 1}

    function Workspace(N::Integer)

        ref = zeros(Float64, N)

        return new(ref, copy(ref), copy(ref), copy(ref), copy(ref), copy(ref))
    end

end


