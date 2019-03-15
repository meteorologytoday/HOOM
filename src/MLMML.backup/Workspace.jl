
mutable struct Workspace
    taux   :: Array{Float64}
    tauy   :: Array{Float64}
    fric_u :: Array{Float64}
    hflx   :: Array{Float64}
    swflx  :: Array{Float64}
    ifrac  :: Array{Float64}

    function Workspace(Nx::Integer, Ny::Integer)

        ref = zeros(Float64, Nx, Ny)

        return new(ref, copy(ref), copy(ref), copy(ref), copy(ref), copy(ref))
    end

end


