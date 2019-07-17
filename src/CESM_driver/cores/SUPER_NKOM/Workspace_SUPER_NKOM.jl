
mutable struct Workspace
    taux   :: Array{Float64}
    tauy   :: Array{Float64}
    fric_u :: Array{Float64}
    weighted_fric_u :: Array{Float64}
    nswflx :: Array{Float64}
    swflx  :: Array{Float64}
    sumflx :: Array{Float64}
    ifrac  :: Array{Float64}
    frwflx :: Array{Float64}
    qflx   :: Array{Float64}

    function Workspace(Nx::Integer, Ny::Integer, Nz::Integer)

        ref = zeros(Float64, Nx, Ny)

        return new(
            ref,
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref),
            copy(ref)
        )
    end

end


