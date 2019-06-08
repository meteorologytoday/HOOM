
mutable struct Workspace
    τx   :: Array{Float64}
    τy   :: Array{Float64}
    nswflx :: Array{Float64}
    swflx  :: Array{Float64}
    frwflx :: Array{Float64}
    ifrac :: Array{Float64}
    SST :: Array{Float64}
    sumflx  :: Array{Float64}

    function Workspace(Nx::Integer, Ny::Integer)

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
        )
    end

end


