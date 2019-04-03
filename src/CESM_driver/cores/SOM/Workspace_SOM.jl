
mutable struct Workspace
    nswflx :: Array{Float64}
    swflx  :: Array{Float64}
    tfdiv  :: Array{Float64}
    eflx   :: Array{Float64}


    function Workspace(Nx::Integer, Ny::Integer)

        ref = zeros(Float64, Nx, Ny)

        return new(
            ref,
            copy(ref),
            copy(ref),
            copy(ref),
        )
    end

end


