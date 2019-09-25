mutable struct AccumulativeVariables

    T_hadvs :: AbstractArray
    T_vadvs :: AbstractArray
    S_hadvs :: AbstractArray
    S_vadvs :: AbstractArray

    ∇∇T      :: AbstractArray
    ∇∇S      :: AbstractArray

    dTdt_ent :: AbstractArray
    dSdt_ent :: AbstractArray


    function AccumulativeVariables(Nx, Ny, Nz)
        return new(
            zeros(Nz, Nx, Ny),
            zeros(Nz, Nx, Ny),
            zeros(Nz, Nx, Ny),
            zeros(Nz, Nx, Ny),
            zeros(Nz, Nx, Ny),
            zeros(Nz, Nx, Ny),
            zeros(Nx, Ny),
            zeros(Nx, Ny),
        )
    end

end


