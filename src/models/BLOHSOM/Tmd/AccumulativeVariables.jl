mutable struct AccumulativeVariables

    XFLUX_CONV  :: AbstractArray

    XFLUX_DEN_x :: AbstractArray
    XFLUX_DEN_y :: AbstractArray
    XFLUX_DEN_z :: AbstractArray

    XFLUX_bot   :: AbstractArray
    XFLUX_top   :: AbstractArray


    function AccumulativeVariables(Nx, Ny, Nz, NX)
        return new(
            zeros(Nz,   Nx, Ny,   NX),

            zeros(Nz,   Nx, Ny,   NX),
            zeros(Nz,   Nx, Ny+1, NX),
            zeros(Nz+1, Nx, Ny,   NX),

            zeros(Nx, Ny, NX),
            zeros(Nx, Ny, NX),
        )
    end

end


