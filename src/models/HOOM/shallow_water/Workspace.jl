mutable struct Workspace

    fT :: Array
    fU :: Array
    fV :: Array

    cT :: Array
    cU :: Array
    cV :: Array

    sT :: Array
    sU :: Array
    sV :: Array


    ptr :: Dict


    function Workspace(;
        Nx :: Int64,
        Ny :: Int64,
        Nz_f :: Int64,
        Nz_c :: Int64,
        fT :: Int64=0,
        fU :: Int64=0,
        fV :: Int64=0,
        cT :: Int64=0,
        cU :: Int64=0,
        cV :: Int64=0,
        sT :: Int64=0,
        sU :: Int64=0,
        sV :: Int64=0,
    )


        _fT = []
        _fU = []
        _fV = []

        _cT = []
        _cU = []
        _cV = []

        _sT = []
        _sU = []
        _sV = []

        for i=1:fT
            push!(_fT, zeros(Float64, Nx, Ny, Nz_f))
        end 

        for i=1:fU
            push!(_fU, zeros(Float64, Nx, Ny, Nz_f))
        end 

        for i=1:fV
            push!(_fV, zeros(Float64, Nx, Ny+1, Nz_f))
        end 

        for i=1:cT
            push!(_cT, zeros(Float64, Nx, Ny, Nz_c))
        end 

        for i=1:cU
            push!(_cU, zeros(Float64, Nx, Ny, Nz_c))
        end 

        for i=1:cV
            push!(_cV, zeros(Float64, Nx, Ny+1, Nz_c))
        end 

        for i=1:sT
            push!(_sT, zeros(Float64, Nx, Ny))
        end 

        for i=1:sU
            push!(_sU, zeros(Float64, Nx, Ny))
        end 

        for i=1:sV
            push!(_sV, zeros(Float64, Nx, Ny+1))
        end 


        ptr = Dict(
            :fT => 1,
            :fU => 1,
            :fV => 1,
            :cT => 1,
            :cU => 1,
            :cV => 1,
            :sT => 1,
            :sU => 1,
            :sV => 1,
        )


        return new(
            _fT, _fU, _fV,
            _cT, _cU, _cV,
            _sT, _sU, _sV,
            ptr,
        )

    end
    

end

function getSpace!(
    wksp :: Workspace,
    grid :: Symbol;
    flat :: Bool = false
)
    i = wksp.ptr[grid]
    list = getfield(wksp, grid)
    wksp.ptr[grid] += 1

    if i > length(list)
        throw(ErrorException("Running out of workspace of " * string(grid)))
    end

    return ( flat ) ? view(list[i], :) : list[i]
end

function reset!(
    wksp :: Workspace,
    grid :: Symbol=:ALL,
)
    if grid == :ALL
        for k in keys(wksp.ptr)
            wksp.ptr[k] = 1
        end
    else
        wksp.ptr[grid] = 1
    end

end

