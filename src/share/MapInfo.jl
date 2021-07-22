module ModelMap

    export MapInfo

    using NCDatasets

    float = Union{Float64, Float32}

    mutable struct MapInfo{T <: float}
        nx    :: Integer
        ny    :: Integer
        lsize :: Integer
     
        xc    :: Array{T, 2}
        yc    :: Array{T, 2}
     
        xv    :: Array{T, 3}
        yv    :: Array{T, 3}

       
        mask  :: Array{T, 2}
        area  :: Array{T, 2}
        frac  :: Array{T, 2}
        
        missing_value :: T

        function MapInfo{T}(
            filename::String;
            missing_value::T = 1e20
        ) where T <: float
        
            ds = Dataset(filename, "r")
            _mask = ds["mask"][:, :]
            _area = ds["area"][:, :]
            _frac = ds["frac"][:, :]
            _xc  = ds["xc"][:, :]
            _yc  = ds["yc"][:, :]
            _xv  = ds["xv"][:]
            _yv  = ds["yv"][:]

            _nx  = ds.dim["ni"]
            _ny  = ds.dim["nj"]
            close(ds)

            return new{T}(
                _nx, _ny, _nx * _ny,
                _xc, _yc, _xv, _yv,
                _mask, _area, _frac,
                missing_value,
            )

        end

    end
end
