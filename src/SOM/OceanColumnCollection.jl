mutable struct OceanColumnCollection
    Nx       :: Integer           # Number of columns in i direction
    Ny       :: Integer           # Number of columns in j direction
    
    Kh_T     :: Float64           # Horizontal diffusion coe of temperature

    mask     :: AbstractArray{Float64, 2}
    mask_idx :: Any

    T_ML     :: AbstractArray{Float64, 2}
    h_ML     :: AbstractArray{Float64, 2}
    qflx2atm :: AbstractArray{Float64, 2} # The energy flux to atmosphere if freezes



    # Derived quantities
    N_ocs  :: Integer           # Number of columns

    function OceanColumnCollection(;
        Nx     :: Integer,
        Ny     :: Integer,
        Kh_T   :: Float64,
        T_ML   :: Float64,
        h_ML   :: Float64,
        mask   :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
    )

        if mask == nothing
            mask = SharedArray{Float64}(Nx, Ny)
            mask .+= 1.0
        else
            mask = copy(mask)
        end

        mask_idx = (mask .== 0.0)


        _T_ML     = SharedArray{Float64}(Nx, Ny)
        _h_ML     = SharedArray{Float64}(Nx, Ny)
        qflx2atm  = SharedArray{Float64}(Nx, Ny)

        _T_ML .= T_ML
        _h_ML .= h_ML
       
        occ = new(
            Nx, Ny, 
            Kh_T,
            mask, mask_idx,
            _T_ML, _h_ML, qflx2atm,
            Nx * Ny
        )

        return occ
    end

end

function copyOCC!(fr_occ::OceanColumnCollection, to_occ::OceanColumnCollection)

    if (fr_occ.Nx, fr_occ.Ny) != (to_occ.Nx, to_occ.Ny)
        throw(ErrorException("These two OceanColumnCollection have different dimensions."))
    end
    
    to_occ.Kh_T = fr_occ.Kh_T
    to_occ.T_ML[:, :]      = fr_occ.T_ML
    to_occ.h_ML[:, :]      = fr_occ.h_ML
    to_occ.qflx2atm[:, :]  = fr_occ.qflx2atm
end


function copyOCC(occ::OceanColumnCollection)
    occ2 = makeBlankOceanColumnCollection(occ.Nx, occ.Ny; mask=mask)
    copyOCC!(occ, occ2)
    
    return occ2
end


function makeBlankOceanColumnCollection(
    Nx   :: Integer,
    Ny   :: Integer;
    mask :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
)

    return OceanColumnCollection(;
        Nx   = Nx,
        Ny   = Ny,
        Kh_T = 0.0,
        T_ML = 0.0,
        h_ML = 0.0,
        mask = mask,
    )
end
