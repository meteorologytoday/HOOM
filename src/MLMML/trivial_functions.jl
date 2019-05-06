
#=
"""
This function checks if CFL criteria is satisfied which is required by Euler Forward Scheme. Explicitly,

 K Δt       1
------  <= ---
(Δz)^2      2

for every layer. This function returns true every layer is satisfied, returns false if any of the layers is not.

"""
function checkAllDiffusionStability(;
    Δzs:: AbstractArray{Float64, 1},
    K  :: Float64,
    Δt :: Float64,
)

    return all( Δzs .>= √(2.0 * K * Δt) )
end

function checkDiffusionStability(;
    Δz :: Float64,
    K  :: Float64,
    Δt :: Float64,
)

    return Δz >= √(2.0 * K * Δt)
end

function checkDiffusionStability(oc::OceanColumn; Δt)
    return checkAllDiffusionStability(Δzs=oc.Δzs, K=oc.K, Δt=Δt)
end


function minΔz(;
    K :: Float64,
    Δt:: Float64,
)
    return √(2.0 * K * Δt)
end

=#

function boundMLD(h_ML::Float64; h_ML_min::Float64, h_ML_max::Float64)
    return max(min(h_ML, h_ML_max), h_ML_min)
end


"""
    getTKE(fric_u)

# Description
This function returns the TKE (turbulent kinetic energy) `k = 0.5 * (v'^2)` of ML. This parameterization is given by Kim 1976: "A Generalized Bulk Model of the Oceanic Mixed Layer" in its equation (11)

"""
function getTKE(;
    fric_u :: Float64
)
    cm = max(3e-2, 3.0 * fric_u)
    return  0.5 * cm^2.0
end


function updateFLDO!(
        occ :: OceanColumnCollection,
    )
    for i=1:occ.Nx, j=1:occ.Ny
        OC_updateFLDO!(occ, i, j)
    end
end


function OC_updateFLDO!(
        occ :: OceanColumnCollection,
        i   :: Integer,
        j   :: Integer,
    )
    occ.FLDO[i, j] = getFLDO(zs=occ.zs, h_ML=occ.h_ML[i, j])
end

"""

    Returns the FLDO. If mixed-layer depth is equal to the total depth
    of ocean column, -1 will be returned.

"""
function getFLDO(;
    zs   :: AbstractArray{Float64,1},
    h_ML :: Float64,
)
    for i = 1:length(zs)-1
        #println("h:", h, "; Δzs= ", zs[1] - zs[i+1])
        if h_ML < (zs[1] - zs[i+1])  # I don't use equality in order to avoid Δb = 0 during some initialization
            return i
        end
    end

    # 
    return -1
    #throw(ErrorException("h_ML cannot be equal or greather than -z[end]"))
end

function getWindStress(;
    u10::Float64
)

    return u10 * 1e-3 * ( (u10 < 25.0) 
                    ? 2.7 + 0.142 * u10 + 0.0764 * u10^2.0
                    : u10 * (2.16 + 0.5406 * (1.0 - exp(- (u10 - 25.0) / 7.5)))
    )

end

function getFricU(;
    ua::Float64
)
    return √(getWindStress(u10=ua) / ρ)
end

function TS2b(T::Float64, S::Float64)
    return g * (α * (T - T_ref) - β * (S - S_ref))
end

function OC_updateB!(
    occ :: OceanColumnCollection,
    i   :: Integer,
    j   :: Integer,
)

    occ.b_ML[i, j] = TS2b(occ.T_ML[i, j], occ.S_ML[i, j])
    for k=1:occ.Nz
        occ.bs[i, j, k] = TS2b(occ.Ts[i, j, k], occ.Ss[i, j, k])
    end

end

function updateB!(occ::OceanColumnCollection)

    for i=1:occ.Nx, j=1:occ.Ny
        OC_updateB!(occ, i, j)
    end

end
