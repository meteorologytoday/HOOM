mutable struct OceanColumnCollection
    N_ocs    :: Integer           # Number of columns
    mask     :: Array{Float64}
    mask_idx :: Any

    b_ML     :: Array{Float64}
    b_DO     :: Array{Float64}

    h_ML     :: Array{Float64}
    Q_ML     :: Array{Float64}
    we       :: Array{Float64}

    qflx2atm :: Array{Float64}

    t        :: Array{Float64}
    Δt       :: Array{Float64}
    period   :: Float64
    
    function OceanColumnCollection(;
        N_ocs    :: Integer,
        b_ML     :: Array{Float64},
        b_DO     :: Array{Float64},
        h_ML     :: Array{Float64},
        Q_ML     :: Array{Float64},
        we       :: Array{Float64},
        qflx2atm :: Array{Float64},
        t        :: Array{Float64},
        period   :: Float64,
        mask     :: Union{Array{Float64}, Nothing} = nothing,
    )
        
        if mask == nothing
            mask = zeros(Float64, N_ocs)
            mask .+= 1.0
        end

        mask_idx = (mask .== 0.0)

        Δt = zeros(Float64, length(t))
        Δt[1:end-1] = t[2:end] - t[1:end-1]
        Δt[end] = period - (t[end] - t[1])
        Δt = (circshift(Δt, 1) + Δt)  / 2.0

        return new(N_ocs, mask, mask_idx, b_ML, b_DO, h_ML, Q_ML, we, qflx2atm, t, Δt, period)
    end

end

function makeBlankOceanColumnCollection(;
    N_ocs    :: Integer,
    period_n :: Integer,
    period   :: Float64,
    t        :: Array{Float64},
    mask     :: Union{Array{Float64}, Nothing} = nothing,
)

    if period == nothing
        period = convert(Float64, period_n)
    end
    
    return OceanColumnCollection(
        N_ocs    = N_ocs,
        b_ML     = zeros(Float64, N_ocs),
        b_DO     = zeros(Float64, N_ocs),
        h_ML     = zeros(Float64, N_ocs, period_n),
        Q_ML     = zeros(Float64, N_ocs, period_n),
        we       = zeros(Float64, N_ocs, period_n),
        qflx2atm = zeros(Float64, N_ocs),
        t        = t,
        period   = period,
        mask     = mask
    )

end 



