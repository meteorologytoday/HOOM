mutable struct OceanColumnCollection
    N_ocs    :: Integer           # Number of columns
    mask     :: Array{Float64}
    mask_idx :: Any

    b_ML     :: Array{Float64}
    b_DO     :: Array{Float64}

    h_ML     :: Array{Float64}
    Q_ML     :: Array{Float64}

    wksp   :: Workspace

    function OceanColumnCollection(;
        N_ocs  :: Integer,
        b_ML   :: Array{Float64},
        b_DO   :: Array{Float64},
        h_ML   :: Array{Float64},
        Q_ML   :: Array{Float64},
        mask   :: Union{Array{Float64}, Nothing} = nothing,
    )
        
        if mask == nothing
            mask = zeros(Float64, N_ocs)
            mask .+= 1.0
        end

        mask_idx = (mask .== 0.0)
        wksp = Workspace(N_ocs)

        return new(N_ocs, mask, mask_idx, b_ML, b_DO, h_ML, Q_ML, wksp)
    end
end



