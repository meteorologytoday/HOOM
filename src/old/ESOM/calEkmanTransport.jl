
@inline function _helperEkmanTransport(
    τ_x :: Float64,
    τ_y :: Float64,
    ϵ   :: Float64,
    f   :: Float64,
)

    return (ϵ * τ_x + f * τ_y) / (ϵ^2 + f^2), (ϵ * τ_y - f * τ_x) / (ϵ^2 + f^2)

end

function calEkmanTransport!(
    occ :: OceanColumnCollection,
    τ_x :: AbstractArray{Float64, 2},
    τ_y :: AbstractArray{Float64, 2},
    M_x :: AbstractArray{Float64, 2},
    M_y :: AbstractArray{Float64, 2},
    f   :: AbstractArray{Float64, 2},
    ϵ   :: AbstractArray{Float64, 2},
)

    for i = 1:occ.Nx, j = 1:occ.Ny    
        M_x[i, j], M_y[i, j] = _helperEkmanTransport(
            τ_x[i, j],
            τ_y[i, j],
            ϵ[i, j],
            f[i, j]
        )
    end
 
end


