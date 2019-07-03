

function LinearRegression(
    x :: AbstractArray{T, 1},
    y :: AbstractArray{T, 1}
) where T <: AbstractFloat 

    N = length(x)

    ϕ = zeros(N, 2)
    ϕ[:, 1] .= 1.0
    ϕ[:, 2] = x

    # ϕ β = y => β = ϕ \ y

    return ϕ \ y

end

function detrend(
    x :: AbstractArray{T, 1},
    y :: AbstractArray{T, 1},
) where T <: AbstractFloat
    β = LinearRegression(x, y)
    return y .- ( β[1] .+ β[2] * x )


end
