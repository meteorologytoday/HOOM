using Statistics: mean
using Printf

function LinearRegression(
    x :: Array{T, 1},
    y :: Array{T, 1}
) where T <: AbstractFloat 

N = length(x)

ϕ = zeros(N, 2)
ϕ[:, 1] .= 1.0
ϕ[:, 2] = x

# ϕ β = y => β = ϕ \ y

return ϕ \ y

end

function rmSeasonalCycle(;
    x,
    period
)

    N = length(x)
    N_PERIODS = Int(N / period)

    t = collect(eltype(x), 1:N)
    x_detrend = detrend(x=x, t=t)

    x_seasonal_sig = mean(reshape(x_detrend, period, :), dims=2)[:, 1]
    x_detrend -= repeat(x_seasonal_sig, outer=(N_PERIODS,))

    return x_detrend
end

function detrend(;
    x,
    t
)
    β = LinearRegression(t, x) 

    return x - (β[1] .+ β[2] * t)
end

function missing2nan!(x)
    x[x .== missing_value] .= NaN
end

function nan2missing!(x)
    x[isnan.(x)] .= missing_value
end

function readModelVar(filename, varname, range_tuple=())
    local var
    @printf("# Reading var [%s]...", varname)
    ds = Dataset(filename,"r")
    v = ds[varname]

    if range_tuple == ()
        range_tuple = repeat([Colon(),], outer=(length(size(v)),))
    end

    v = v[range_tuple...]
    v = nomissing(v, NaN)
    var = convert(Array{dtype}, v)
    close(ds)

    println("done.")

    return var
end

