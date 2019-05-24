using Statistics: mean




function nanmean(
    x;
    dims :: Tuple,
)

    y = copy(x)

    nan_idx    = isnan.(y)
    finite_idx = isfinite.(y)

    y[nan_idx] .= 0.0

    return sum(y, dims=dims) ./ sum(finite_idx, dims=dims)
end
