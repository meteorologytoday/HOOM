struct VerticalAverager
    Nz_f    :: Int
    Nz_c    :: Int

    wgt_f2c :: AbstractArray{Float64, 1}
    wgt_c2s :: AbstractArray{Float64, 1}
    
    rng_f2c :: AbstractArray

    function VerticalAverager(;
        H_c :: AbstractArray{Float64, 1},
        H_f :: AbstractArray{Float64, 1},
        height_level_counts :: AbstractArray{Int64, 1},
    )

        Nz_f = length(H_f) 
        Nz_c = length(H_c)
 
        wgt_f2c = zeros(Float64, Nz_f)
        _tmp = 1
        rng_f2c = []
        for (k_c, cnt) in enumerate(height_level_counts)
            rng = _tmp:_tmp+cnt-1
            push!(rng_f2c, rng)
            wgt_f2c[rng] .= H_f[rng] / H_c[k_c]
            _tmp += cnt
        end

        H_total = sum(H_c)
        wgt_c2s = zeros(Float64, Nz_c)
        for k_c in 1:length(H_c)
            wgt_c2s[k_c] = H_c[k_c] / H_total
        end


        return new(
            Nz_f,
            Nz_c,
            wgt_f2c,
            wgt_c2s,
            rng_f2c,
        )
    end
end

function calAverage_f2c!(
    va    :: VerticalAverager,
    var_f :: AbstractArray{Float64, 3},
    var_c :: AbstractArray{Float64, 3};
    mask  :: AbstractArray{Float64, 2}, 
)
    Nx, Ny = size(mask)
    for i=1:Nx, j=1:Ny
        if mask[i, j] == 0
            continue
        end

        for k_c = 1:va.Nz_c
            var_c[k_c, i, j] = mean(var_f[va.rng_f2c[k_c], i, j])
        end
    end
end

function calAverage_c2s!(
    va    :: VerticalAverager,
    var_c :: AbstractArray{Float64, 3},
    var_s :: AbstractArray{Float64, 2};
    mask  :: AbstractArray{Float64, 2}, 
)
    Nx, Ny = size(mask)
    for i=1:Nx, j=1:Ny
        if mask[i, j] == 0
            continue
        end

        var_s[i, j] = mean(var_c[:, i, j])
    end
end

function projVertical_c2f!(
    va    :: VerticalAverager,
    var_c :: AbstractArray{Float64, 3},
    var_f :: AbstractArray{Float64, 3};
    mask  :: AbstractArray{Float64, 2}, 
)
    Nx, Ny = size(mask)
    for i=1:Nx, j=1:Ny
        if mask[i, j] == 0
            continue
        end

        for k_c = 1:va.Nz_c
            var_f[va.rng_f2c[k_c], i, j] .= var_c[k_c, i, j]
        end
    end
end
