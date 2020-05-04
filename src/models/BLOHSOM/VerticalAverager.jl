struct VerticalAverager
    Nz_f    :: Int
    Nz_c    :: Int

    wgt_f2c :: AbstractArray{Float64, 1}
    wgt_c2s :: AbstractArray{Float64, 1}
    
    rng_f2c :: AbstractArray

    function VerticalAverager(;
        z_bnd_f :: AbstractArray{Float64, 1},
        height_level_counts :: AbstractArray{Int64, 1},
    )

        H_f, H_c = _helper_calLayerThickness(z_bnd_f, height_level_counts)

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
#    shape :: Symbol = :xyz,
)

    Nx, Ny, _ = size(var_c)

    for i=1:Nx, j=1:Ny
        for k_c = 1:va.Nz_c
            tmp = 0.0
            
            for k_f in va.rng_f2c[k_c]
                tmp += var_f[i, j, k_f] * va.wgt_f2c[k_f]
            end
            var_c[i, j, k_c] = tmp
        end
    end
end

function calAverage_c2s!(
    va    :: VerticalAverager,
    var_c :: AbstractArray{Float64, 3},
    var_s :: AbstractArray{Float64, 2},
)
    Nx, Ny = size(var_s)

    #println(size(var_s), ";" , size(var_c))

    for i=1:Nx, j=1:Ny
        tmp = 0.0
        for k_c=1:va.Nz_c
            tmp += var_c[i, j, k_c] * va.wgt_c2s[k_c]
        end
        var_s[i, j] = tmp
    end
end

function projVertical_c2f!(
    va    :: VerticalAverager,
    var_c :: AbstractArray{Float64, 3},
    var_f :: AbstractArray{Float64, 3}
)



    Nx, Ny, _ = size(var_c)
    for i=1:Nx, j=1:Ny
        for k_c = 1:va.Nz_c
            var_f[i, j, va.rng_f2c[k_c]] .= var_c[i, j, k_c]
        end
    end
end

function _helper_calLayerThickness(
    z_bnd_f :: AbstractArray{Float64, 1},
    height_level_counts :: AbstractArray{Int64, 1},
)
    H_f  = z_bnd_f[1:end-1] - z_bnd_f[2:end]
    H_c  = zeros(Float64, length(height_level_counts))
   
    if sum(height_level_counts) != length(H_f)
        throw(ErrorException("sum(height_level_counts) != length(z_bnd_f) - 1"))
    end

    idx = 1
    for (k, cnt) in enumerate(height_level_counts)
        H_c[k] = sum(H_f[idx:idx+cnt-1])
        idx += cnt
    end

    return H_f, H_c

    return
end
