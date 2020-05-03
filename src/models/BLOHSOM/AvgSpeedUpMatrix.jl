@inline function flat_i(
    k :: Int64,
    i :: Int64,
    j :: Int64,
    Nz :: Int64,
    Nx :: Int64,
    Ny :: Int64,
)
    return k + (i-1) * Nz + (j-1) * Nz * Nx
end

@inline function cyc(i::Int64, N::Int64)
    return mod(i-1, N) + 1
end

mutable struct AvgSpeedUpMatrix

    cT_wavg_fT    :: AbstractArray{Float64, 2}  # average from fine to coarse grid
    sT_wavg_cT    :: AbstractArray{Float64, 2}  # average from coase grid to slab

    cU_wavg_fU    :: AbstractArray{Float64, 2}  # average from fine to coarse grid
    sU_wavg_cU    :: AbstractArray{Float64, 2}  # average from coase grid to slab

    cV_wavg_fV    :: AbstractArray{Float64, 2}  # average from fine to coarse grid
    sV_wavg_cV    :: AbstractArray{Float64, 2}  # average from coase grid to slab

    sT_wsum_cT    :: AbstractArray{Float64, 2}  # average from coase grid to slab

    function AvgSpeedUpMatrix(;
        Nx             :: Int64,
        Ny             :: Int64,
        Nz_f           :: Int64,
        Nz_c           :: Int64,
        mask2          :: AbstractArray{Float64, 2},
        height_level_counts :: AbstractArray{Int64, 1},
        H_f            :: AbstractArray{Float64, 1},
        H_c            :: AbstractArray{Float64, 1},
    )

        if Nz_f < Nz_c

            throw(ErrorException("Nz_f < Nz_c"))

        end 



        num_fT = 






        elm_max = Nz_f*Nx*(Ny+1) * Nz_c 
        I = zeros(Int64, elm_max)
        J = zeros(Int64, elm_max)
        V = zeros(Float64, elm_max)
        idx = 0

        function add!(i::Int64, j::Int64, v::Float64)
            idx += 1
            I[idx] = i
            J[idx] = j
            V[idx] = v
        end

        function getSparse!(m::Int64, n::Int64)
            s = sparse(view(I, 1:idx), view(J, 1:idx), view(V, 1:idx), m, n)
            idx = 0
            return s 
        end


        wgts_f2c = zeros(Float64, Nz_f)
        _tmp = 1
        for (k_c, cnt) in enumerate(height_level_counts)
            rng = _tmp:_tmp+cnt-1
            wgts_f2c[rng] .= H_f[rng] / H_c[k_c]
            _tmp += cnt
        end

        # Making C_wavg_F
        for i=1:Nx, j=1:Ny

            if mask2[i, j] == 0.0
                continue
            end

            for k_c=1:Nz_c
                for k_f=1:Nz_f
                    i_c = flat_i(k_c, i, j, Nz_c, Nx, Ny)
                    i_f = flat_i(k_f, i, j, Nz_f, Nx, Ny)
                
                    add!(i_c, i_f, wgts_f2c[k_f])
                end
            end
        end
        C_wavg_F = getSparse!(Nz_c * Nx * Ny, Nz_f * Nx * Ny)


        H_total = sum(H_c)
        wgts_c2s = zeros(Float64, Nz_c)
        for k_c in 1:length(H_c)
            wgts_c2s[k_c] = H_c[k_c] / H_total
        end
        # Making S_wavg_C
        for i=1:Nx, j=1:Ny

            if mask2[i, j] == 0.0
                continue
            end

            for k_c=1:Nz_c
                i_c = flat_i(k_c, i, j, Nz_c, Nx, Ny)
                i_s = flat_i(  1, i, j,    1, Nx, Ny)
                
                add!(i_s, i_c, wgts_c2s[k_c])
            end
        end
        S_wavg_C = getSparse!(1 * Nx * Ny, Nz_c * Nx * Ny)

        S_wsum_C = H_total * S_wavg_C

        return new(
            C_wavg_F,
            S_wavg_C,
            S_wsum_C,
        )
    end
end
