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

    C_avg_F    :: AbstractArray{Float64, 2}  # average from fine to coarse grid
    S_avg_C    :: AbstractArray{Float64, 2}  # average from coase grid to slab

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

        elm_max = Nz_f*Nx*(Ny+1) * 2 
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
        for (k_c, cnt) in height_level_counts
            rng = _tmp:_tmp+cnt-1
            wgts_f2c[rng] .= H_f[rng] / H_c[k_c]
            _tmp += cnt
        end

        # Making C_avg_F
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
        C_avg_F = getSparse!(Nz_c * Nx * Ny, Nz_f * Nx * Ny)


        H_total = sum(H_c)
        wgts_c2s = zeros(Float64, Nz_c)
        _tmp = 1
        for k_c=1:Nz_c
            wgts_c2s[rng] = H_c[k_c] / H_total
        end
        # Making S_avg_C
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
        S_avg_C = getSparse!(1 * Nx * Ny, Nz_c * Nx * Ny)



        return new(
            C_avg_F,
            S_avg_C,
        )
    end
end
