

module DisplacedPoleCoordinate


using LinearAlgebra: ⋅, normalize, normalize!, norm

"""

This function calculate the angle between the real north and the
north defined by vertices. I call it T-north.

Each T-grid is defined by four vertices. As shown below:

[4] --- [3]
 |       | 
 |       | 
[1] --- [2]

T-north for this grid is defined by the average of [1]->[4] and
[2]->[3]. Likewise, T-east is the east defined by vertices, the
average of [1]->[2], [4]->[3]

Positive angle α means the local T-north and T-east will align
with true north and east if rotating true north counterclockwise
by angle α. 

"""

d2r = π / 180.0
r2d = 180 / π


mutable struct GridInfo
    Nx    :: Integer
    Ny    :: Integer
    α     :: AbstractArray{Float64, 2}
    dx_w  :: AbstractArray{Float64, 2}
    dx_c  :: AbstractArray{Float64, 2}
    dx_e  :: AbstractArray{Float64, 2}
    dy_s  :: AbstractArray{Float64, 2}
    dy_c  :: AbstractArray{Float64, 2}
    dy_n  :: AbstractArray{Float64, 2}
    dσ    :: AbstractArray{Float64, 2}
   

    function GridInfo(
        Nx     :: Integer,
        Ny     :: Integer,
        c_lon  :: AbstractArray{Float64, 2},  # center longitude
        c_lat  :: AbstractArray{Float64, 2},  # center latitude
        vs_lon :: AbstractArray{Float64, 3},  # vertices longitude (4, Nx, Ny)
        vs_lat :: AbstractArray{Float64, 3},  # vertices latitude  (4, Nx, Ny)
    )


    
        α  = zeros(Float64, Nx, Ny)
        dx_w = zeros(Float64, Nx, Ny)
        dx_c = zeros(Float64, Nx, Ny)
        dx_e = zeros(Float64, Nx, Ny)
        dy_s = zeros(Float64, Nx, Ny)
        dy_c = zeros(Float64, Nx, Ny)
        dy_n = zeros(Float64, Nx, Ny)
        dσ = zeros(Float64, Nx, Ny)
    
        ps = zeros(Float64, 4, 3)

        true_north  = zeros(Float64, 3)
        true_east   = zeros(Float64, 3)
        true_upward = zeros(Float64, 3)


        c_lon_rad = d2r * c_lon
        c_lat_rad = d2r * c_lat

        vs_lon_rad = d2r * vs_lon
        vs_lat_rad = d2r * vs_lat

        for i = 1:Nx, j = 1:Ny

            for k = 1:4

                λ = vs_lon_rad[k, i, j]
                θ = vs_lat_rad[k, i, j]

                ps[k, 1] = cos(θ) * cos(λ)
                ps[k, 2] = cos(θ) * sin(λ)
                ps[k, 3] = sin(θ)

            end

            u1 = ps[2, :] - ps[1, :]
            u2 = ps[3, :] - ps[2, :]
            u3 = ps[3, :] - ps[4, :]
            u4 = ps[4, :] - ps[1, :]


            dx_c[i, j] = (norm(u1) + norm(u3)) / 2.0
            dy_c[i, j] = (norm(u2) + norm(u4)) / 2.0
            dσ[i, j] = dx_c[i, j] * dy_c[i, j]

            grid_east  = u1+u3
            grid_north = u2+u4

            λc = c_lon_rad[i, j]
            θc = c_lat_rad[i, j]

            true_north[:]  = [ - sin(θc) * cos(λc), - sin(θc) * sin(λc),   cos(θc) ]
            true_east[:]   = [ - sin(λc)          ,   cos(λc)          ,   0.0     ]
            true_upward[:] = [   cos(θc) * cos(λc),   cos(θc) * sin(λc),   sin(θc) ]

            grid_north -= (grid_north ⋅ true_upward) * true_upward

            cos_α = (grid_north ⋅ true_north) / norm(grid_north) / norm(true_north)
            cos_β = (grid_north ⋅ true_east)  / norm(grid_north)

            if abs(cos_α) > 1.0 && abs(cos_α - 1 < 1e-3)
                cos_α = (cos_α > 0) ? 1.0 : -1.0
            end

            α[i, j] = acos(cos_α)

            if cos_β < 0.0   # Flip angle
                α[i, j] = 2*π - α[i, j]
            end

        end
        
        for i = 1:Nx, j = 1:Ny

            i_w = i-1
            i_e = i+1
            j_s = j-1
            j_n = j+1

            if i_w == 0
                i_w = Nx
            end

            if i_e == Nx+1
                i_e = 1
            end

            if j_s == 0
                continue
            end

            if j_n == Ny+1
                continue
            end

            dx_w[i, j] = (dx_c[i_w, j] + dx_c[i,   j]) / 2.0
            dx_e[i, j] = (dx_c[i,   j] + dx_c[i_e, j]) / 2.0

            dy_s[i, j] = (dy_c[i, j_s] + dy_c[i, j  ]) / 2.0
            dy_n[i, j] = (dy_c[i, j  ] + dy_c[i, j_n]) / 2.0

        end

        α .*= r2d

        return new(
            Nx, Ny,
            α,
            dx_w, dx_c, dx_e,
            dy_s, dy_c, dy_n,
            dσ,
        )
    
 
    end
end


end
