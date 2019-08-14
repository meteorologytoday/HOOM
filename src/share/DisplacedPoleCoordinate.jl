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


function getCyclicNeighbors(Nx, Ny, i, j)
    return mod(i-1 - 1, Nx) + 1, mod(i+1 - 1, Nx) + 1, mod(j-1 - 1, Ny) + 1, mod(j+1 - 1, Ny) + 1
end


function extend(
    a         :: AbstractArray{Float64, 2};
    cyclic_x  :: Bool,
    cyclic_y  :: Bool,
)
    Nx, Ny = size(a)

    aa = Array{Float64}(undef, Nx + 2, Ny + 2)
    fill!(aa, NaN)

    aa[2:end-1, 2:end-1] = a
    
    if cyclic_x
        aa[1,   2:end-1] = a[end, :]
        aa[end, 2:end-1] = a[1,   :]
    end

    if cyclic_y
        aa[2:end-1,   1] = a[:, end]
        aa[2:end-1, end] = a[:,   1]
    end

    return aa 

end


struct GridInfo

    R     :: Float64

    Nx    :: Integer
    Ny    :: Integer

    α     :: AbstractArray{Float64, 2}
    cosα  :: AbstractArray{Float64, 2}
    sinα  :: AbstractArray{Float64, 2}
    dx_w  :: AbstractArray{Float64, 2}
    dx_c  :: AbstractArray{Float64, 2}
    dx_e  :: AbstractArray{Float64, 2}
    dy_s  :: AbstractArray{Float64, 2}
    dy_c  :: AbstractArray{Float64, 2}
    dy_n  :: AbstractArray{Float64, 2}
    
    ds1   :: AbstractArray{Float64, 2}
    ds2   :: AbstractArray{Float64, 2}
    ds3   :: AbstractArray{Float64, 2}
    ds4   :: AbstractArray{Float64, 2}

    dσ    :: AbstractArray{Float64, 2}
   
    mask  :: AbstractArray{Float64, 2}        # For now any point next to sigular point should be skipped

    function GridInfo(
        R      :: Float64,
        Nx     :: Integer,
        Ny     :: Integer,
        c_lon  :: AbstractArray{Float64, 2},  # center longitude
        c_lat  :: AbstractArray{Float64, 2},  # center latitude
        vs_lon :: AbstractArray{Float64, 3},  # vertices longitude (4, Nx, Ny)
        vs_lat :: AbstractArray{Float64, 3};  # vertices latitude  (4, Nx, Ny)
        mask   :: Union{Nothing, AbstractArray{Float64, 2}} = nothing,
        angle_unit :: Symbol = :deg,
    )
   
        α    = zeros(Float64, Nx, Ny)
        cosα = zeros(Float64, Nx, Ny)    
        sinα = zeros(Float64, Nx, Ny)   
        dx_w = zeros(Float64, Nx, Ny)
        dx_c = zeros(Float64, Nx, Ny)
        dx_e = zeros(Float64, Nx, Ny)
        dy_s = zeros(Float64, Nx, Ny)
        dy_c = zeros(Float64, Nx, Ny)
        dy_n = zeros(Float64, Nx, Ny)
        dσ = zeros(Float64, Nx, Ny)
        
        ds1 = zeros(Float64, Nx, Ny)
        ds2 = zeros(Float64, Nx, Ny)
        ds3 = zeros(Float64, Nx, Ny)
        ds4 = zeros(Float64, Nx, Ny)
    
        ps = zeros(Float64, 3, 4)

        true_north  = zeros(Float64, 3)
        true_east   = zeros(Float64, 3)
        true_upward = zeros(Float64, 3)


        c_lon_rad = copy(c_lon)
        c_lat_rad = copy(c_lat)

        vs_lon_rad = copy(vs_lon)
        vs_lat_rad = copy(vs_lat)

        if angle_unit == :deg

            c_lon_rad .*= d2r
            c_lat_rad .*= d2r

            vs_lon_rad .*= d2r 
            vs_lat_rad .*= d2r
 
        elseif angle_unit == :rad
            
            # do nothing

        else
            throw(ErrorException("Unknown `angle_unit`: " * angle_unit))
        end


        # North and sounth boundaries are next to singular points
        # Here we simply skip these points
        if mask == nothing
            mask = zeros(Float64, Nx, Ny)
            mask .= 1.0
            mask[:,   1] .= 0.0
            mask[:, end] .= 0.0
        end


        for i = 1:Nx, j = 1:Ny

            for k = 1:4

                λ = vs_lon_rad[k, i, j]
                θ = vs_lat_rad[k, i, j]

                ps[1, k] = cos(θ) * cos(λ)
                ps[2, k] = cos(θ) * sin(λ)
                ps[3, k] = sin(θ)

            end

            ps .*= R

            u1 = ps[:, 2] - ps[:, 1]
            u2 = ps[:, 3] - ps[:, 2]
            u3 = ps[:, 3] - ps[:, 4]
            u4 = ps[:, 4] - ps[:, 1]

            ds1[i, j] = norm(u1)
            ds2[i, j] = norm(u2)
            ds3[i, j] = norm(u3)
            ds4[i, j] = norm(u4)

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

            cosα[i, j] = cos(α[i, j])
            sinα[i, j] = sin(α[i, j])

        end
        
        for i = 1:Nx, j = 1:Ny

            if mask[i, j] == 0.0 
                continue
            end

            i_w, i_e, j_s, j_n = getCyclicNeighbors(Nx, Ny, i, j)

            dx_w[i, j] = (dx_c[i_w, j] + dx_c[i,   j]) / 2.0
            dx_e[i, j] = (dx_c[i,   j] + dx_c[i_e, j]) / 2.0

            dy_s[i, j] = (dy_c[i, j_s] + dy_c[i, j  ]) / 2.0
            dy_n[i, j] = (dy_c[i, j  ] + dy_c[i, j_n]) / 2.0

        end


        return new(
            R, Nx, Ny,
            α, cosα, sinα,
            dx_w, dx_c, dx_e,
            dy_s, dy_c, dy_n,
            ds1, ds2, ds3, ds4,
            dσ,
            mask,
        )
 
    end
end


function project!(
    gi    :: GridInfo,
    ivf_e :: AbstractArray{Float64, 2},     # input vector field east
    ivf_n :: AbstractArray{Float64, 2},     # input vector field north
    ovf_e :: AbstractArray{Float64, 2},    # output vector field east
    ovf_n :: AbstractArray{Float64, 2};    # output vector field north
    direction = :Forward,
)

    if direction == :Forward   # from outside world onto dispalced pole grid

        for i=1:gi.Nx, j=1:gi.Ny
            ovf_e[i, j] =   ivf_e[i, j] * gi.cosα[i, j] - ivf_n[i, j] * gi.sinα[i, j]
            ovf_n[i, j] =   ivf_e[i, j] * gi.sinα[i, j] + ivf_n[i, j] * gi.cosα[i, j]
        end
 
    else                       # from displaced pole grid onto outside world

        for i=1:gi.Nx, j=1:gi.Ny
            ovf_e[i, j] =   ivf_e[i, j] * gi.cosα[i, j] + ivf_n[i, j] * gi.sinα[i, j]
            ovf_n[i, j] = - ivf_e[i, j] * gi.sinα[i, j] + ivf_n[i, j] * gi.cosα[i, j]
        end

    end

end


function DIV!(
    gi   :: GridInfo,
    vf_e :: AbstractArray{Float64, 2},
    vf_n :: AbstractArray{Float64, 2},
    div  :: AbstractArray{Float64, 2},
)

    for i=1:gi.Nx, j=1:gi.Ny

        #if gi.mask[i, j] == 0.0
        #    div[i, j] = NaN
        #    continue
        #end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)

        flux_e = (gi.mask[i_e, j  ] == 0.0) ? 0.0 : vf_e[i_e, j] * gi.ds2[i, j]
        flux_w = (gi.mask[i_w, j  ] == 0.0) ? 0.0 : vf_e[i_w, j] * gi.ds4[i, j] 
        flux_n = (gi.mask[i  , j_n] == 0.0) ? 0.0 : vf_n[i, j_n] * gi.ds3[i, j]
        flux_s = (gi.mask[i  , j_s] == 0.0) ? 0.0 : vf_n[i, j_s] * gi.ds1[i, j]

        div[i, j] =  (  flux_e - flux_w + flux_n - flux_s ) / gi.dσ[i, j]

    end
        
end



function ∇∇!(
    gi   :: GridInfo,
    f    :: AbstractArray{Float64, 2},
    ∇∇f  :: AbstractArray{Float64, 2},
)


    for i=1:gi.Nx, j=1:gi.Ny

        if gi.mask[i, j] == 0.0
            ∇∇f[i, j] = NaN
            continue
        end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)

        f_c = f[i, j]

        flux_e = (gi.mask[i_e, j  ] == 0.0) ? 0.0 : (f[i_e, j] - f_c) / gi.dx_e[i, j]
        flux_w = (gi.mask[i_w, j  ] == 0.0) ? 0.0 : (f_c - f[i_w, j]) / gi.dx_w[i, j] 
        flux_n = (gi.mask[i  , j_n] == 0.0) ? 0.0 : (f[i, j_n] - f_c) / gi.dy_n[i, j]
        flux_s = (gi.mask[i  , j_s] == 0.0) ? 0.0 : (f_c - f[i, j_s]) / gi.dy_s[i, j]

        ∇∇f[i, j] = ( flux_e - flux_w ) / gi.dx_c[i, j] + ( flux_n - flux_s ) / gi.dy_c[i, j]

    end
        
end


function hadv_upwind!(
    gi    :: GridInfo,
    hadvs :: AbstractArray{Float64, 2},
    us    :: AbstractArray{Float64, 2},
    vs    :: AbstractArray{Float64, 2},
    qs    :: AbstractArray{Float64, 2},
)

    for i=1:gi.Nx, j=1:gi.Ny

        if gi.mask[i, j] == 0.0
            hadvs[i, j] = NaN
            continue
        end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)
        
        hadv = 0.0
        u = us[i, j]
        v = vs[i, j]

        if u > 0 && gi.mask[i_w, j] != 0.0           # use point on the west
            hadv += - u * ( qs[i, j] - qs[i_w, j] ) / gi.dx_w[i, j]
        elseif u <= 0 && gi.mask[i_e, j] != 0.0      # use point on the east
            hadv += - u * ( qs[i_e, j] - qs[i, j] ) / gi.dx_e[i, j]
        end

        if v > 0 && gi.mask[i, j_s] != 0.0           # use point on the south
            hadv += - v * ( qs[i, j] - qs[i, j_s] ) / gi.dy_s[i, j]
        elseif v <= 0 && gi.mask[i, j_n] != 0.0      # use point on the north
            hadv += - v * ( qs[i, j_n] - qs[i, j] ) / gi.dy_n[i, j]
        end

        hadvs[i, j] = hadv
    end

end

end
