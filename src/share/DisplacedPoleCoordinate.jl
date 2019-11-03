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

    c_lon :: AbstractArray{Float64, 2}
    c_lat :: AbstractArray{Float64, 2}

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

    weight_e :: AbstractArray{Float64, 2}    # ( Nx+1 , Ny   )
    weight_n :: AbstractArray{Float64, 2}    # ( Nx   , Ny+1 )
    DX    :: AbstractArray{Float64, 2}       # ( Nx   , Ny+1 )   The X-side grid size.
    DY    :: AbstractArray{Float64, 2}       # ( Nx+1 , Ny   )   The Y-side grid size.

  
    function GridInfo(
        R      :: Float64,
        Nx     :: Integer,
        Ny     :: Integer,
        c_lon  :: AbstractArray{Float64, 2},  # center longitude
        c_lat  :: AbstractArray{Float64, 2},  # center latitude
        vs_lon :: AbstractArray{Float64, 3},  # vertices longitude (4, Nx, Ny)
        vs_lat :: AbstractArray{Float64, 3},  # vertices latitude  (4, Nx, Ny)
        area   :: AbstractArray{Float64, 2};  # area in radian^2
        angle_unit :: Symbol = :deg,
        sub_yrng :: Union{Colon, UnitRange} = Colon(),
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
        
        weight_e = zeros(Float64, Nx+1, Ny)
        weight_n = zeros(Float64, Nx, Ny+1)

        DX = zeros(Float64, Nx, Ny+1)
        DY = zeros(Float64, Nx+1, Ny)



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
            dσ[i, j] = R^2.0 * area[i, j]  #dx_c[i, j] * dy_c[i, j]

            grid_east  = u1+u3
            grid_north = u2+u4

            λc = c_lon_rad[i, j]
            θc = c_lat_rad[i, j]

            true_north[:]  = [ - sin(θc) * cos(λc), - sin(θc) * sin(λc),   cos(θc) ]
            true_east[:]   = [ - sin(λc)          ,   cos(λc)          ,   0.0     ]
            true_upward[:] = [   cos(θc) * cos(λc),   cos(θc) * sin(λc),   sin(θc) ]

            #grid_north -= (grid_north ⋅ true_upward) * true_upward

            cos_α = (grid_north ⋅ true_north) / norm(grid_north) #/ norm(true_north)
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

            i_w, i_e, j_s, j_n = getCyclicNeighbors(Nx, Ny, i, j)

            dx_w[i, j] = (dx_c[i_w, j] + dx_c[i,   j]) / 2.0
            dx_e[i, j] = (dx_c[i,   j] + dx_c[i_e, j]) / 2.0

            dy_s[i, j] = (dy_c[i, j_s] + dy_c[i, j  ]) / 2.0
            dy_n[i, j] = (dy_c[i, j  ] + dy_c[i, j_n]) / 2.0

        end

        #
        # weight_e[i, ?] is the relative portion of grid weighting
        # to the west of the boundary of eastward vectors
        #
        #                     bnd
        #                     [i]
        #     |                |                  |
        #     |    grid[i-1]   |    grid[i]       |
        #     |                |                  |
        #     | <-- dx[i-1]--> | <--- dx[i] --->  |
        #     |                |                  |
        #
        #   weight_e[i, ?] = dx[i-1] / (dx[i-1] + dx[i])
        #
        #   weight_n would be the same idea with 
        #   portion represent the south grid
        #

        for j = 1:Ny
            weight_e[1, j] = weight_e[Nx+1, j] = dx_c[Nx, j] / (dx_c[Nx, j] + dx_c[1, j])
        end
        for i = 2:Nx, j = 1:Ny
            weight_e[i, j] = dx_c[i-1, j] / (dx_c[i-1, j] + dx_c[i, j])
        end

        # Ignore the northest and the southest because information
        # is unknown
        for i = 1:Nx, j = 2:Ny
            weight_n[i, j] = dy_s[i, j-1] / ( dy_c[i, j-1] + dy_c[i, j] )
        end


        # Calculate DX
        for i = 1:Nx, j = 1:Ny
            DX[i, j] = ds1[i, j]
        end
        DX[:, Ny+1] = ds3[:, Ny]

        # Calculate DY
        for i = 1:Nx, j = 1:Ny
            DY[i, j] = ds4[i, j]
        end
        DY[Nx+1, :] = DY[1, :]

        if sub_yrng == Colon()
            sub_yrng = 1:Ny
        end

        new_Ny = length(sub_yrng)
        sub_yrng_ext = sub_yrng[1]:sub_yrng[end]+1
       
        return new(
            R,
            Nx,
            new_Ny,
            c_lon_rad[:, sub_yrng],
            c_lat_rad[:, sub_yrng],
            α[:, sub_yrng],
            cosα[:, sub_yrng],
            sinα[:, sub_yrng],
            dx_w[:, sub_yrng],
            dx_c[:, sub_yrng],
            dx_e[:, sub_yrng],
            dy_s[:, sub_yrng],
            dy_c[:, sub_yrng],
            dy_n[:, sub_yrng],
            ds1[:, sub_yrng],
            ds2[:, sub_yrng],
            ds3[:, sub_yrng],
            ds4[:, sub_yrng],
            dσ[:, sub_yrng],
            weight_e[:, sub_yrng],
            weight_n[:, sub_yrng_ext],
            DX[:, sub_yrng_ext],
            DY[:, sub_yrng],
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

function GRAD!(
    gi     :: GridInfo,
    scalar :: AbstractArray{Float64, 2},
    vf_e   :: AbstractArray{Float64, 2},
    vf_n   :: AbstractArray{Float64, 2},
    mask   :: AbstractArray{Float64, 2},
)

    for i=1:gi.Nx, j=1:gi.Ny

        if mask[i, j] == 0.0
            vf_e[i, j] = NaN
            vf_n[i, j] = NaN
            continue
        end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)

        s_c = scalar[i, j]
        
        if mask[i_e, j] != 0
            if mask[i_w, j] != 0
                # i_e = 1, i_w = 1 (Both side derivative)
                vf_e[i, j] = ( scalar[i_e, j] - scalar[i_w, j] ) / (2.0 * gi.dx_c[i, j])
            else
                # i_e = 1, i_w = 0 (East side derivative)
                vf_e[i, j] = ( scalar[i_e, j] - s_c ) / gi.dx_c[i, j]
            end
        elseif mask[i_w, j] != 0
                # i_e = 0, i_w = 1 (West side derivative)
                vf_e[i, j] = ( s_c - scalar[i_w, j] ) / gi.dx_c[i, j]
        else
                # i_e = 0, i_w = 0 (0 derivative)
                vf_e[i, j] = 0.0
        end

        if mask[i, j_n] != 0
            if mask[i, j_s] != 0
                # j_n = 1, j_s = 1 (Both side derivative)
                vf_n[i, j] = ( scalar[i, j_n] - scalar[i, j_s] ) / (2.0 * gi.dy_c[i, j])
            else
                # j_n = 1, j_s = 0 (North side derivative)
                vf_n[i, j] = ( scalar[i, j_n] - s_c ) / gi.dy_c[i, j]
            end
        elseif mask[i, j_s] != 0
                # j_n = 0, j_s = 1 (South side derivative)
                vf_n[i, j] = ( s_c - scalar[i, j_s] ) / gi.dy_c[i, j]
        else
                # j_n = 0, j_s = 0 (0 derivative)
                vf_n[i, j] = 0.0
        end

    end
        
end





function DIV!(
    gi   :: GridInfo,
    vf_e :: AbstractArray{Float64, 2},
    vf_n :: AbstractArray{Float64, 2},
    div  :: AbstractArray{Float64, 2},
    mask :: AbstractArray{Float64, 2},
)

    for i=1:gi.Nx, j=1:gi.Ny

        #if gi.mask[i, j] == 0.0
        #    div[i, j] = NaN
        #    continue
        #end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)

        vf_e_c = vf_e[i, j]
        vf_n_c = vf_n[i, j]

        flux_e = (mask[i_e, j  ] == 0.0) ? 0.0 : ( vf_e[i_e, j] + vf_e_c ) / 2.0 * gi.ds2[i, j]
        flux_w = (mask[i_w, j  ] == 0.0) ? 0.0 : ( vf_e[i_w, j] + vf_e_c ) / 2.0 * gi.ds4[i, j] 
        flux_n = (mask[i  , j_n] == 0.0) ? 0.0 : ( vf_n[i, j_n] + vf_n_c ) / 2.0 * gi.ds3[i, j]
        flux_s = (mask[i  , j_s] == 0.0) ? 0.0 : ( vf_n[i, j_s] + vf_n_c ) / 2.0 * gi.ds1[i, j]

        div[i, j] =  (  flux_e - flux_w + flux_n - flux_s ) / gi.dσ[i, j]

    end
        
end





function ∇∇!(
    gi   :: GridInfo,
    f    :: AbstractArray{Float64, 2},
    ∇∇f  :: AbstractArray{Float64, 2},
    mask :: AbstractArray{Float64, 2},
)


    for i=1:gi.Nx, j=1:gi.Ny

        if mask[i, j] == 0.0
            ∇∇f[i, j] = NaN
            continue
        end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)

        f_c = f[i, j]

        flux_e = (mask[i_e, j  ] == 0.0) ? 0.0 : (f[i_e, j] - f_c) / gi.dx_e[i, j]
        flux_w = (mask[i_w, j  ] == 0.0) ? 0.0 : (f_c - f[i_w, j]) / gi.dx_w[i, j] 
        flux_n = (mask[i  , j_n] == 0.0) ? 0.0 : (f[i, j_n] - f_c) / gi.dy_n[i, j]
        flux_s = (mask[i  , j_s] == 0.0) ? 0.0 : (f_c - f[i, j_s]) / gi.dy_s[i, j]

        ∇∇f[i, j] = ( flux_e - flux_w ) / gi.dx_c[i, j] + ( flux_n - flux_s ) / gi.dy_c[i, j]

    end
        
end

function hadv_upwind!(
    gi    :: GridInfo,
    hadvs :: AbstractArray{Float64, 2},
    us    :: AbstractArray{Float64, 2},
    vs    :: AbstractArray{Float64, 2},
    qs    :: AbstractArray{Float64, 2},
    mask  :: AbstractArray{Float64, 2},
)

    for i=1:gi.Nx, j=1:gi.Ny

        if mask[i, j] == 0.0
            hadvs[i, j] = NaN
            continue
        end

        i_w, i_e, j_s, j_n = getCyclicNeighbors(gi.Nx, gi.Ny, i, j)
        
        hadv = 0.0
        u = us[i, j]
        v = vs[i, j]

        if u > 0 && mask[i_w, j] != 0.0           # use point on the west
            hadv += - u * ( qs[i, j] - qs[i_w, j] ) / gi.dx_w[i, j]

        #    if qs[i, j] != qs[i_w, j]
        #        println("[1] i,j = ", i, ", ", j)
        #        println("qs[  i,j] = ", qs[i, j])
        #        println("qs[i_w,j] = ", qs[i_w, j])
        #        throw(ErrorException("!!!"))
        #    end

        elseif u <= 0 && mask[i_e, j] != 0.0      # use point on the east
            hadv += - u * ( qs[i_e, j] - qs[i, j] ) / gi.dx_e[i, j]

        #    if qs[i, j] != qs[i_e, j]
        #        println("[2] i,j = ", i, ", ", j)
        #        throw(ErrorException("!!!"))
        #    end


        end

        if v > 0 && mask[i, j_s] != 0.0           # use point on the south
            hadv += - v * ( qs[i, j] - qs[i, j_s] ) / gi.dy_s[i, j]
        elseif v <= 0 && mask[i, j_n] != 0.0      # use point on the north
            hadv += - v * ( qs[i, j_n] - qs[i, j] ) / gi.dy_n[i, j]
        end

        hadvs[i, j] = hadv
    end

end


end
