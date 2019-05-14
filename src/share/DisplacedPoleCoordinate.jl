

module DisplacedPoleCoordinate


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

mutable struct GridInfo
    Nx    :: Integer
    Ny    :: Integer
    α     :: AbstractArray{Float64, 2}
    dx    :: AbstractArray{Float64, 2}
    dy    :: AbstractArray{Float64, 2}
    dσ     :: AbstractArray{Float64, 2}
   

    function GridInfo(
        Nx     :: Integer,
        Ny     :: Integer,
        c_lon  :: AbstractArray{Float64, 2},  # center longitude
        c_lat  :: AbstractArray{Float64, 2},  # center latitude
        vs_lon :: AbstractArray{Float64, 3},  # vertices longitude
        vs_lat :: AbstractArray{Float64, 3},  # vertices latitude
    )

        d2r = π / 180.0
    
        α  = zeros(Float64, Nx, Ny)
        dx = zeros(Float64, Nx, Ny)
        dy = zeros(Float64, Nx, Ny)
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

                λ = vs_lon_rad[i, j, k]
                θ = vs_lat_rad[i, j, k]

                ps[k, 1] = cos(θ) * cos(λ)
                ps[k, 2] = cos(θ) * sin(λ)
                ps[k, 3] = sin(θ)

            end

            u1 = ps[2, :] - ps[1, :]
            u2 = ps[3, :] - ps[2, :]
            u3 = ps[3, :] - ps[4, :]
            u4 = ps[4, :] - ps[1, :]

            deast  = u1+u3
            dnorth = u2+u4

            λc = c_lon_rad
            θc = c_lat_rad

            true_north[:]  = [ - sin(θc) * cos(λc), - sin(θc) * sin(λc),   cos(θc) ]
            true_east[:]   = [ - sin(λc)          ,   cos(λc)          ,   0.0     ]
            true_upward[:] = [   cos(θc) * cos(λc),   cos(θc) * sin(λc),   sin(θc) ]
            
            
        end
        

        return new(Nx, Ny, α, dx, dy, dσ)
    
 
    end
end


function calGridInformation!(
    Nx     :: Integer,
    Ny     :: Integer,
    vs_lon :: AbstractArray{Float64, 3},
    vs_lat :: AbstractArray{Float64, 3},
)


    # Calculate angle between true north and grid north
    for i = 1:Nx, j = 1:Ny
        
        
    end

end


function calDisplacedAngle(
    nx :: Integer,
    ny :: Integer,
    
)







end
