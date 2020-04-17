mutable struct DynamicAdv    # Adam Bashford

    G_bt :: AbstractArray{Float64, 4}  # G term for barotropic mode.   (Nx, Ny,       UV, time)
    G_bc :: AbstractArray{Float64, 5}  # G term for baroclinic mode.   (Nx, Ny, Nz_c, UV, time)

    GRAD_b_f :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nx, Ny, Nz_f)
    GRAD_b   :: AbstractArray{Float64, 3}   # Horizontal gradient of buoyancy  (Nx, Ny, Nz_c)
    
    Cof      :: AbstractArray{Float64, 4}   # Coriolis force                   (Nx, Ny, Nz_c, UV)
    V_GRAD   :: AbstractArray{Float64, 3}   # Horizontal gradient of velocity  (Nx, Ny, Nz_c, UV)
    V_ADV    :: AbstractArray{Float64, 4}   # Advection of velocity            (Nx, Ny, Nz_c, UV)

    

    
    function DynamicAdv(env, state, datakind)
        



        return new(
            
        )
    end
end


