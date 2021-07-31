mutable struct Field 

    _X_      :: AbstractArray{Float64, 4}
    _X       :: AbstractArray{Float64, 1}

    _vel      :: AbstractArray{Float64, 1}
    _u        :: AbstractArray{Float64, 1}
    _v        :: AbstractArray{Float64, 1}
    _w        :: AbstractArray{Float64, 1}


    HMXL     :: AbstractArray{Float64, 2}

    τx       :: AbstractArray{Float64, 2}
    τy       :: AbstractArray{Float64, 2}

    function Field(
        ev :: Env,
    )
        
        
        Nx, Ny, Nz = ev.gd.Nx, ev.gd.Ny, ev.gd.Nz
        
        U_pts = Nx * Ny * Nz
        V_pts = Nx * (Ny+1) * Nz
        W_pts = Nx * Ny * (Nz+1)

        _X_ = zeros(Float64, Nz, Ny, Nx, 2)
        _X  = view(_X_, :)


        _vel = zeros(Float64, U_pts + V_pts + W_pts)

        idx = 0;
        _u = view(_vel, (idx+1):(idx+U_pts)) ; idx+=U_pts
        _v = view(_vel, (idx+1):(idx+V_pts)) ; idx+=V_pts
        _w = view(_vel, (idx+1):(idx+W_pts)) ; idx+=W_pts


        HMXL = zeros(Float64, Nx, Ny)
        τx = zeros(Float64, Nx, Ny)
        τy = zeros(Float64, Nx, Ny)

        return new(

            _X_,
            _X,

            _vel,
            _u,
            _v,
            _w,


            HMXL,

            τx,
            τy,
        )

    end

end


mutable struct SugarView
    
    f :: Field

    TEMP     :: AbstractArray{Float64, 3}
    SALT     :: AbstractArray{Float64, 3}

end

