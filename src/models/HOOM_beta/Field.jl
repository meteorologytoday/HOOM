mutable struct Field 

    _X_      :: AbstractArray{Float64, 2}
    _X       :: AbstractArray{Float64, 1}

    _vel      :: AbstractArray{Float64, 1}
    _u        :: AbstractArray{Float64, 1}
    _v        :: AbstractArray{Float64, 1}
    _w        :: AbstractArray{Float64, 1}


    HMXL     :: AbstractArray{Float64, 3}

    τx       :: AbstractArray{Float64, 3}
    τy       :: AbstractArray{Float64, 3}

    # sugar view
    sv       :: Dict 

    function Field(
        ev :: Env,
    )

        Nz, Nx, Ny = ev.Nz, ev.Nx, ev.Ny
        
        T_pts = Nx * Ny * Nz
        U_pts = Nx * Ny * Nz
        V_pts = Nx * (Ny+1) * Nz
        W_pts = Nx * Ny * (Nz+1)

        _X_ = zeros(Float64, T_pts, 2)
        _X  = view(_X_, :)


        _vel = zeros(Float64, U_pts + V_pts + W_pts)

        idx = 0;
        _u = view(_vel, (idx+1):(idx+U_pts)) ; idx+=U_pts
        _v = view(_vel, (idx+1):(idx+V_pts)) ; idx+=V_pts
        _w = view(_vel, (idx+1):(idx+W_pts)) ; idx+=W_pts


        HMXL = zeros(Float64, 1, Nx, Ny)
        τx = zeros(Float64, 1, Nx, Ny)
        τy = zeros(Float64, 1, Nx, Ny)
        println("shape of _X_: ", size(_X_))
        println("shape : ", Nz, ", ", Nx, "," , Ny)
        sv = Dict(
            :TEMP => reshape(view(_X_, :, 1), Nz, Nx, Ny),
            :SALT => reshape(view(_X_, :, 2), Nz, Nx, Ny),
            :UVEL => reshape(_u, Nz, Nx, Ny),
            :VVEL => reshape(_v, Nz, Nx, Ny+1),
            :WVEL => reshape(_w, Nz+1, Nx, Ny),
        )


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

            sv,
        )

    end

end


mutable struct SugarView
    
    f :: Field

    TEMP     :: AbstractArray{Float64, 3}
    SALT     :: AbstractArray{Float64, 3}

end

