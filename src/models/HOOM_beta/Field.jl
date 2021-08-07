mutable struct Field 

    _X_      :: AbstractArray{Float64, 2}
    _X       :: AbstractArray{Float64, 1}

    _vel      :: AbstractArray{Float64, 1}
    _u        :: AbstractArray{Float64, 1}
    _v        :: AbstractArray{Float64, 1}
    _w        :: AbstractArray{Float64, 1}

    _Xflx_U_  :: AbstractArray{Float64, 2}
    _Xflx_V_  :: AbstractArray{Float64, 2}
    _Xflx_W_  :: AbstractArray{Float64, 2}

    _ADVX_    :: AbstractArray{Float64, 2}


    _CHKX_   :: AbstractArray{Float64, 2}
    _CHKX    :: AbstractArray{Float64, 1}

    _TMP_CHKX_   :: AbstractArray{Float64, 2} # Used to store the ∫ X dz before steping
    _TMP_CHKX    :: AbstractArray{Float64, 1}

    _TMP_SUBSTEP_BUDGET_   :: AbstractArray{Float64, 2} # Used to store substep tracer budget


    HMXL     :: AbstractArray{Float64, 3}

    SWFLX    :: AbstractArray{Float64, 3}
    NSWFLX   :: AbstractArray{Float64, 3}
 
    TAUX         :: AbstractArray{Float64, 3}
    TAUY         :: AbstractArray{Float64, 3}
   
    TAUX_east    :: AbstractArray{Float64, 3}
    TAUY_north   :: AbstractArray{Float64, 3}

    QFLX2ATM     :: AbstractArray{Float64, 3}


    # sugar view
    sv       :: Union{Nothing, Dict} 


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

        sT_pts = Nx * Ny * 1

        _CHKX_ = zeros(Float64, sT_pts, 2)
        _CHKX  = view(_CHKX_, :)
 
        _TMP_CHKX_ = zeros(Float64, sT_pts, 2)
        _TMP_CHKX  = view(_TMP_CHKX_, :)
    
        _TMP_SUBSTEP_BUDGET_ = zeros(Float64, sT_pts, 2)
        
        _vel = zeros(Float64, U_pts + V_pts + W_pts)

        idx = 0;
        _u = view(_vel, (idx+1):(idx+U_pts)) ; idx+=U_pts
        _v = view(_vel, (idx+1):(idx+V_pts)) ; idx+=V_pts
        _w = view(_vel, (idx+1):(idx+W_pts)) ; idx+=W_pts

        _Xflx_U_ = zeros(Float64, U_pts, 2)
        _Xflx_V_ = zeros(Float64, V_pts, 2)
        _Xflx_W_ = zeros(Float64, W_pts, 2)
        
        _ADVX_ = zeros(Float64, T_pts, 2)

        HMXL = zeros(Float64, 1, Nx, Ny)
        SWFLX = zeros(Float64, 1, Nx, Ny)
        NSWFLX = zeros(Float64, 1, Nx, Ny)
        TAUX_east = zeros(Float64, 1, Nx, Ny)
        TAUY_north = zeros(Float64, 1, Nx, Ny)
 
        TAUX = zeros(Float64, 1, Nx, Ny)
        TAUY = zeros(Float64, 1, Nx, Ny)
        
        QFLX2ATM = zeros(Float64, 1, Nx, Ny)
        
        sv = Dict(
            :TEMP => reshape(view(_X_, :, 1), Nz, Nx, Ny),
            :SALT => reshape(view(_X_, :, 2), Nz, Nx, Ny),
            :UVEL => reshape(_u, Nz, Nx, Ny),
            :VVEL => reshape(_v, Nz, Nx, Ny+1),
            :WVEL => reshape(_w, Nz+1, Nx, Ny),
            :CHKTEMP => reshape(view(_CHKX_, :, 1), 1, Nx, Ny),
            :CHKSALT => reshape(view(_CHKX_, :, 2), 1, Nx, Ny),
            :ADVT => reshape(view(_ADVX_, :, 1), Nz, Nx, Ny),
        )

        fi = new(

            _X_,
            _X,

            _vel,
            _u,
            _v,
            _w,

            _Xflx_U_,
            _Xflx_V_,
            _Xflx_W_,

            _ADVX_,

            _CHKX_,
            _CHKX,

            _TMP_CHKX_,
            _TMP_CHKX,

            _TMP_SUBSTEP_BUDGET_,
            HMXL,

            SWFLX,
            NSWFLX,

            TAUX,
            TAUY,

            TAUX_east,
            TAUY_north,

            QFLX2ATM,

            nothing,
        )

        fi.sv = getSugarView(ev, fi)

        return fi

    end

end

function getSugarView(
    ev :: Env,
    fi :: Field,
)

    Nx, Ny, Nz = ev.Nx, ev.Ny, ev.Nz

    return sv = Dict(
        :TEMP => reshape(view(fi._X_, :, 1), Nz, Nx, Ny),
        :SALT => reshape(view(fi._X_, :, 2), Nz, Nx, Ny),
        :UVEL => reshape(fi._u, Nz, Nx, Ny),
        :VVEL => reshape(fi._v, Nz, Nx, Ny+1),
        :WVEL => reshape(fi._w, Nz+1, Nx, Ny),
        :CHKTEMP => reshape(view(fi._CHKX_, :, 1), 1, Nx, Ny),
        :CHKSALT => reshape(view(fi._CHKX_, :, 2), 1, Nx, Ny),
        :ADVT => reshape(view(fi._ADVX_, :, 1), Nz, Nx, Ny),
    )

end
